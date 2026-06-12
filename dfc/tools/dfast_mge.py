#! /usr/bin/env python
# coding: UTF8

import os
import json
from Bio.SeqFeature import (FeatureLocation, ExactPosition,
                            BeforePosition, AfterPosition)
from ..models.bio_feature import ExtendedFeature
from .base_tools import ContigAnnotationTool

# JSON short type code -> INSDC /mobile_element_type value（mobile_element として登録）
_MOBILE_ELEMENT_TYPE = {
    "is": "insertion sequence",
    "mite": "MITE",
    "tn": "transposon",            # unit transposon
    "cn": "transposon",            # composite transposon (DB hit)
    "integron": "integron",
    "retrotransposon": "retrotransposon",
}

# PredictionEvidence.PUTATIVE の JSON 値（io.py の `r["evidence"] == 2` に対応）
PUTATIVE_EVIDENCE = 2

# JSON short type code -> 読みやすい正式名（report/amr_summary 用）
_TYPE_NAME = {
    "is": "insertion sequence",
    "mite": "MITE",
    "tn": "unit transposon",
    "cn": "composite transposon",
    "ice": "integrative conjugative element",
    "aice": "actinomycete integrative conjugative element",
    "ime": "integrative mobilizable element",
    "cime": "cis-mobilizable element",
    "mic": "mobile insertion cassette",
    "iscr": "IS common region",
    "integron": "integron",
    "retrotransposon": "retrotransposon",
    "other": "mobile genetic element",
}


def classify_mge(type_code, evidence):
    """MEF の type 短縮コードと evidence から (feature_key, mobile_element_type, is_putative) を返す。

    - putative composite transposon (type=="cn" かつ evidence==2) は misc_feature。
    - INSDC /mobile_element_type に明確対応する型は mobile_element。
    - それ以外（ICE/IME/CIME/MIC 等・未知）は misc_feature。
    """
    is_putative = (type_code == "cn" and evidence == PUTATIVE_EVIDENCE)
    if is_putative:
        return ("misc_feature", None, True)
    met = _MOBILE_ELEMENT_TYPE.get(type_code)
    if met is not None:
        return ("mobile_element", met, False)
    return ("misc_feature", None, False)


# misc_feature として登録する型の説明ラベル
_MISC_LABEL = {
    "ice": "integrative conjugative element (ICE)",
    "aice": "actinomycete integrative conjugative element (AICE)",
    "ime": "integrative mobilizable element (IME)",
    "cime": "cis-mobilizable element (CIME)",
    "mic": "mobile insertion cassette",
    "iscr": "IS common region (ISCR)",
    "other": "mobile genetic element",
}


def build_qualifiers(entry):
    """MEF JSON の result エントリから GenBank qualifiers(dict) を作る。"""
    name = entry["name"]
    type_code = entry["type"]
    evidence = entry.get("evidence", 1)
    feature_key, met, is_putative = classify_mge(type_code, evidence)

    identity = float(entry.get("identity", 0)) * 100
    coverage = float(entry.get("coverage", 0)) * 100
    accession = (entry.get("template") or {}).get("accession", "")

    detail = "MobileElementFinder: {name}; identity:{id:.1f}%, coverage:{cov:.1f}%".format(
        name=name, id=identity, cov=coverage)
    if accession:
        detail += "; similar to {acc} (MGEdb)".format(acc=accession)
    notes = [detail]
    if is_putative:
        notes.append("putative composite transposon (predicted by MobileElementFinder)")
    elif feature_key == "misc_feature":
        label = _MISC_LABEL.get(type_code, "mobile genetic element")
        notes.append("{label}: {name}".format(label=label, name=name))

    qualifiers = {"note": notes}
    if met is not None:
        qualifiers["mobile_element_type"] = ["{met}:{name}".format(met=met, name=name)]
    return qualifiers


def _location(start, end, strand, trunc_5p, trunc_3p):
    """1-based の start/end と strand, truncation から FeatureLocation を作る。

    biopython は 0-based half-open。truncation はストランドを考慮して
    5'/3' を genomic な左右端の partial(Before/After)に対応させる。

    MEF の trunc_5p は参照配列のアラインメント開始位置(1-based)で、5' 端まで
    完全なら 1。trunc_3p は 3' 側の未アラインメント塩基数で、完全なら 0。
    したがって 5' truncated は trunc_5p > 1、3' truncated は trunc_3p > 0。
    """
    left = int(start) - 1
    right = int(end)
    is_5p_trunc = int(trunc_5p) > 1
    is_3p_trunc = int(trunc_3p) > 0
    # strand=+1: 左端=5', 右端=3' / strand=-1: 左端=3', 右端=5'
    if strand == -1:
        left_trunc, right_trunc = is_3p_trunc, is_5p_trunc
    else:
        left_trunc, right_trunc = is_5p_trunc, is_3p_trunc
    left_pos = BeforePosition(left) if left_trunc else ExactPosition(left)
    right_pos = AfterPosition(right) if right_trunc else ExactPosition(right)
    return FeatureLocation(left_pos, right_pos, strand=strand)


def entry_to_feature(entry, index):
    """MEF JSON の result エントリ1件を (seq_id, ExtendedFeature) に変換する。"""
    seq_id = entry["contig"].split()[0]
    feature_key, _met, _putative = classify_mge(entry["type"], entry.get("evidence", 1))
    location = _location(entry["start"], entry["end"], entry["strand"],
                         entry.get("trunc_5p", 0), entry.get("trunc_3p", 0))
    feature = ExtendedFeature(location=location, type=feature_key,
                              id="MGE_{0}".format(index), seq_id=seq_id)
    feature.qualifiers = build_qualifiers(entry)
    return seq_id, feature


def parse_mge_results(data):
    """MEF JSON 全体を {seq_id: [ExtendedFeature, ...]} に変換する。"""
    result = {}
    for index, entry in enumerate(data.get("result", []), start=1):
        seq_id, feature = entry_to_feature(entry, index)
        result.setdefault(seq_id, []).append(feature)
    return result


class MobileElementFinder(ContigAnnotationTool):
    """MobileElementFinder

    Tool type: contig annotation (feature)
    URL: https://bitbucket.org/mhkj/mge_finder
    REF: Johansson et al., 2021
    """
    version = None
    TYPE = "feature"
    NAME = "MobileElementFinder"
    # mefinder は pkg_resources の DeprecationWarning を stderr に出す。
    # base_tools.setVersion() は stderr があると "not found" 扱いするため、
    # バージョンチェック時のみ stderr を捨てる（実行時は returncode で判定されるため影響なし）。
    VERSION_CHECK_CMD = ["mefinder", "--version", "2>/dev/null"]
    VERSION_PATTERN = r"(\d+\.\d+\.\d+)"

    def __init__(self, options=None, workDir="OUT"):
        if options is None:
            options = {}
        super(MobileElementFinder, self).__init__(options, workDir)
        self.cmd_options = options.get("cmd_options", "")
        self.db_path = options.get("db_path", "")
        self.output_directory = os.path.join(self.workDir, "contig_annotation", "mge_finder")
        self.result_file = os.path.join(self.output_directory, "mge.json")
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)

    def getCommand(self):
        prefix = os.path.join(self.output_directory, "mge")
        cmd = ["mefinder", "find",
               "-c", self.genomeFasta,
               "--json",
               "-t", str(self.options.get("cpu", 1)),
               "--temp-dir", self.output_directory]
        if self.db_path:
            cmd += ["--db-path", self.db_path]
        if self.cmd_options:
            cmd += self.cmd_options.split()
        cmd += [prefix]
        return cmd

    def getFeatures(self):
        """{seq_id: [ExtendedFeature, ...]} を返す。"""
        with open(self.result_file) as fh:
            data = json.load(fh)
        return parse_mge_results(data)

    def getResult(self):
        """ContigAnnotationTool 契約: (source_notes, report)。

        MEF は座標付き feature を getFeatures() で返すため source_notes は空。
        report には amr_summary.tsv の ## 行用の MGE 要約を入れる（PlasmidFinder 同形式）。
        type は読みやすい正式名（_TYPE_NAME）に変換して出力する。
        """
        source_notes = {}
        report = {}
        with open(self.result_file) as fh:
            data = json.load(fh)
        for entry in data.get("result", []):
            seq_id = entry["contig"].split()[0]
            type_name = _TYPE_NAME.get(entry["type"], entry["type"])
            summary = "MobileElementFinder: {name} ({type}), identity:{id:.1f}%, coverage:{cov:.1f}%, {s}..{e}".format(
                name=entry["name"], type=type_name,
                id=float(entry["identity"]) * 100, cov=float(entry["coverage"]) * 100,
                s=entry["start"], e=entry["end"])
            report.setdefault(seq_id, []).append(summary)
        return source_notes, report
