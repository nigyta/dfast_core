#!/usr/bin/env python3

import json
import re
from dataclasses import dataclass, field
from typing import List, Dict, Optional

# ========================================== 
# Experimental implementation. Will be removed in the future version.
# DFASTの結果をjson形式に変換するスクリプト。将来、モジュール化を行うい別のスクリプトに統合される予定
# ========================================== 


BOOL_QUALIFIERS = ["pseudo", "environmental_sample", "ribosomal_slippage", "circular_RNA", "proviral", "focus", "germline", "macronuclear", "circular"]  # circular is only for topology
# the list above may need to be updated. See https://www.ddbj.nig.ac.jp/ddbj/qualifiers.html


DATA_MODEL_VERSION = "0.1"


# ==========================================
# Data model (from part_ann.py in the prototype)
@dataclass
class Sequence:
    id: str
    seq: str

    def to_fasta(self, width=60, separator=False) -> str:
        fasta = []
        fasta.append(f">{self.id}")
        for i in range(0, len(self.seq), width):
            fasta.append(self.seq[i:i+width])
        fasta = "\n".join(fasta)
        if separator:
            fasta += "\n//"
        return fasta + "\n"

@dataclass
class Feature:
    """DDBJアノテーションのFeatureを表すクラス

    Attributes:
        type (str): Featureのタイプ（例：'source', 'CDS', 'gene'など）
        id (int | str): Featureの一意な識別子
        location (str): 配列上の位置情報（例：'1..300'）。COMMONタイプの場合は空文字列
        qualifiers (Dict[str, List[str | bool]]): Featureの修飾子（qualifier）を格納する辞書
            - キー：修飾子の名前（例：'product', 'note'など）
            - 値：修飾子の値のリスト。真偽値の修飾子の場合はbool型
    """
    type: str
    id: int # | str
    location: str = ""  # COMMONタイプの場合は空文字列
    # qualifiers: Dict[str, List[str | bool]] = field(default_factory=dict)
    qualifiers: Dict = field(default_factory=dict)

    def to_dict(self):
        """Convert the feature to a dictionary"""
        d = {
            "id": self.id,
            "type": self.type,
            "location": self.location,
            "qualifiers": self.qualifiers
        }
        locus_tag = self.qualifiers.get("locus_tag", [""])[0]
        if locus_tag:
            if "_" not in locus_tag:
                raise ValueError(f"Invalid locus_tag: {locus_tag}")
            del d["qualifiers"]["locus_tag"]
            d["locus_tag_id"] = locus_tag.split("_", 1)[1]
        return d

    def to_tsv(self) -> List[List[str]]:
        # ５列TSV形式に変換
        tsv = []
        for key, values in self.qualifiers.items():
            for value in values:
                if isinstance(value, bool):
                    tsv.append(["", "", "", key, ""])
                else:
                    tsv.append(["", "", "", key, value])
        # location がある場合はlocationを追加 (MSSファイルではlocationを持たないfeatureもある e.g. TOPOLOGYやDDBJ登録に関するもの)
        if self.location:
            tsv[0][2] = self.location
        tsv[0][1] = self.type
        return tsv

    def show(self):
        """Print the feature in a human-readable format"""
        indent = 2
        print(" " * indent + f"Feature {self.id}:")
        indent += 2
        print(" " * indent + f"Type: {self.type}")
        if self.location:
            print(" " * indent + f"Location: {self.location}")
        for key, values in self.qualifiers.items():
            for value in values:
                print(" " * indent + f"{key}: {value}")

@dataclass
class Entry:
    id: int # | str
    name: str = ""
    features: List[Feature] = field(default_factory=list)

    def to_tsv(self) -> List[List[str]]:
        # ５列TSV形式に変換
        tsv = []
        for feature in self.features:
            tsv.extend(feature.to_tsv())
        if tsv:
            tsv[0][0] = self.id
        return tsv

    def show(self):
        """Print the entry in a human-readable format"""
        print(f"Entry ID: {self.id} Name: {self.name}")
        print("Features:")
        for feature in self.features:
            feature.show()

@dataclass
class MSS:
    ann_file: str
    seq_file: str
    entries: List[Entry] = field(default_factory=list)
    sequences: List[Sequence] = field(default_factory=list)

    @staticmethod
    def parse(ann_file, seq_file) -> "MSS":
        """
        MSSのannファイルとseqファイルをパースしてMSSインスタンスを作成する       
        """
        current_entry = None
        current_feature: Optional[Feature] = None
        feature_counter = 1  # 連番カウンター
        mss = MSS(ann_file=ann_file, seq_file=seq_file)
        with open(mss.ann_file, "r") as f:
            for line in f:
                line = line.strip("\n")
                if not line:
                    continue
                # ５列に分割. entry, feature, location, qualifier_key, qualifier_value
                col_entry, col_feature, col_location, col_qkey, col_qvalue = line.split("\t", 4)
                
                if col_entry:
                    # 列1に値を持つ場合は新しいEntryを作成
                    current_entry = Entry(id=col_entry, name=col_entry)
                    mss.entries.append(current_entry)

                if col_feature:
                    # 列2に値を持つ場合は新しいFeatureを作成
                    current_feature = Feature(type=col_feature, id=f"feature_{feature_counter}")
                    feature_counter += 1
                    current_entry.features.append(current_feature)

                if col_location:
                    # 列3に値を持つ場合はcurrent_featureのlocationを設定
                    current_feature.location = col_location

                if col_qkey:
                    # current_featureにqualifierを追加
                    if col_qkey in BOOL_QUALIFIERS and col_qvalue == "":
                        # 値が空の場合はTrueを追加
                        current_feature.qualifiers.setdefault(col_qkey, []).append(True)
                    else:
                        current_feature.qualifiers.setdefault(col_qkey, []).append(col_qvalue)
        mss.parse_seq()
        return mss

    def parse_seq(self) -> None:
        """Parse the DDBJ sequence file"""
        current_seq = None
        
        with open(self.seq_file, "r") as f:
            entries = f.read()
            entries = entries.replace("/", "") # remove trailing slash (// in DDBJ format)
            entries = entries.split(">")
            entries = entries[1:]  # remove the first empty entry
            for entry in entries:
                lines = entry.split("\n")
                seq_id = lines[0]
                seq = "".join(lines[1:])
                self.sequences.append(Sequence(id=seq_id, seq=seq))



    def to_tsv(self) -> str:
        """Convert the parsed data to TSV format"""
        tsv = []
        for entry in self.entries:
            tsv.extend(entry.to_tsv())
        return "\n".join(["\t".join(row) for row in tsv])

    def to_fasta(self, separator=False) -> str:
        """Convert the parsed sequence data to FASTA format"""
        fasta = []
        for sequence in self.sequences:
            fasta.append(sequence.to_fasta(separator=separator))
        return "".join(fasta)

    def show(self):
        """Print the parsed data in a human-readable format"""
        print("\nEntries:")
        for entry in self.entries:
            entry.show()
        
        print("\nSequences:")
        for seq_id, sequence in self.sequences.items():
            print(f"\n{seq_id}:")
            print(f"  Length: {len(sequence)}")
            if len(sequence) > 10:
                print(f"  Sequence: {sequence[:10]}...")
            else:
                print(f"  Sequence: {sequence}")

    @staticmethod
    def from_json(json_file: str) -> "MSS":
        """
        JSONファイルからMSSインスタンスを作成する
        """

        with open(json_file) as f:
            data = json.load(f)

        mss = MSS(ann_file="", seq_file="")  # 空のMSSインスタンスを作成
        
        # COMMON エントリーの作成
        common_entry = Entry(id="COMMON", name="COMMON")
        feature_counter = 1
        
        # DBLINKの作成
        if "DBLINK" in data["COMMON"]:
            feature = Feature(
                id=f"feature_{feature_counter}",
                type="DBLINK"
            )
            feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                for k, v in data["COMMON"]["DBLINK"].items()}
            common_entry.features.append(feature)
            feature_counter += 1
        
        # SUBMITTERの作成
        if "SUBMITTER" in data["COMMON"]:
            feature = Feature(
                id=f"feature_{feature_counter}",
                type="SUBMITTER"
            )
            feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                for k, v in data["COMMON"]["SUBMITTER"].items()}
            common_entry.features.append(feature)
            feature_counter += 1
        
        # REFERENCEの作成
        if "REFERENCE" in data["COMMON"]:
            for ref in data["COMMON"]["REFERENCE"]:
                feature = Feature(
                    id=f"feature_{feature_counter}",
                    type="REFERENCE"
                )
                feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                    for k, v in ref.items()}
                common_entry.features.append(feature)
                feature_counter += 1
        
        # COMMENTの作成
        if "COMMENT" in data["COMMON"]:
            for comment in data["COMMON"]["COMMENT"]:
                feature = Feature(
                    id=f"feature_{feature_counter}",
                    type="COMMENT"
                )
                feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                    for k, v in comment.items()}
                common_entry.features.append(feature)
                feature_counter += 1
        
        # ST_COMMENTの作成
        if "ST_COMMENT" in data["COMMON"]:
            feature = Feature(
                id=f"feature_{feature_counter}",
                type="ST_COMMENT"
            )
            feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                for k, v in data["COMMON"]["ST_COMMENT"].items()}
            common_entry.features.append(feature)
            feature_counter += 1
        
        # DATEの作成
        if "DATE" in data["COMMON"]:
            feature = Feature(
                id=f"feature_{feature_counter}",
                type="DATE"
            )
            feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                for k, v in data["COMMON"]["DATE"].items()}
            common_entry.features.append(feature)
            feature_counter += 1
        
        # DATATYPEの作成
        if "DATATYPE" in data["COMMON"]:
            feature = Feature(
                id=f"feature_{feature_counter}",
                type="DATATYPE"
            )
            feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                for k, v in data["COMMON"]["DATATYPE"].items()}
            common_entry.features.append(feature)
            feature_counter += 1
        
        # KEYWORDの作成
        if "KEYWORD" in data["COMMON"]:
            feature = Feature(
                id=f"feature_{feature_counter}",
                type="KEYWORD"
            )
            feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                for k, v in data["COMMON"]["KEYWORD"].items()}
            common_entry.features.append(feature)
        
        # common_entryを先頭に追加
        mss.entries.append(common_entry)
        
        # COMMON_SOURCEの内容を取得
        common_source = data.get("COMMON_SOURCE", {})
        
        # エントリーの作成
        for entry_data in data["ENTRIES"]:
            entry_id = entry_data["id"]
            entry = Entry(id=entry_id, name=entry_data["name"])
            features = []
            
            # source featureを作成し、COMMON_SOURCEの内容を追加
            for feature_data in entry_data["features"]:
                if feature_data["type"] == "source":
                    feature_id = feature_data["id"]
                    feature = Feature(
                        id=feature_id,
                        type="source",
                        location=feature_data["location"]
                    )
                    # COMMON_SOURCEの内容を追加
                    qualifiers = dict(common_source)
                    # エントリー固有のqualifierを追加
                    qualifiers.update({k: v for k, v in feature_data["qualifiers"].items()})
                    feature.qualifiers = {k: [v] if not isinstance(v, list) else v 
                                       for k, v in qualifiers.items()}
                    features.append(feature)
                    break
            
            # topologyがcircularの場合、TOPOLOGY featureを追加
            if entry_data.get("_topology") == "circular":
                topology_feature = Feature(
                    id=f"topology_{entry_id}",
                    type="TOPOLOGY",
                    location=""
                )
                topology_feature.qualifiers = {"circular": [True]}
                features.append(topology_feature)
            
            # その他のfeatureを追加
            for feature_data in entry_data["features"]:
                feature_id = feature_data["id"]
                if feature_data["type"] != "source":  # sourceは既に追加済み
                    feature = Feature(
                        id=feature_id,
                        type=feature_data["type"],
                        location=feature_data["location"]
                    )
                    feature.qualifiers = feature_data["qualifiers"]
                    features.append(feature)
            
            # 作成したfeaturesをentryに設定
            entry.features = features
            
            # シーケンスの作成
            sequence = Sequence(id=entry_id, seq=entry_data["sequence"])
            mss.sequences.append(sequence)
            mss.entries.append(entry)
        
        return mss
# ==========================================


def common_to_dict(mss):
    # COMMON entryの情報を取得して辞書で返す
    qualifiers_with_multi_values = ["line", "ab_name", "keyword", "sequence read archive"]
    features_with_multi_values = ["COMMENT", "REFERENCE"]
    common_entry = [entry for entry in mss.entries if entry.id == "COMMON"]
    len_common_entry = len(common_entry)
    if len_common_entry == 0:
        raise ValueError("COMMON entry not found.")
    elif len_common_entry > 1:
        raise ValueError("Multiple COMMON entries found.")
    common_entry = common_entry[0]
    D = {}
    for feature in common_entry.features:
        inner_dict = {}
        for key, value in feature.qualifiers.items():
            if key in qualifiers_with_multi_values:
                inner_dict[key] = value
            else:
                assert len(value) == 1
                inner_dict[key] = value[0]
        if feature.type in features_with_multi_values:
            D.setdefault(feature.type, []).append(inner_dict)
        else:
            D[feature.type] = inner_dict
    return {"COMMON": D}

def source_to_dict(mss):
    """
    source featureの情報を取得して、共通する項目をCOMMON_SOURCEとして辞書で返す
    COMMON_SOURCEに含まれるqualifierは各entryのsource featureからは削除する
    """
    qualifiers_with_multi_values = ["note"]
    ignore_qualifiers = ["submitter_seqid", "ff_definition"]
    source_features = [feature for entry in mss.entries for feature in entry.features if feature.type == "source"]

    # すべてのenrtyについて等しい値を持つqualifierを抽出する
    common_qualifiers = {}
    for feature in source_features:
        if feature.type != "source":
            continue
        for key, value in feature.qualifiers.items():
            if key not in ignore_qualifiers:
                value = "//".join(value)  # 一時的にリストを文字列に変換
                common_qualifiers.setdefault(key, []).append(value)
    common_qualifiers = {key: value for key, value in common_qualifiers.items() if len(value) == len(source_features) and len(set(value))==1}

    # 各entryのsource featureからCOMMON_SOURCEに含まれるqualifierを削除
    common_qualifiers_keys = set(common_qualifiers.keys())
    for source_feature in source_features:
        for key in common_qualifiers_keys:
            if key in source_feature.qualifiers:
                del source_feature.qualifiers[key]

    # COMMON_SOURCEとして返す辞書を作成
    D = {}
    for key, value in common_qualifiers.items():
        value = value[0]
        if key in qualifiers_with_multi_values:
            D[key] = value.split("//")
        else:
            D[key] = value
    return {"COMMON_SOURCE": D}

def infer_meta_from_ann(mss):
    """
    annファイルに記載された情報から、
    trad_submission_category: GNM; WGS
    seq_type: chromosome; plasmid; unplaced; other
    seq_topology: linear; circular
    を決定する
    common_meta と seq_info の2つの辞書を返す
    注意: DFASTが出力するannファイルにのみ対応
    """

    def _get_locus_tag_prefix(entries): # -> str|None:
        """
        locus_tag_prefix の決定
        entry, featureのqualifierにlocus_tagが含まれる場合、その値の _ より前の部分をlocus_tag_prefixとする
        """
        for entry in entries:
            for feature in entry.features:
                if "locus_tag" in feature.qualifiers:
                    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                    if locus_tag:
                        return locus_tag.split("_")[0]
        return None

    def _get_dfast_version(entries): # -> str|None:
        """
        DFASTのバージョンを取得する
        COMMENT		line	Annotated by DFAST v.1.3.4 https://dfast.ddbj.nig.ac.jp/
        の1.3.4の部分を取得する
        """
        dfast_version_pat = re.compile(r"DFAST (ver\.|ver|version|v|v\.)?\s?(\d+\.\d+\.\d+)")
        for entry in entries:
            if entry.id == "COMMON":
                for feature in entry.features:
                    if feature.type == "COMMENT":
                        for qualifier_value in feature.qualifiers.get("line", []):
                            m = dfast_version_pat.search(qualifier_value)
                            if m:
                                dfast_version = m.group(2)
                                # print(f"dfast_version: {dfast_version}, match1 '{m.group(1)}', match2 '{m.group(2)}'")  # debug
                                return dfast_version
        return None

    dfast_version = _get_dfast_version(mss.entries)
    if dfast_version:
        common_meta = {"dfast_version": dfast_version, "division": "BCT"}
    else:
        common_meta = {"division": "BCT"}  # For DFAST, the division is always BCT

    # trad_submission_category (WGS, GNM) かの判定
    # COMMON entryのDATATYPE featureのtypeがWGSであればWGSとする
    # それ以外はGNMとする 
    trad_submission_category = "GNM"
    for entry in mss.entries:
        if entry.id == "COMMON":
            for feature in entry.features:
                if feature.type == "DATATYPE" and "WGS" in feature.qualifiers.get("type", []):
                    trad_submission_category = "WGS"
                    break
            break
    common_meta["trad_submission_category"] = trad_submission_category
    # print(f"_trad_submission_category: {trad_submission_category}")  # debug

    # seq_prefixの決定 (wgsの場合のみ)
    # sequence01 contig01 等から数字の部分を取り除いたもの。
    seq_prefix = ""
    if trad_submission_category == "WGS":
        for entry in mss.entries:
            if entry.id == "COMMON":
                continue
            seq_name = entry.name
            seq_prefix = seq_name.rstrip("0123456789")
            break
    if seq_prefix:
        common_meta["seq_prefix"] = seq_prefix
    # print(f"seq_prefix: {seq_prefix}")  # debug

    # locus_tag_prefix の決定
    locus_tag_prefix = _get_locus_tag_prefix(mss.entries)
    if locus_tag_prefix:
        common_meta["locus_tag_prefix"] = locus_tag_prefix

    # seq_infoの決定
    # gnm (complete) の場合、 seq_type, seq_topology
    # wgs (draft) の場合、不要
    seq_info = {}
    if trad_submission_category == "GNM":
        for entry in mss.entries:
            source_or_tpoology = [feature for feature in entry.features if feature.type in ["source", "TOPOLOGY"]]
            seq_topology, seq_type = "linear", "other"
            for feature in source_or_tpoology:
                if feature.type == "source":
                    ff_definition = feature.qualifiers.get("ff_definition", [""])[0]
                    if "plasmid" in ff_definition:
                        seq_type = "plasmid"
                    elif "unplaced sequence" in ff_definition:
                        seq_type = "unplaced"
                    elif "complete genome" in ff_definition:
                        seq_type = "chromosome"

                # seq_topologyの決定 TOPOLOGY featureが記載されていればcircularとする
                if feature.type == "TOPOLOGY" and "circular" in feature.qualifiers:
                    seq_topology = "circular"
            seq_info[entry.id] = {"seq_type": seq_type, "seq_topology": seq_topology}


    return {"COMMON_META": common_meta}, seq_info

def make_entry_dict(mss, seq_info):
    """
    mssオブジェクトを受け取り、各entryの情報を辞書で返す
    """
    # D = {}
    L = []
    sequence_dict = {sequence.id: sequence.seq for sequence in mss.sequences}
    for entry in mss.entries:
        if entry.id == "COMMON":
            continue
        sequence = sequence_dict[entry.id]  # note: sequence.id and entry.id is always identical, but entry.name may be different.
        # sequence = sequence[:10] + "..."  # for debug
        seq_info_dict = seq_info.get(entry.id, {})
        seq_type, seq_topology = seq_info_dict.get("seq_type", "other"), seq_info_dict.get("seq_topology", "linear")
        features = [feature.to_dict() for feature in entry.features if feature.type != "TOPOLOGY"]
        # features = {feature.id: feature.to_dict() for feature in entry.features if feature.type != "TOPOLOGY"}
        # D[entry.id] = {"name": entry.name, "_type": seq_type, "_topology": seq_topology, "sequence": sequence, "features": features}
        L.append({"id": entry.id, "name": entry.name, "type": seq_type, "topology": seq_topology, "sequence": sequence, "features": features})
    return {"ENTRIES": L}


def ann2json_for_dfast(ann_file, seq_file, out_json_file=None):
    """
    DFAST が生成する ann ファイルと seq ファイルを読み込み、json形式のファイルに変換する

    todo: CDS featureのtranslationが含まれていないので, jsonにtranslationを含めるか、jsonの情報を元にtranslationを生成する機能が必要
    """
    dat = {"schema_version": DATA_MODEL_VERSION}

    # MSSインスタンスを作成してパース
    mss = MSS.parse(ann_file, seq_file)
    

    common_dict = (common_to_dict(mss))

    common_source = source_to_dict(mss)

    common_meta, seq_info = infer_meta_from_ann(mss)
    trad_submission_category = common_meta["COMMON_META"].pop("trad_submission_category")
    common_dict["COMMON"]["trad_submission_category"] = trad_submission_category

    dict_entry = make_entry_dict(mss, seq_info)

    dat.update(common_dict)
    dat.update(common_source)
    dat.update(common_meta)
    dat.update(dict_entry)

    if out_json_file is None:
        print(json.dumps(dat, indent=4))
    else:
        with open(out_json_file, "w") as f:
            json.dump(dat, f, indent=4)



if __name__ == "__main__":
    import sys

    ann_file = sys.argv[1]
    seq_file = sys.argv[2]

    if len(sys.argv) > 3:
        out_json_file = sys.argv[3]
    else:
        out_json_file = None
    ann2json_for_dfast(ann_file, seq_file, out_json_file)