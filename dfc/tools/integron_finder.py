#! /usr/bin/env python
# coding: UTF8

import csv
import glob
import os
from Bio.SeqFeature import FeatureLocation
from ..models.bio_feature import ExtendedFeature
from .base_tools import ContigAnnotationTool

# Integron types reported by IntegronFinder -> description for /note
_INTEGRON_TYPE_LABEL = {
    "complete": "complete integron",
    "In0": "In0 (integron integrase without attC sites)",
    "CALIN": "CALIN (cluster of attC sites lacking integron integrase)",
}


def read_integrons_file(file_name):
    """Read the .integrons table of IntegronFinder. Returns a list of dicts (one per element).

    The table starts with '#' comment lines; a file with no integron has only comments.
    """
    with open(file_name) as f:
        lines = [line for line in f if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def group_integrons(rows):
    """Group element rows by (replicon, integron id) and summarize each integron."""
    integrons = {}
    for row in rows:
        key = (row["ID_replicon"], row["ID_integron"])
        integron = integrons.setdefault(key, {
            "seq_id": row["ID_replicon"], "id": row["ID_integron"], "type": row["type"],
            "start": int(row["pos_beg"]), "end": int(row["pos_end"]), "attC": [], "intI": False})
        integron["start"] = min(integron["start"], int(row["pos_beg"]))
        integron["end"] = max(integron["end"], int(row["pos_end"]))
        if row["type_elt"] == "attC":
            integron["attC"].append((int(row["pos_beg"]), int(row["pos_end"]), int(row["strand"])))
        elif row["annotation"] == "intI":
            integron["intI"] = True
    return list(integrons.values())


def integron_summary(integron):
    label = _INTEGRON_TYPE_LABEL.get(integron["type"], integron["type"])
    return "{label}; attC sites: {n}; integrase: {intI}".format(
        label=label, n=len(integron["attC"]), intI="intI" if integron["intI"] else "none")


def integron_to_features(integron):
    """Return [mobile_element, misc_recomb (attC), ...] for one integron."""
    note = "IntegronFinder: {0}".format(integron_summary(integron))
    integron_feature = ExtendedFeature(
        location=FeatureLocation(integron["start"] - 1, integron["end"]), type="mobile_element",
        id="{0}_{1}".format(integron["seq_id"], integron["id"]), seq_id=integron["seq_id"],
        qualifiers={"mobile_element_type": ["integron"], "note": [note]})
    features = [integron_feature]
    for i, (start, end, strand) in enumerate(integron["attC"], 1):
        features.append(ExtendedFeature(
            location=FeatureLocation(start - 1, end, strand=strand), type="misc_recomb",
            id="{0}_{1}_attC{2}".format(integron["seq_id"], integron["id"], i), seq_id=integron["seq_id"],
            qualifiers={"note": ["attC site (IntegronFinder)"]}))
    return features


class IntegronFinder(ContigAnnotationTool):
    """IntegronFinder 2

    Tool type: contig annotation (feature)
    URL: https://github.com/gem-pasteur/Integron_Finder
    REF: Neron et al., 2022
    """
    version = None
    TYPE = "feature"
    NAME = "IntegronFinder"
    # integron_finder may write matplotlib cache messages to stderr, which setVersion() treats as "not found".
    VERSION_CHECK_CMD = ["integron_finder", "--version", "2>/dev/null"]
    VERSION_PATTERN = r"integron_finder version (\d+\.\d+\.\d+)"
    VERSION_ERROR_MSG = "IntegronFinder 2 is required for --modules integron. Install it with 'conda install -c bioconda integron_finder'."

    def __init__(self, options=None, workDir="OUT"):
        if options is None:
            options = {}
        super(IntegronFinder, self).__init__(options, workDir)
        self.cmd_options = options.get("cmd_options", "")
        self.output_directory = os.path.join(self.workDir, "contig_annotation", "integron_finder")
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)

    def getCommand(self):
        # ponytail: all replicons are searched as linear, so an integron spanning the origin of a circular
        # replicon is split. Pass a topology file (--topology-file) if this matters.
        cmd = ["integron_finder", "--linear", "--cpu", str(self.options.get("cpu", 1)),
               "--outdir", self.output_directory]
        if self.cmd_options:
            cmd += self.cmd_options.split()
        cmd += [self.genomeFasta]
        return cmd

    def getIntegrons(self):
        result_files = glob.glob(os.path.join(self.output_directory, "Results_Integron_Finder_*", "*.integrons"))
        rows = []
        for result_file in result_files:
            rows += read_integrons_file(result_file)
        return group_integrons(rows)

    def getFeatures(self):
        """Return {seq_id: [ExtendedFeature, ...]}."""
        result = {}
        for integron in self.getIntegrons():
            result.setdefault(integron["seq_id"], []).extend(integron_to_features(integron))
        return result

    def getResult(self):
        """ContigAnnotationTool contract: (source_notes, report). report holds the ## lines of amr_summary.tsv."""
        report = {}
        for integron in self.getIntegrons():
            report.setdefault(integron["seq_id"], []).append("IntegronFinder: {0}, {1}..{2}".format(
                integron_summary(integron), integron["start"], integron["end"]))
        return {}, report

    def getSummaryRows(self):
        return [{"module": "integron", "sequence": integron["seq_id"],
                 "location": "{0}..{1}".format(integron["start"], integron["end"]),
                 "gene/system": _INTEGRON_TYPE_LABEL.get(integron["type"], integron["type"]),
                 "note": "attC sites: {0}; integrase: {1}".format(len(integron["attC"]),
                                                                 "intI" if integron["intI"] else "none")}
                for integron in self.getIntegrons()]
