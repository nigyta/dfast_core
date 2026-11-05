import glob
import json
import os.path
import re
from Bio.SeqIO.InsdcIO import _insdc_location_string
from ..models.hit import defense_system_label

# One table of the accessory annotations: one row per gene, region or system.
HEADER = ["module", "sequence", "location", "locus_tag", "gene/system", "identity", "note"]
MODULE_ORDER = ["amr", "vf", "plasmid", "mge", "integron", "defense"]
NUCL_DB_MODULE = {"CARD": "amr", "VFDB": "vf"}


def nucl_hit_rows(genome, work_dir):
    """AMR/VF genes annotated by NuclSearch or NuclMapping (nucl_hits.json of each instance)"""
    rows = []
    for component in ("NuclSearch", "NuclMapping"):
        for json_file in glob.glob(os.path.join(work_dir, component + "*", "*.json")):
            for feature_id, hit in json.load(open(json_file)).items():
                feature = genome.features.get(feature_id)
                if feature is None:  # removed after the annotation
                    continue
                record = genome.seq_records[feature.seq_id]
                rows.append({"module": NUCL_DB_MODULE.get(hit["db_name"], hit["db_name"].lower()),
                             "sequence": feature.seq_id,
                             "location": _insdc_location_string(feature.location, len(record)),
                             "locus_tag": feature.qualifiers.get("locus_tag", [""])[0],
                             "gene/system": hit["model"]["gene"],
                             "identity": f"{hit['identity']:.1f}",
                             "note": f"{hit['model']['product']}; {hit['model']['note']}"})
    return rows


def defense_rows(genome, work_dir):
    """Systems found by DefenseFinder (systems.json): one row per system"""
    rows = []
    for json_file in glob.glob(os.path.join(work_dir, "DefenseFinder*", "systems.json")):
        for system in json.load(open(json_file)).values():
            genes = [(genome.features[fid], name) for fid, name in zip(system["features"], system["genes"])
                     if fid in genome.features]
            if not genes:
                continue
            features = [f for f, _ in genes]
            rows.append({"module": "defense", "sequence": features[0].seq_id,
                         "location": "{0}..{1}".format(min(int(f.location.start) for f in features) + 1,
                                                       max(int(f.location.end) for f in features)),
                         "locus_tag": ",".join(f.qualifiers.get("locus_tag", [f.id])[0] for f in features),
                         "gene/system": defense_system_label(system["type"], system["subtype"], system["activity"]),
                         "note": "genes: " + ", ".join(name for _, name in genes)})
    return rows


def summarize_accessory(genome, work_dir, contig_rows, output_file):
    """Write accessory_summary.tsv. contig_rows are the rows of the contig annotation tools (getSummaryRows)."""
    rows = nucl_hit_rows(genome, work_dir) + contig_rows + defense_rows(genome, work_dir)
    if not rows:
        return
    seq_order = {seq_id: i for i, seq_id in enumerate(genome.seq_records)}

    def sort_key(row):
        module = MODULE_ORDER.index(row["module"]) if row["module"] in MODULE_ORDER else len(MODULE_ORDER)
        return module, row["module"], seq_order.get(row["sequence"], len(seq_order)), int(re.search(r"\d+", row["location"]).group())

    with open(output_file, "w") as f:
        f.write("#" + "\t".join(HEADER) + "\n")
        for row in sorted(rows, key=sort_key):
            f.write("\t".join(str(row.get(key, "")) for key in HEADER) + "\n")
