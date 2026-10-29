import os.path
import glob
import json
from Bio.SeqIO.InsdcIO import _insdc_location_string
from ..models.hit import defense_system_label

def summarize_amr(genome, work_dir, contig_annotation_report, output_file):

    nucl_hits_json_files = []
    for component in ("NuclSearch", "NuclMapping"):
        nucl_hits_json_files += glob.glob(os.path.join(work_dir, component + "*", "*.json"))

    nucl_hits = {}
    features_with_hits = {}
    for json_file in nucl_hits_json_files:
        dat = json.load(open(json_file))
        for feature_id, hit in dat.items():
            feature = genome.features.get(feature_id)
            if feature is None:
                # feature is deleted (partial, pseudo, or something else)
                # todo: output log
                continue

            seq_id = feature.seq_id
            if feature_id not in features_with_hits.setdefault(seq_id, []):
                features_with_hits[seq_id].append(feature_id)
            nucl_hits.setdefault(feature_id, []).append(hit)

    # Defense systems found by DefenseFinder: one ## line per system on the sequence of its genes
    for json_file in glob.glob(os.path.join(work_dir, "DefenseFinder*", "systems.json")):
        for system in json.load(open(json_file)).values():
            features = [genome.features[fid] for fid in system["features"] if fid in genome.features]
            if not features:
                continue
            label = defense_system_label(system["type"], system["subtype"], system["activity"])
            genes = ", ".join(f.qualifiers.get("locus_tag", [f.id])[0] for f in features)
            contig_annotation_report.setdefault(features[0].seq_id, []).append(f"DefenseFinder: {label}, genes: {genes}")

    header = ["locus", "location", "hit_accession", "gene", "product", "identity", "q_cov", "s_cov", "e_value", "note"]
    output_buffer = ""

    for rec in genome.seq_records.values():
        for contig_annotation in contig_annotation_report.get(rec.id, []):
            output_buffer += f"## {rec.id} {contig_annotation}\n"
        for feature_id in features_with_hits.get(rec.id, []):
            for nucl_hit in nucl_hits[feature_id]:
                # print(feature_id)
                feature = genome.features.get(feature_id)
                if feature is None:
                    # feature is deleted (partial, pseudo, or something else)
                # todo: output log
                    continue
                locus_tag = feature.qualifiers.get("locus_tag", [""])[0] or feature_id
                location = rec.id + ":" + _insdc_location_string(feature.location, len(rec))
                accession = nucl_hit["model"]["accession"]
                gene = nucl_hit["model"]["gene"]
                product = nucl_hit["model"]["product"]
                identity, q_cov, s_cov = f"{nucl_hit['identity']:.1f}", f"{nucl_hit['q_cov']:.1f}", f"{nucl_hit['s_cov']:.1f}"
                e_value = f"{nucl_hit['e_value']}"
                note = nucl_hit["model"]["note"]
                out_list = [locus_tag, location, accession, gene, product, identity, q_cov, s_cov, e_value, note]
                output_buffer += "\t".join(out_list) + "\n"
            # print(feature_id, nucl_hits[feature_id])
            # amr_hit = nucl_hits[]
            # out_list = [rec.id, ]
            # output_buffer += f"{amr_hit}\n"

    if output_buffer:
        output_buffer = "#" + "\t".join(header) + "\n" + output_buffer

        with open(output_file, "w") as f:
            f.write(output_buffer)

if __name__ == '__main__':
    import pickle
    import dfc

    # with open('OUT/genome.pickle', 'rb') as f:
    #       genome = pickle.load(f)
    # summarize_pseudo(genome, "OUT/pseudogene_summary.tsv")