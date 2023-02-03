import os
import json
from datetime import datetime

def make_json_for_st_comments(dfast_version, config, tools, features):

    # initialization
    D = {}
    cnt_features = {}
    cnt_rrna = {"5S": 0, "16S": 0 , "23S": 0}
    gap_open = 0
    gap_length = 0

    # check annotated features
    target_features = ["CDS", "rRNA", "tRNA", "CRISPR"]
    annotation_method = []
    features_annotated = []
    for tool in tools:
        if tool.TYPE in target_features:
            annotation_method.append(f"{tool.__class__.NAME} v{tool.version}")
            features_annotated.append(tool.TYPE)
            cnt_features[tool.TYPE] = 0

    dt_now = datetime.now()
    D["Annotation Date"] = dt_now.strftime('%Y/%m/%d %H:%M:%S')
    if annotation_method:
        D["Annotation Method"] = "; ".join(annotation_method)
    D["Annotation Pipeline"] = "DDBJ Fast Annotation and Submission Tool (DFAST)"
    D["Annotation Provider"] = "DDBJ Center"
    D["Annotation Software revision"] = dfast_version
    if features_annotated:
        D["Features Annotated"] = "; ".join(features_annotated)

    # count annotated features
    for feature in features.values():
        if feature.type == "assembly_gap":
            gap_open += 1
            gap_length += len(feature)

        if feature.type in target_features:
            if feature.type == "rRNA":
                for rrna_type in cnt_rrna:
                    if rrna_type in feature.qualifiers.get("product", "")[0]:
                        cnt_rrna[rrna_type] += 1
            else:
                cnt_features[feature.type] += 1

    # Output feature counts
    for f_type in target_features:
        if f_type not in cnt_features:
            continue
        if f_type == "CRISPR":
            D["CRISPR Arrays"] = cnt_features[f_type]
        elif f_type == "rRNA":
            # i.e. 'rRNAs': '4, 4, 4 (5S, 16S, 23S)'
            D["rRNAs"] = ", ".join([str(cnt) for cnt in cnt_rrna.values()]) + " (" + ", ".join(cnt_rrna.keys()) + ")"
        else:
            D[f_type + "s"] = cnt_features[f_type]
    if gap_open > 0:
        D["Assembly Gaps"] = f"{gap_open} (Total: {gap_length} nt)"

    # output to json
    out_json_file = os.path.join(config.WORK_DIR, "st_comments_annotation.json")
    with open(out_json_file, "w") as f:
        json.dump(D, f, indent=4)

    # Debug
    # for key, value in D.items():
        # print(key, value, sep=":\t")
    # print("DEBUG")
    # for f in features.values():
    #     if f.type not in target_features:
    #         print(f)
