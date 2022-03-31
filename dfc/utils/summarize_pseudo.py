from dfc.models.hit import ProteinHit, PseudoGene

header = ["# locus_id", "locus_tag", "location", "classification", "internal_stop", "indel",
          "ref_id", "description", "query_cov", "ref_cov", "identity"]

def summarize_pseudo(genome, output_file):
    output_buffer = [header]
    for key, feature in genome.features.items():
        # print(key, feature)
        hits = feature.secondary_hits + [feature.primary_hit]
        pseudogene_hits = [h for h in hits if isinstance(h, PseudoGene)]
        
        if pseudogene_hits:
            locus_tag = feature.qualifiers.get("locus_tag", ["-"])[0]
            assert len(pseudogene_hits) == 1
            pseudogene = pseudogene_hits[0]
            protein_hits = [h for h in hits if isinstance(h, ProteinHit) and h.id == pseudogene.ref_id ]
            # assert len(protein_hits) == 1
            p_hit = protein_hits[0]
            # print(locus_tag)
            # print(protein_hits)
            # for idx, hit in enumerate(protein_hits):
                # print(hit)
            # print(protein_hit)
            # print(p_hit.id, p_hit.description)
            # print(feature.location.start, feature.location.end, feature.location.strand)
            # print(round(p_hit.q_cov,1), p_hit.s_cov, p_hit.identity)
            strand = "+" if feature.location.strand == 1 else "-"
            internal_stop = [f"{x[0]+1}..{x[1]}({x[2]})".format(x[0] + 1, x[1], x[2]) for x in pseudogene.stop_codon]
            internal_stop = ",".join(internal_stop)
            # As of 1.2.16, insertion and deletion are integrated as indel
            indel = ",".join(map(str, pseudogene.indel))
            # insertion = ",".join(map(str, pseudogene.insertion))
            # deletion = ",".join(map(str, pseudogene.deletion))
            loc_string = f"{feature.seq_id}:{feature.location.start+1}..{feature.location.end}({strand})"
            classification = []
            if indel:
                classification.append("frameshift")
            if internal_stop:
                classification.append("internal_stop_codon")
            if (not classification) and "partial hit" in p_hit.flag:
                classification.append("partial")
            if classification:
                classification = ",".join(classification)
                ret = [feature.id, locus_tag, loc_string, classification, internal_stop, indel, p_hit.id, p_hit.description, round(p_hit.q_cov, 1), round(p_hit.s_cov, 1), round(p_hit.identity, 1)]
                ret = list(map(str, ret))
                # print("\t".join(map(str, ret)))
                output_buffer.append(ret)
        # print(feature.primary_hit)
        # print([str(hit) for hit in feature.secondary_hits])
        # for hit in feature.secondary_hits:
        #     print(type(hit))

    with open(output_file, "w") as f:
        for row in output_buffer:
            f.write("\t".join(row) + "\n")

if __name__ == '__main__':
    import pickle
    import dfc

    with open('OUT/genome.pickle', 'rb') as f:
          genome = pickle.load(f)
    summarize_pseudo(genome, "OUT/pseudogene_summary.tsv")