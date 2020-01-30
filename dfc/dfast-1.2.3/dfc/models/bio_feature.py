#! /usr/bin/env python
# coding: UTF8


from Bio.SeqFeature import SeqFeature


class ExtendedFeature(SeqFeature):
    def __init__(self, location=None, type='', location_operator='', strand=None, id='<unknown id>', qualifiers=None,
                 sub_features=None, ref=None, ref_db=None, primary_hit=None, seq_id=None, secondary_hits=None,
                 annotations=None):

        super(ExtendedFeature, self).__init__(location, type, location_operator, strand, id, qualifiers,
                                              sub_features, ref, ref_db)
        self.seq_id = seq_id
        self.primary_hit = primary_hit
        if not secondary_hits:
            secondary_hits = []
        self.secondary_hits = secondary_hits
        # self.partial_flag = flags  # value must be either of 00 01 10 11. 1 means partial, 0 means complete
        if annotations is None:
            annotations = {}
        self.annotations = annotations

    def assign_hit(self, verbosity=2):
        if self.type == "CRISPR":
            self.type = "repeat_region"
        if self.primary_hit:
            self.primary_hit.assign(self, verbosity)
        # if verbosity >= 2:
        for hit in self.secondary_hits:
            hit.assign_as_note(self, verbosity)
        if verbosity >= 3:
            self.qualifiers.setdefault("note", []).append(self.id)

if __name__ == '__main__':
    pass
