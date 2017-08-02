#! /usr/bin/env python
# coding: UTF8

class CDDmodel(object):

    def __init__(self, id_, accession, short_name, description, length):
        self.id, self.accession, self.short_name, self.description, self.length = (
            id_, accession, short_name, description, length)

    @staticmethod
    def read_cdd_definition(file_name):
        D = {}
        with open(file_name) as f:
            data = f.readlines()
            for row in data:
                id_, accession, short_name, description, length = row.strip("\n").split("\t")
                cdd_model = CDDmodel(id_, accession, short_name, description, length)
                D[cdd_model.id] = cdd_model
        return D

if __name__ == '__main__':
    pass
