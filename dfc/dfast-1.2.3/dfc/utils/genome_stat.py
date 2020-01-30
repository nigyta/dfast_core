#! /usr/bin/env python
# coding: UTF8

import os
from Bio import SeqIO
from collections import OrderedDict
import json
from ..genome import Genome

class GenomeStat(object):

    def __init__(self):
        self.C = 0
        self.G = 0
        self.A = 0
        self.T = 0
        self.N = 0
        self.lengthList = []
        self.cds = 0
        self.cdsLength = 0
        self.rRNA = 0
        self.tRNA = 0
        self.CRISPR = 0

    def update(self, record):
        seq = record.seq.upper()
        self.C += seq.count("C")
        self.G += seq.count("G")
        self.A += seq.count("A")
        self.T += seq.count("T")
        self.N += seq.count("N")
        self.lengthList.append(len(seq))

        cdsList = [len(f) for f in record.features if f.type == "CDS"]
        self.cds += len(cdsList)
        self.cdsLength += sum(cdsList)

        self.rRNA += len([f for f in record.features if f.type == "rRNA"])
        self.tRNA += len([f for f in record.features if f.type == "tRNA"])
        self.CRISPR += len([f for f in record.features if (f.type == "repeat_region" or f.type == "CRISPR")])

    def finalize(self):
        self.totalLength = sum(self.lengthList)
        self.seqNumber = len(self.lengthList)
        self.lengthList.sort(reverse=True)
        self.longest = self.lengthList[0]
        tmpLength = 0
        for x in self.lengthList:
            tmpLength += x
            if tmpLength >= int(self.totalLength / 2):
                self.N50 = x
                break
        if self.totalLength > 0:
            self.Npercent = 100.0 * self.N / self.totalLength
            self.codingRatio = 100.0 * self.cdsLength / self.totalLength
        else:
            self.Npercent = 0
            self.codingRatio = 0
        if self.totalLength - self.N > 0:
            self.GCcontent = 100.0 * (self.G + self.C) / (self.totalLength - self.N)
        else:
            self.GCcontent  = 0
        if self.cds > 0:
            self.aveProtLength = 1.0 * self.cdsLength / self.cds / 3
        else:
            self.aveProtLength = 0

    def output(self, outputFileName):
        with open(outputFileName, "w") as f:
            f.write("%s\t%d\n" % ("Total Sequence Length (bp)", self.totalLength))
            f.write("%s\t%d\n" % ("Number of Sequences", self.seqNumber))
            f.write("%s\t%d\n" % ("Longest Sequences (bp)", self.longest))
            f.write("%s\t%d\n" % ("N50 (bp)", self.N50))
            f.write("%s\t%f\n" % ("Gap Ratio (%)", self.Npercent))
            f.write("%s\t%.1f\n" % ("GCcontent (%)", self.GCcontent))
            f.write("%s\t%d\n" % ("Number of CDSs", self.cds))
            f.write("%s\t%.1f\n" % ("Average Protein Length", self.aveProtLength))
            f.write("%s\t%.1f\n" % ("Coding Ratio (%)", self.codingRatio))
            f.write("%s\t%d\n" % ("Number of rRNAs", self.rRNA))
            f.write("%s\t%d\n" % ("Number of tRNAs", self.tRNA))
            f.write("%s\t%d\n" % ("Number of CRISPRs", self.CRISPR))

    def toJson(self, outputFileName):
        D = OrderedDict()
        D["total_sequence_length"] = self.totalLength
        D["total_sequence_number"] = self.seqNumber
        D["longest_sequence_length"] = self.longest
        D["n50"] = self.N50
        D["gc_content"] = self.GCcontent
        D["gap_percent"] = self.Npercent
        D["number_of_proteins"] = self.cds
        D["average_protein_length"] = self.aveProtLength
        D["coding_ratio"] = self.codingRatio
        D["number_of_rrna"] = self.rRNA
        D["number_of_trna"] = self.tRNA
        D["number_of_crisprs"] = self.CRISPR
        with open(outputFileName, "w") as f:
            json.dump(D, f, indent=4, sort_keys=False, separators=(',', ': '))

    def print2screen(self):
        print("%s\t%d" % ("Assembled_Total_Length", self.totalLength))
        print("%s\t%d" % ("Assembled_Sequence_Number", self.seqNumber))
        print("%s\t%d" % ("Assembled_LongestSequence", self.longest))
        print("%s\t%d" % ("N50", self.N50))
        print("%s\t%f" % ("N_percent", self.Npercent))
        print("%s\t%.1f" % ("GCcontent", self.GCcontent))
        print("%s\t%d" % ("Number_of_Proteins", self.cds))
        print("%s\t%.1f" % ("Average_Protein_Length", self.aveProtLength))
        print("%s\t%.1f" % ("CodingRatio", self.codingRatio))
        print("%s\t%d" % ("Number_of_rRNA", self.rRNA))
        print("%s\t%d" % ("Number_of_tRNA", self.tRNA))
        print("%s\t%d" % ("Number_of_CRISPRs", self.CRISPR))

    @staticmethod
    def from_file(genbankFileName, output_file):
        gs = GenomeStat()
        for record in SeqIO.parse(genbankFileName, "genbank"):
            gs.update(record)
        gs.finalize()
        gs.output(output_file)

    @staticmethod
    def execute(genome):
        output_file = os.path.join(genome.workDir, "statistics.txt")
        gs = GenomeStat()
        for record in genome.seq_records.values():
            gs.update(record)
        gs.finalize()
        gs.output(output_file)


if __name__ == '__main__':
    pass
