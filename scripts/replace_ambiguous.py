# replace_ambiguous
# Replace ambiguous bases with Ns. [i.e. R (=A or G), Y (=T or C), and so on.]

from __future__ import print_function
import sys
import os

if not len(sys.argv) == 2:
    sys.stderr.write("\n\tUsage: python replace_ambiguous.py file_name.fasta > output_file_name.fasta\n\n")
    exit()

# ambiguous bases
# http://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
targets = list("RYKMSWBDHVX-")

file_name = sys.argv[1]
if not os.path.exists(file_name):
    sys.stderr.write("File not found: " + file_name + "\n")
    exit()

for line in open(file_name):
    if line.startswith(">"):
        print(line, end="")
    else:
        line = line.upper()
        for base in targets:
            line = line.replace(base, "N")
        print(line, end="")

with open(file_name) as f:
    dat = f.readlines()
    dat = "".join([line for line in dat if not line.startswith(">")]).upper()
    for base in targets:
        cnt = dat.count(base)
        if cnt > 0:
            sys.stderr.write("Replaced {} {}s with N.\n".format(cnt, base))
