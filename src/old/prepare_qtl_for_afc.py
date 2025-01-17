# Extract from gzipped FastQTL output file the first (gene ID) and
# seventh (variant ID) column. Then extract chromosome and location
# from variant IDs and include those too. Include column indices to
# handle other eQTL files.

import sys
import gzip
import re

infile = sys.argv[1]
# outfile = sys.argv[2]
if len(sys.argv) > 2:
    gene_col = int(sys.argv[2]) - 1
    snp_col = int(sys.argv[3]) - 1
else:
    gene_col = 0
    snp_col = 6

with gzip.open(infile, "rt") as f:
    # with open(outfile, "w") as out:
    lines = f.read().splitlines()
    # out.write("pid\tsid\tsid_chr\tsid_pos\n")
    print("pid\tsid\tsid_chr\tsid_pos")
    for line in lines[1:]:
        items = line.split("\t")
        m = re.match("^chr(\\d+):(\\d+)$", items[snp_col])
        chrom, pos = m.group(1), m.group(2)
        # out.write("{}\t{}\t{}\t{}\n".format(items[0], items[6], chrom, pos))
        print("{}\t{}\t{}\t{}".format(items[gene_col], items[snp_col], chrom, pos))
