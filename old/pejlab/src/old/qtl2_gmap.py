import sys
import os
import pandas as pd
# import pysam
from bisect import bisect_left

# VCF = sys.argv[1]
SNPS = sys.argv[1]
OUTFILE = sys.argv[2]
MAP_DIR = sys.argv[3]


def genetic_pos(chrmap: pd.DataFrame, pos: int) -> float:
    r = bisect_left(chrmap["pos"], pos)
    if r == len(chrmap["pos"]):
        return chrmap["cm"][r - 1]
    elif chrmap["pos"][r] == pos or r == 0:
        return chrmap["cm"][r]
    else:
        # Interpolate the genetic position.
        p_lo = chrmap["pos"][r - 1]
        p_hi = chrmap["pos"][r]
        g_lo = chrmap["cm"][r - 1]
        g_hi = chrmap["cm"][r]
        rel = (pos - p_lo) / (p_hi - p_lo)
        return g_lo + rel * (g_hi - g_lo)


maps = {}
for chrom in range(1, 21):
    filename = os.path.join(MAP_DIR, "MAP4chr{}.txt".format(chrom))
    maps[chrom] = pd.read_table(filename, sep=" ", names=["pos", "ratio", "cm"])

# vcf = pysam.VariantFile(VCF)
IDs = open(SNPS, "r").read().splitlines()

with open(OUTFILE, "w") as out:
    out.write("marker,chr,pos\n")
    # for rec in vcf.fetch():
    for ID in IDs:
        # ID = rec.id if rec.id is not None else "{}:{}".format(rec.contig, rec.pos)
        chrom, pos = tuple(ID.replace("chr", "").split(":"))
        # pos = genetic_pos(maps[int(rec.contig)], rec.pos)
        gpos = genetic_pos(maps[int(chrom)], int(pos))
        # out.write("{},{},{}\n".format(ID, rec.contig, round(pos, 6)))
        out.write("{},{},{}\n".format(ID, chrom, round(gpos, 6)))
