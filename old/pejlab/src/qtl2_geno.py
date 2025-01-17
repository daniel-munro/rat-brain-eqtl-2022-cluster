import sys
import os
import pysam
import pandas as pd
from bisect import bisect_left


def genotype_code(gt: tuple, founder: bool = False) -> str:
    if gt == (None, None):
        return "-"
    elif gt == (0, 0):
        return "A"
    elif (gt[0] == 0 and gt[1] > 0) or (gt[0] > 0 and gt[1] == 0):
        return "-" if founder else "H"
    elif gt[0] > 0 and gt[1] > 0:
        return "B"
    else:
        raise ValueError("GT not recognized: {}".format(gt))


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


SAMPLE_VCF = sys.argv[1]
FOUNDER_VCF = sys.argv[2]
SNPS = sys.argv[3]
GMAP_DIR = sys.argv[4]
GENO_OUT = sys.argv[5]
FOUNDER_OUT = sys.argv[6]
PMAP_OUT = sys.argv[7]
GMAP_OUT = sys.argv[8]

ID_list = open(SNPS, "r").read().splitlines()
IDs = set(ID_list)

maps = {}
for chrom in range(1, 21):
    filename = os.path.join(GMAP_DIR, "MAP4chr{}.txt".format(chrom))
    maps[chrom] = pd.read_table(filename, sep=" ", names=["pos", "ratio", "cm"])

vcf = pysam.VariantFile(SAMPLE_VCF)
samples = list(vcf.header.samples)
genos = {}
refs = {}
for rec in vcf.fetch():
    ID = rec.id if rec.id is not None else "{}:{}".format(rec.contig, rec.pos)
    if IDs is None or ID in IDs:
        gt = [s["GT"] for s in rec.samples.values()]
        labels = [genotype_code(g, founder=False) for g in gt]
        genos[ID] = labels
        refs[ID] = rec.ref

ID_list = [x for x in ID_list if x in genos.keys()]
IDs = set(ID_list)

vcf = pysam.VariantFile(FOUNDER_VCF)
strains = list(vcf.header.samples)
founder_genos = {}
ref_mismatch = 0
remove = set()
for rec in vcf.fetch():
    ID = rec.id if rec.id is not None else "{}:{}".format(rec.contig, rec.pos)
    if IDs is None or ID in IDs:
        gt = [s["GT"] for s in rec.samples.values()]
        labels = [genotype_code(g, founder=True) for g in gt]
        founder_genos[ID] = labels
        # assert rec.ref == refs[ID]
        if rec.ref != refs[ID]:
            print("Reference mismatch: {} {} vs. {}".format(ID, rec.ref, refs[ID]))
            ref_mismatch += 1
            remove.add(ID)
            del genos[ID]
            del founder_genos[ID]

if ref_mismatch > 0:
    print("{} SNPs removed due to reference mismatch.".format(ref_mismatch))
    ID_list = [ID for ID in ID_list if ID not in remove]

with open(GENO_OUT, "w") as out:
    out.write("id,{}\n".format(",".join(samples)))
    for ID in ID_list:
        out.write("{},{}\n".format(ID, ",".join(genos[ID])))

with open(FOUNDER_OUT, "w") as out:
    out.write("id,{}\n".format(",".join(strains)))
    for ID in ID_list:
        out.write("{},{}\n".format(ID, ",".join(founder_genos[ID])))

with open(PMAP_OUT, "w") as out:
    out.write("marker,chr,pos\n")
    for ID in ID_list:
        chrom, pos = tuple(ID.replace("chr", "").split(":"))
        pos = int(pos) / 1e6  # Units are Mbp.
        out.write("{},{},{}\n".format(ID, chrom, pos))

with open(GMAP_OUT, "w") as out:
    out.write("marker,chr,pos\n")
    for ID in ID_list:
        chrom, pos = tuple(ID.replace("chr", "").split(":"))
        gpos = genetic_pos(maps[int(chrom)], int(pos))
        out.write("{},{},{}\n".format(ID, chrom, round(gpos, 6)))
