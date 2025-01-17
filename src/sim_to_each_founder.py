import argparse
import sys
import os
import pysam
import pandas as pd
import numpy as np


def assemble_genotypes(gtypes):
    for gtrow in gtypes:
        for gt in gtrow:
            for i, g in enumerate(gt):
                if g is None:
                    gt[i] = -1
    gtypes = np.array(gtypes, dtype=int)
    return np.sort(gtypes, axis=2)


parser = argparse.ArgumentParser(description="Calculate similarity between individuals and founder strains")
parser.add_argument("pop_vcf", help="Population VCF file")
parser.add_argument("founder_vcf", help="Founder VCF file")
parser.add_argument("output", help="Output file (TSV)")
args = parser.parse_args()

IDs = set()

print("Reading population genotypes.", flush=True)
vcf = pysam.VariantFile(args.pop_vcf)
# vcf = pysam.VariantFile("data/genotype/phased/P50.rnaseq.88.unpruned.vcf.gz")
individuals = list(vcf.header.samples)
genos = {}
refs = {}
for rec in vcf.fetch():
    ID = rec.id if rec.id is not None else f"{rec.contig}:{rec.pos}"
    IDs.add(ID)
    # genos[ID] = [list(s["GT"]) for s in rec.samples.values()]
    genos[ID] = [list(rec.samples[ind]["GT"]) for ind in individuals]
    refs[ID] = rec.ref
    if len(IDs) % 1e5 == 0:
        print(ID, flush=True)

print("Reading founder genotypes.", flush=True)
vcf = pysam.VariantFile(args.founder_vcf)
# vcf = pysam.VariantFile("data/genotype/founders.vcf.gz")
strains = list(vcf.header.samples)
founder_genos = {}
ref_mismatch = 0
for rec in vcf.fetch():
    ID = rec.id if rec.id is not None else f"{rec.contig}:{rec.pos}"
    if ID in IDs:
        if rec.ref != refs[ID]:
            ref_mismatch += 1
        else:
            # founder_genos[ID] = [list(s["GT"]) for s in rec.samples.values()]
            founder_genos[ID] = [list(rec.samples[strain]["GT"]) for strain in strains]
            if len(founder_genos) % 1e5 == 0:
                print(ID, flush=True)

IDs.intersection_update(founder_genos.keys())

if ref_mismatch > 0:
    print(f"{ref_mismatch} SNPs removed due to reference mismatch.", flush=True)

print("Calculating similarity.", flush=True)
ID_list = list(IDs)
gt_pop = assemble_genotypes([genos[ID] for ID in ID_list])
gt_founder = assemble_genotypes([founder_genos[ID] for ID in ID_list])

simil = np.zeros((len(individuals), len(strains)))
for i, ind in enumerate(individuals):
    for j, strain in enumerate(strains):
        include = np.logical_and(
            np.sum(gt_pop[:, i, :] == -1, axis=-1) == 0,
            np.sum(gt_founder[:, j, :] == -1, axis=-1) == 0
        )
        same = np.sum(gt_pop[include, i, :] == gt_founder[include, j, :])
        simil[i, j] = same / (np.sum(include) * 2)

df = pd.DataFrame(simil.T, index=strains, columns=individuals)
df.to_csv(args.output, sep="\t", index_label="ID", float_format="%g")
# df.to_csv("data/tensorqtl/sim_to_founders.txt", sep="\t", index_label="ID", float_format="%g")
