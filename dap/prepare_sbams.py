import argparse
import os
import pandas as pd
import pysam


def genotype_code(gt: tuple) -> str:
    if gt == (None, None):
        return "NA"
    else:
        return str(gt[0] + gt[1])


p = argparse.ArgumentParser(description="Read VCF and expression BED and output sbams files for DAP-G.")
p.add_argument("vcf", help="VCF file")
p.add_argument("bed", help="Expression BED file")
# p.add_argument("gene", help="Name of gene in BED file for which to prepare sbams file")
p.add_argument("outdir", help="Name of directory to write output files ({gene}.sbams.dat)")
p.add_argument("--covar", help="Covariates file (TSV, one row per covariate)")
args = p.parse_args()

vcf = pysam.VariantFile(args.vcf)
bed = pd.read_csv(args.bed, sep="\t", dtype=str)
if args.covar is not None:
    covar = pd.read_csv(args.covar, sep="\t", dtype=str)
samples = bed.columns[4:]
for i in range(bed.shape[0]):
    gene = bed["gene_id"].iloc[i]
    out = open(os.path.join(args.outdir, f"{gene}.sbams.dat"), "w")
    # bed = bed.loc[bed["gene_id"] == args.gene, :]
    # assert len(bed["gene_id"]) == 1
    tss = (bed["#chr"].iloc[i], int(bed["start"].iloc[i]))
    exprs = "\t".join(bed.iloc[i, 4:])
    out.write(f"pheno\t{gene}\tgroup\t{exprs}\n")

    for rec in vcf.fetch(tss[0], max(0, tss[1] - 1e6), tss[1] + 1e6):
        geno = [rec.samples[sample] for sample in samples]
        geno = [genotype_code(s["GT"]) for s in geno]
        genos = "\t".join(geno)
        out.write(f"geno\t{rec.id}\tgroup\t{genos}\n")

    if args.covar is not None:
        for i in range(covar.shape[0]):
            values = [covar[sample].iloc[i] for sample in samples]
            values = "\t".join(values)
            out.write(f"controlled\t{covar.ID[i]}\tgroup\t{values}\n")
