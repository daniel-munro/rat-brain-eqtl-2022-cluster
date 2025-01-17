import sys
import pysam
import pandas as pd


def genotype_code(gt: tuple) -> str:
    if gt == (None, None):
        return "NA"
    else:
        return str(gt[0] + gt[1])
        # raise ValueError("GT not recognized: {}".format(gt))


VCF = sys.argv[1]
BED = sys.argv[2]
GENOOUT = sys.argv[3]
PHENOOUT = sys.argv[4]
GENESOUT = sys.argv[5]
SNPSOUT = sys.argv[6]
SAMPLESOUT = sys.argv[7]

bed = pd.read_csv(BED, sep="\t", dtype=str)
# bed = bed.loc[bed["#chr"] == "12", :]
# samples = bed.columns[:]
# print(samples)
# quit()

vcf = pysam.VariantFile(VCF)
samples = [s for s in vcf.header.samples if s in bed.columns]

with open(GENOOUT, "w") as out:
    with open(SNPSOUT, "w") as snps:
        for rec in vcf.fetch():
            # if rec.contig != "12":
            #     continue
            ID = rec.id if rec.id is not None else "{}:{}".format(rec.contig, rec.pos)
            out.write(f"{ID},{rec.ref},{rec.alts[0]},")
            # geno = [genotype_code(s["GT"]) for s in rec.samples.values()]
            # geno = [g for i, g in enumerate(geno) if vcf.header.samples[i] in samples]
            geno = [rec.samples[sample] for sample in samples]
            geno = [genotype_code(s["GT"]) for s in geno]
            out.write(",".join(geno) + "\n")
            snps.write(f"{ID},{rec.pos},{rec.contig}\n")

bed[samples].T.to_csv(PHENOOUT, sep="\t", header=False, index=False)

# with open(GENESOUT, "w") as out:
#     for gene in bed["gene_id"]:
#         out.write(gene + "\n")

# Include gene locations to allow cis-window testing:
bed = bed[["gene_id", "#chr", "start"]]
bed = bed.rename(columns={"#chr": "chr", "start": "tss"})
bed.to_csv(GENESOUT, sep="\t", index=False)

with open(SAMPLESOUT, "w") as out:
    for sample in samples:
        out.write(sample + "\n")
