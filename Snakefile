import pandas as pd

tissues = ["IL", "LHb", "NAcc", "OFC", "PL"]
fastqs = pd.read_csv("data/fastq/fastq_files_listing.fixed.txt", sep="\t")

# sample_ids = {}
# for tissue in tissues:
#     fastqs_tissue = fastqs.loc[fastqs["brain_region"] == tissue]
#     sam_tissue = fastqs_tissue["library"].unique()
#     sample_ids[tissue] = [x for x in sam_tissue if x not in exclude[tissue]]
#     sample_ids[tissue] = sam_tissue
samples = pd.read_csv("data/samples.txt", sep="\t")
sample_ids = {}
for tissue in tissues:
    sam_tissue = samples.loc[samples["brain_region"] == tissue, :]
    sample_ids[tissue] = sam_tissue["library"].to_numpy()


localrules:
    individual_vcf,
    observed_snp_list,
    vcf_chr_list,
    index_vcf,
    index_bam,

wildcard_constraints:
    rat_id = "[A-Z0-9]+",
    tissue = "[A-Za-z]+",
    method = "[a-z0-9]+",
    chr = "chr[0-9]+",
    code = "[A-Za-z]+",

include: "align.smk"
include: "ASE.smk"
include: "tensorqtl.smk"
# include: "qtl2.smk"
include: "afc.smk"
include: "splice.smk"


rule all:
    input:
        # expand("data/markdup_out/IL/{sample_id}.bam.bai", sample_id=sample_ids["IL"]),
        # expand("data/markdup_out/LHb/{sample_id}.bam.bai", sample_id=sample_ids["LHb"]),
        # expand("data/markdup_out/NAcc/{sample_id}.bam.bai", sample_id=sample_ids["NAcc"]),
        # expand("data/phaser_out/{sample_id}.gene_ae.txt", sample_id=sample_ids),
        # expand("data/phaser_pop_out/{tissue}.expr_matrix.gw_phased.bed.gz", tissue=tissues),
        # expand("data/tensorqtl/{code}.cis_qtl.txt.gz", code=["IQCT", "LQCT", "NQCT", "OQCT", "PQCT"]),
        # expand("data/tensorqtl/{code}.cis_independent_qtl.txt.gz", code=["IQCT", "LQCT", "NQCT", "OQCT", "PQCT"]),
        # expand("data/tensorqtl/{code}.cis_qtl_signif.txt.gz", code=["IQCT", "LQCT", "NQCT", "OQCT", "PQCT"]),
        # expand("data/tensorqtl/{code}.trans_qtl_pairs.txt.gz", code=["IQCT", "LQCT", "NQCT", "OQCT", "PQCT"]),
        # expand("data/afc/{code}.ASE_aFC.txt", code=["IQCT", "LQCT", "NQCT", "OQCT", "PQCT"]),
        # expand("data/splice/IL/junc/{sample_id}.junc.gz", sample_id=sample_ids["IL"][:10]),
        # expand("data/splice/{tissue}/{tissue}.leafcutter.bed.gz", tissue=["IL"]),
        expand("data/splice/{tissue}/{tissue}_splice.cis_qtl.txt.gz", tissue=["IL", "LHb", "NAcc", "OFC", "PL"]),
        expand("data/splice/{tissue}/{tissue}_splice.cis_independent_qtl.txt.gz", tissue=["IL", "LHb", "NAcc", "OFC", "PL"]),

rule individual_vcf:
    input:
        "data/genotype/P50.rnaseq.88.unpruned.vcf.gz"
    output:
        "data/genotype/individual/{rat_id}.vcf.gz"
    shell:
        "bcftools view -s {wildcards.rat_id} --min-ac=1 -O z -o {output} {input}"

rule observed_snp_list:
    input:
        expand("data/genotype/imputing/chr{chr}.observed.snplist.txt", chr=range(1, 21))
    output:
        "data/genotype/imputing/observed.snplist.txt"
    shell:
        "cat {input} > {output}"

rule vcf_chr_list:
    input:
        vcf = "{base}.vcf.gz",
        vcfi = "{base}.vcf.gz.tbi"
    output:
        "{base}.chrlist.txt"
    shell:
        "tabix --list-chroms {input.vcf} > {output}"

rule index_vcf:
    input:
        "{base}.vcf.gz"
    output:
        "{base}.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"

rule index_bam:
    input:
        "{base}.bam"
    output:
        "{base}.bam.bai"
    conda:
        "../envs/bioinfo.yaml"
    shell:
        "samtools index {input}"

