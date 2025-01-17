# rule plink_to_vcf:
#     input:
#         plink = expand("data/genotype/{{base}}.{ext}", ext=["bed", "bim", "fam"]),
#         ids = "rat_ids.txt"
#     params:
#         prefix = "data/genotype/{base}"
#     output:
#         "data/genotype/{base}.vcf"
#     shell:
#         """
#         tools/plink2 --bfile {params.prefix} \
#         --recode vcf id-paste=iid \
#         --keep {input.ids} \
#         --out {params.prefix}
#         """

rule individual_vcf:
    input:
        # "data/genotype/round8_unpruned.vcf.gz"
        "data/genotype/phased/P50.rnaseq.88.unpruned.vcf.gz"
    output:
        "data/vcfs_for_star/{rat_id}.vcf.gz"
    shell:
        "bcftools view -s {wildcards.rat_id} --min-ac=1 -O z -o {output} {input}"

# rule bgzip_vcf:
#     input: "data/genotype/{base}.vcf"
#     output: "data/genotype/{base}.vcf.gz"
#     shell:
#         "bgzip --threads 8 {input}"

# rule qtl2_genotype:
#     input:
#         # vcf = "data/genotype/round8_unpruned.vcf.gz"
#         vcf = "data/genotype/phased/P50.rnaseq.88.unpruned.vcf.gz",
#         snps = "data/genotype/imputing/observed.snplist.txt"
#     output:
#         "data/qtl2/geno.csv"
#     shell:
#         # This only outputs ACGT alleles, which aren't allowed.
#         # "gatk VariantsToTable -V {input} -F ID -GF GT -O {output}"
#         "python3 src/qtl2_geno.py {input.vcf} {output} {input.snps}"

# rule qtl2_founder_genotype:
#     input:
#         vcf = "data/genotype/founders.vcf.gz",
#         # geno = "data/qtl2/geno.csv"
#         snps = "data/genotype/imputing/observed.snplist.txt"
#     output:
#         "data/qtl2/founder_geno.csv"
#     shell:
#         # This only outputs ACGT alleles, which aren't allowed.
#         # "gatk VariantsToTable -V {input} -F ID -GF GT -O {output}"
#         """
#         python3 src/qtl2_geno.py {input.vcf} {output} {input.snps}
#         """

# rule qtl2_pmap:
#     input:
#         # "data/genotype/round8_unpruned.vcf.gz"
#         "data/genotype/imputing/observed.snplist.txt"
#     output:
#         "data/qtl2/pmap.csv"
#     shell:
#         "python3 src/qtl2_pmap.py {input} {output}"

# rule qtl2_gmap:
#     input:
#         # "data/genotype/round8_unpruned.vcf.gz"
#         "data/genotype/imputing/observed.snplist.txt"
#     output:
#         "data/qtl2/gmap.csv"
#     shell:
#         "python3 src/qtl2_gmap.py {input} {output} data/qtl2/genetic_map"


rule run_qtl2_error_prob:
    input:
        geno = "data/qtl2/geno.csv",
        founder_geno = "data/qtl2/founder_geno.csv",
        pheno = "data/qtl2/{region}_pheno.csv",
        covar = "data/qtl2/covar.csv",
        pmap = "data/qtl2/pmap.csv",
        gmap = "data/qtl2/gmap.csv",
        control = "data/qtl2/{region}.yaml",
        tss = "data/qtl2/gene_tss.tsv"
    output:
        "data/qtl2/error_{error}/{region}_gene_var_lod.tsv.gz"
    threads: 8
    shell:
        """
        Rscript src/run_qtl2.R {input.control} {input.tss} {output} nom \
        --error_prob {wildcards.error}
        """

rule run_qtl2_perm_error_prob:
    input:
        geno = "data/qtl2/geno.csv",
        founder_geno = "data/qtl2/founder_geno.csv",
        pheno = "data/qtl2/{region}_pheno.csv",
        covar = "data/qtl2/covar.csv",
        pmap = "data/qtl2/pmap.csv",
        gmap = "data/qtl2/gmap.csv",
        control = "data/qtl2/{region}.yaml",
        tss = "data/qtl2/gene_tss.tsv"
    output:
        "data/qtl2/error_{error}/{region}_gene_lod_perm.tsv.gz"
    threads: 8
    shell:
        """
        Rscript src/run_qtl2.R {input.control} {input.tss} {output} perm \
        --error_prob {wildcards.error}
        """

rule run_qtl2_coef_error_prob:
    input:
        geno = "data/qtl2/geno.csv",
        founder_geno = "data/qtl2/founder_geno.csv",
        pheno = "data/qtl2/{region}_pheno.csv",
        covar = "data/qtl2/covar.csv",
        pmap = "data/qtl2/pmap.csv",
        gmap = "data/qtl2/gmap.csv",
        control = "data/qtl2/{region}.yaml",
        tss = "data/qtl2/gene_tss.tsv"
    output:
        "data/qtl2/error_{error}/{region}_gene_var_coefs.tsv.gz"
    threads: 8
    shell:
        """
        Rscript src/run_qtl2.R {input.control} {input.tss} {output} coef \
        --error_prob {wildcards.error}
        """

rule qtl2_pvals_error_prob:
    input:
        pairs = "data/qtl2/error_{error}/{region}_gene_var_lod.tsv.gz",
        perm = "data/qtl2/error_{error}/{region}_gene_lod_perm.tsv.gz",
        coefs = "data/qtl2/error_{error}/{region}_gene_var_coefs.tsv.gz"
    output:
        "data/qtl2/error_{error}/{region}_gene_var_pval.tsv"
    shell:
        "Rscript src/qtl2_pvalues.R {input.pairs} {input.perm} {input.coefs} {output}"
