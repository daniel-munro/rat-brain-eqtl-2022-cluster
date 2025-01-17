localrules:
    qtl2_control_file,
    qtl2_pheno,
    qtl2_tss,
    qtl2_covar,
    qtl2_pvals,

rule qtl2_control_file:
    output:
        "data/qtl2/{code}.yaml"
    conda:
        "envs/qtl2.yaml"
    shell:
        """
        Rscript -e 'qtl2::write_control_file( \
        output_file = "{output}", \
        crosstype = "hs", \
        geno_file = "geno.csv", \
        geno_transposed = TRUE, \
        founder_geno_file = "founder_geno.csv", \
        founder_geno_transposed = TRUE, \
        gmap_file = "gmap.csv", \
        pmap_file = "pmap.csv", \
        pheno_file = "{wildcards.code}.pheno.csv", \
        pheno_transposed = TRUE, \
        covar_file = "{wildcards.code}.covar.csv", \
        sex_covar = "sex", \
        sex_codes = c(F = "female", M = "male"), \
        crossinfo_covar = "generations", \
        geno_codes = c(A = 1L, H = 2L, B = 3L), \
        na.strings = "-")'
        """

rule qtl2_geno_files:
    input:
        pop_vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        pop_vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        founder_vcf = "data/genotype/founders.vcf.gz",
        founder_vcfi = "data/genotype/founders.vcf.gz.tbi",
        snps = "data/genotype/imputing/observed.snplist.txt",
        src = "src/qtl2_geno.py",
    output:
        geno = "data/qtl2/geno.csv",
        founder_geno = "data/qtl2/founder_geno.csv",
        pmap = "data/qtl2/pmap.csv",
        gmap = "data/qtl2/gmap.csv"
    params:
        gen_maps_dir = "data/genotype/genetic_map"
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python3 src/qtl2_geno.py \
        {input.pop_vcf} {input.founder_vcf} {input.snps} {params.gen_maps_dir} \
        {output.geno} {output.founder_geno} {output.pmap} {output.gmap}
        """

rule qtl2_pheno:
    input:
        # "data/tensorqtl/{region}/{region}.main.expression.bed.gz"
        # "data/expression/{region}.rsem_expected_count.bed.gz"
        expr_bed
    output:
        "data/qtl2/{code}.pheno.csv"
    shell:
        # "python3 src/qtl2_pheno.py {input} {output}"
        "zcat {input} | cut -f4- | sed 's/\t/,/g' > {output}"

rule qtl2_tss:
    input:
        # "data/tensorqtl/{region}/{region}.main.expression.bed.gz"
        # "data/expression/{region}.rsem_expected_count.bed.gz"
        expr_bed
    output:
        "data/qtl2/{code}.gene_tss.tsv"
    shell:
        # "zcat {input} | cut -f4,1,2 > {output}"
        """
        echo "gene\tchr\ttss" > {output}
        zcat {input} | tail -n+2 | awk '{{ print $4 "\t" $1 "\t" $2 }}' >> {output}
        """

def qtl2_covar_list(wildcards):
    if wildcards.code[1] in ["C", "D", "E"]:
        return "sex"
    else:
        return "batch,sex"

rule qtl2_covar:
    input:
        "data/DemographicInfo_TN_RNASeq_88_LauraSaba.csv"
    output:
        "data/qtl2/{code}.covar.csv"
    params:
       vars = qtl2_covar_list
    shell:
        "python3 src/qtl2_covar.py {input} {output} {params.vars}"

rule run_qtl2:
    input:
        geno = "data/qtl2/geno.csv",
        founder_geno = "data/qtl2/founder_geno.csv",
        pheno = "data/qtl2/{code}.pheno.csv",
        covar = "data/qtl2/{code}.covar.csv",
        pmap = "data/qtl2/pmap.csv",
        gmap = "data/qtl2/gmap.csv",
        control = "data/qtl2/{code}.yaml",
        tss = "data/qtl2/{code}.gene_tss.tsv"
    output:
        "data/qtl2/{code}.gene_var_lod.tsv.gz"
    conda:
        "envs/qtl2.yaml"
    # threads: 8
    resources:
        walltime = 12
    shell:
        "Rscript src/run_qtl2.R {input.control} {input.tss} {output} nom"

rule run_qtl2_perm:
    input:
        geno = "data/qtl2/geno.csv",
        founder_geno = "data/qtl2/founder_geno.csv",
        pheno = "data/qtl2/{code}.pheno.csv",
        covar = "data/qtl2/{code}.covar.csv",
        pmap = "data/qtl2/pmap.csv",
        gmap = "data/qtl2/gmap.csv",
        control = "data/qtl2/{code}.yaml",
        tss = "data/qtl2/{code}.gene_tss.tsv"
    output:
        "data/qtl2/{code}.gene_lod_perm.tsv.gz"
    conda:
        "envs/qtl2.yaml"
    # threads: 8
    resources:
        walltime = 48,
        cpus = 16
    shell:
        """
        Rscript src/run_qtl2.R {input.control} {input.tss} {output} perm \
        --n_perms 1000 \
        --cores 16
        """

rule run_qtl2_coef:
    input:
        geno = "data/qtl2/geno.csv",
        founder_geno = "data/qtl2/founder_geno.csv",
        pheno = "data/qtl2/{code}.pheno.csv",
        covar = "data/qtl2/{code}.covar.csv",
        pmap = "data/qtl2/pmap.csv",
        gmap = "data/qtl2/gmap.csv",
        control = "data/qtl2/{code}.yaml",
        tss = "data/qtl2/{code}.gene_tss.tsv"
    output:
        "data/qtl2/{code}.gene_var_coefs.tsv.gz"
    conda:
        "envs/qtl2.yaml"
    # threads: 8
    resources:
        walltime = 12
    shell:
        "Rscript src/run_qtl2.R {input.control} {input.tss} {output} coef"

rule qtl2_pvals:
    input:
        pairs = "data/qtl2/{code}.gene_var_lod.tsv.gz",
        perm = "data/qtl2/{code}.gene_lod_perm.tsv.gz",
        coefs = "data/qtl2/{code}.gene_var_coefs.tsv.gz"
    output:
        "data/qtl2/{code}.gene_var_pval.tsv.gz"
    conda:
        "envs/bioinfo.yaml"
    shell:
        "Rscript src/qtl2_pvalues.R {input.pairs} {input.perm} {input.coefs} {output}"
