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
        # "data/rsem_expected_count.gct.gz"
        "data/tensorqtl/{region}/{region}.main.expression.bed.gz"
    output:
        "data/qtl2/{region}.pheno.csv"
    shell:
        # "python3 src/qtl2_pheno.py {input} {output}"
        "zcat {input} | cut -f4- | sed 's/\t/,/g' > {output}"

rule qtl2_tss:
    input:
        "data/tensorqtl/{region}/{region}.main.expression.bed.gz"
    output:
        "data/qtl2/{region}.gene_tss.tsv"
    shell:
        # "zcat {input} | cut -f4,1,2 > {output}"
        """
        echo "gene\tchr\ttss" > {output}
        zcat {input} | tail -n+2 | awk '{{ print $4 "\t" $1 "\t" $2 }}' >> {output}
        """

rule qtl2_covar:
    input:
        "data/DemographicInfo_TN_RNASeq_88_LauraSaba.csv"
    output:
        "data/qtl2/covar.csv"
    shell:
        "python3 src/qtl2_covar.py {input} {output}"

rule run_qtl2:
    input:
        geno = "data/qtl2/geno.csv",
        founder_geno = "data/qtl2/founder_geno.csv",
        pheno = "data/qtl2/{region}.pheno.csv",
        covar = "data/qtl2/covar.csv",
        pmap = "data/qtl2/pmap.csv",
        gmap = "data/qtl2/gmap.csv",
        control = "data/qtl2/{region}.yaml",
        tss = "data/qtl2/{region}.gene_tss.tsv"
    output:
        "data/qtl2/{region}.gene_var_lod.tsv.gz"
    conda:
        "envs/qtl2.yaml"
    threads: 8
    shell:
        "Rscript src/run_qtl2.R {input.control} {input.tss} {output} nom"

rule run_qtl2_perm:
    input:
        geno = "data/qtl2/geno.csv",
        founder_geno = "data/qtl2/founder_geno.csv",
        pheno = "data/qtl2/{region}.pheno.csv",
        covar = "data/qtl2/covar.csv",
        pmap = "data/qtl2/pmap.csv",
        gmap = "data/qtl2/gmap.csv",
        control = "data/qtl2/{region}.yaml",
        tss = "data/qtl2/{region}.gene_tss.tsv"
    output:
        "data/qtl2/{region}.gene_lod_perm.tsv.gz"
    conda:
        "envs/qtl2.yaml"
    threads: 8
    shell:
        "Rscript src/run_qtl2.R {input.control} {input.tss} {output} perm --n_perms 200"

rule run_qtl2_coef:
    input:
        geno = "data/qtl2/geno.csv",
        founder_geno = "data/qtl2/founder_geno.csv",
        pheno = "data/qtl2/{region}.pheno.csv",
        covar = "data/qtl2/covar.csv",
        pmap = "data/qtl2/pmap.csv",
        gmap = "data/qtl2/gmap.csv",
        control = "data/qtl2/{region}.yaml",
        tss = "data/qtl2/{region}.gene_tss.tsv"
    output:
        "data/qtl2/{region}.gene_var_coefs.tsv.gz"
    conda:
        "envs/qtl2.yaml"
    threads: 8
    shell:
        "Rscript src/run_qtl2.R {input.control} {input.tss} {output} coef"

rule qtl2_pvals:
    input:
        pairs = "data/qtl2/{region}.gene_var_lod.tsv.gz",
        perm = "data/qtl2/{region}.gene_lod_perm.tsv.gz",
        coefs = "data/qtl2/{region}.gene_var_coefs.tsv.gz"
    output:
        "data/qtl2/{region}.gene_var_pval.tsv.gz"
    conda:
        "envs/bioinfo.yaml"
    shell:
        "Rscript src/qtl2_pvalues.R {input.pairs} {input.perm} {input.coefs} {output}"
