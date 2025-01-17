localrules:
    vcf_and_bed_to_bimbam,
    gemma_make_covar,
    # gemma_top_eqtls,
    gemma_pve,

# rule vcf_to_plink:
#     input:
#         "data/genotype/P50.rnaseq.88.unpruned.vcf.gz"
#     output:
#         multiext("data/gemma/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam")
#     params:
#         prefix = "data/gemma/P50.rnaseq.88.unpruned"
#     shell:
#         """
#         plink2 --make-bed \
#         --vcf {input} \
#         --out {params.prefix}
#         """

rule vcf_and_bed_to_bimbam:
    input:
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        bed = "data/expression/ensembl-gene_inv-quant_ComBat_Acbc.bed.gz"
    output:
        geno = "data/gemma/geno.txt.gz",
        pheno = "data/gemma/pheno.txt",
        genes = "data/gemma/genes.txt",
        snps = "data/gemma/snps.txt",
        samples = "data/gemma/samples.txt", # Just for reference
    params:
        geno = "data/gemma/geno.txt"
    conda:
        "envs/biopython.yaml"
    shell:
        # plink2 \
        # --vcf data/genotype/P50.rnaseq.88.unpruned.vcf.gz \
        # --recode bimbam \
        # --out data/gemma/P50.rnaseq.88.unpruned
        """
        python3 src/vcf_and_bed_to_bimbam.py \
            {input.vcf} {input.bed} \
            {params.geno} {output.pheno} {output.genes} {output.snps} \
            {output.samples}
        gzip {params.geno}
        """

rule gemma_make_covar:
    input:
        covar = "data/tensorqtl/AQCT.combined_covariates.txt",
        samples = "data/gemma/samples.txt",
    output:
        covar = "data/gemma/covar.txt",
        names = "data/gemma/covar_names.txt",
    run:
        samples = list(pd.read_csv(input.samples, header=None)[0])
        d = pd.read_csv(input.covar, sep="\t", dtype=str)
        d = d.rename(columns={"ID": "intercept"}) # just to have in names file
        d["intercept"].to_csv(output.names, index=False)
        # d = d[["intercept"] + samples].T
        d = d[samples].T
        d.insert(0, "intercept", 1)
        d.to_csv(output.covar, sep="\t", index=False, header=False)

rule gemma_kinship_loco:
    input:
        geno = "data/gemma/geno.txt.gz",
        pheno = "data/gemma/pheno.txt",
    output:
        matrix = "data/gemma/kinship/{chrn}.cXX.txt",
        log = "data/gemma/kinship/{chrn}.log.txt",
    params:
        tmpgeno = "output/geno_loco.{chrn}.txt"
    shell:
        """
        zcat {input.geno} | grep -v "^chr{wildcards.chrn}:" > {params.tmpgeno}
        gemma -gk \
            -g {params.tmpgeno} \
            -p {input.pheno} \
            -o {wildcards.chrn}
        mv output/{wildcards.chrn}.cXX.txt {output.matrix}
        mv output/{wildcards.chrn}.log.txt {output.log}
        rm {params.tmpgeno}
        """

# rule gemma:
#     input:
#         geno = "data/gemma/geno.txt",
#         pheno = "data/gemma/pheno.txt",
#         snps = "data/gemma/snps.txt",
#     output:
#         assoc = "data/gemma/output/Acbc.assoc.txt",
#         log = "data/gemma/output/Acbc.log.txt",
#     params:
#         path = "data/gemma",
#         cdgeno = "geno.txt",
#         cdpheno = "pheno.txt",
#         cdsnps = "snps.txt",
#         prefix = "Acbc",
#     shell:
#         """
#         cd {params.path}
#         gemma -g {params.cdgeno} -p {params.cdpheno} -a {params.cdsnps} \
#             -lm 1 -o {params.prefix}
#         """

rule gemma_eqtl_lm:
    input:
        geno = "data/gemma/geno.txt.gz",
        pheno = "data/gemma/pheno.txt",
        snps = "data/gemma/snps.txt",
        genes = "data/gemma/genes.txt",
    output:
        assoc = "data/gemma/lm/Acbc.{chrn}.assoc.txt",
        log = "data/gemma/lm/Acbc.{chrn}.log.txt",
    resources:
        walltime = 12,
        cpus = 16
    shell:
        """
        python3 src/gemma_eqtls.py \
            {input.geno} \
            {input.pheno} \
            {input.snps} \
            {input.genes} \
            {output.assoc} \
            --log {output.log} \
            --chr {wildcards.chrn} \
            --jobs 16
        """

rule gemma_eqtl_lmm:
    input:
        geno = "data/gemma/geno.txt.gz",
        pheno = "data/gemma/pheno.txt",
        snps = "data/gemma/snps.txt",
        genes = "data/gemma/genes.txt",
        kinship = "data/gemma/kinship/{chrn}.cXX.txt"
    output:
        assoc = "data/gemma/lmm/Acbc.{chrn}.assoc.txt",
        log = "data/gemma/lmm/Acbc.{chrn}.log.txt",
    resources:
        walltime = 12,
        mem_mb = 20000,
        cpus = 16
    shell:
        """
        python3 src/gemma_eqtls.py \
            {input.geno} \
            {input.pheno} \
            {input.snps} \
            {input.genes} \
            {output.assoc} \
            --log {output.log} \
            --chr {wildcards.chrn} \
            --kinship {input.kinship} \
            --jobs 16
        """

rule gemma_eqtl_lm_rand:
    input:
        geno = "data/gemma/geno.txt.gz",
        pheno = "data/gemma/pheno_rand.txt",
        snps = "data/gemma/snps.txt",
        genes = "data/gemma/genes.txt",
    output:
        assoc = "data/gemma/lm_rand/Acbc.{chrn}.assoc.txt",
        log = "data/gemma/lm_rand/Acbc.{chrn}.log.txt",
    resources:
        walltime = 12,
        cpus = 16
    shell:
        """
        python3 src/gemma_eqtls.py \
            {input.geno} \
            {input.pheno} \
            {input.snps} \
            {input.genes} \
            {output.assoc} \
            --log {output.log} \
            --chr {wildcards.chrn} \
            --jobs 16
        """

rule gemma_eqtl_lmm_rand:
    input:
        geno = "data/gemma/geno.txt.gz",
        pheno = "data/gemma/pheno_rand.txt",
        snps = "data/gemma/snps.txt",
        genes = "data/gemma/genes.txt",
        kinship = "data/gemma/kinship/{chrn}.cXX.txt"
    output:
        assoc = "data/gemma/lmm_rand/Acbc.{chrn}.assoc.txt",
        log = "data/gemma/lmm_rand/Acbc.{chrn}.log.txt",
    resources:
        walltime = 12,
        cpus = 16
    shell:
        """
        python3 src/gemma_eqtls.py \
            {input.geno} \
            {input.pheno} \
            {input.snps} \
            {input.genes} \
            {output.assoc} \
            --log {output.log} \
            --chr {wildcards.chrn} \
            --kinship {input.kinship} \
            --jobs 16
        """

rule gemma_eqtl_lm_cov:
    input:
        geno = "data/gemma/geno.txt.gz",
        pheno = "data/gemma/pheno.txt",
        covar = "data/gemma/covar.txt",
        snps = "data/gemma/snps.txt",
        genes = "data/gemma/genes.txt",
    output:
        assoc = "data/gemma/lm_cov/Acbc.{chrn}.assoc.txt",
        log = "data/gemma/lm_cov/Acbc.{chrn}.log.txt",
    resources:
        walltime = 12,
        cpus = 16
    shell:
        """
        python3 src/gemma_eqtls.py \
            {input.geno} \
            {input.pheno} \
            {input.snps} \
            {input.genes} \
            {output.assoc} \
            --covar {input.covar} \
            --log {output.log} \
            --chr {wildcards.chrn} \
            --jobs 16
        """

rule gemma_eqtl_lmm_cov:
    input:
        geno = "data/gemma/geno.txt.gz",
        pheno = "data/gemma/pheno.txt",
        covar = "data/gemma/covar.txt",
        snps = "data/gemma/snps.txt",
        genes = "data/gemma/genes.txt",
        kinship = "data/gemma/kinship/{chrn}.cXX.txt"
    output:
        assoc = "data/gemma/lmm_cov/Acbc.{chrn}.assoc.txt",
        log = "data/gemma/lmm_cov/Acbc.{chrn}.log.txt",
    resources:
        walltime = 12,
        cpus = 16
    shell:
        """
        python3 src/gemma_eqtls.py \
            {input.geno} \
            {input.pheno} \
            {input.snps} \
            {input.genes} \
            {output.assoc} \
            --covar {input.covar} \
            --log {output.log} \
            --chr {wildcards.chrn} \
            --kinship {input.kinship} \
            --jobs 16
        """

rule gemma_top_eqtls:
    input:
        expand("data/gemma/{{mode}}/Acbc.{chrn}.assoc.txt", chrn=range(1, 21))
    output:
        "data/gemma/Acbc.{mode}.assoc.txt"
    shell:
        "python3 src/gemma_top_eqtls.py {input} -p p_wald -o {output}"

rule gemma_random_pairs:
    input:
        expand("data/gemma/{mode}/Acbc.{chrn}.assoc.txt",
               mode=["lm", "lmm"],
               chrn=range(1, 21))
    output:
        "data/gemma/Acbc.random_pairs.txt",
    shell:
        "python3 src/gemma_random_eqtls.py"

rule gemma_pve:
    input:
        expand("data/gemma/{{version}}/Acbc.{chrn}.log.txt", chrn=range(1, 21))
    output:
        "data/gemma/Acbc.{version}.pve.txt"
    shell:
        "python3 src/gemma_parse_lmm_logs.py {input} -o {output}"

rule eigenMT:
    input:
        ""
    output:
        ""
    shell:
        """
        python eigenMT.py --chrom {wildcards.chrn} \
        --QTL /gemma_loco/output_loco/OF/OF_gemma_chr${i}_results.txt \
        --GEN OF_genotype_chr${i}.txt \
        --GENPOS chr${i}_gen_pos.txt \
        --PHEPOS rn6_gene_pos.tab \
        --OUT chr${i}_eigenMT_win50_result.txt
        """

# rule phenotype_file_for_plink:
#     input:
#         "data/expression/ensembl-gene_inv-quant_ComBat_Acbc.bed.gz"
#     output:
#         "data/gemma/pheno.txt"
#     run:
#         d = pd.read_csv(input, sep="\t")
