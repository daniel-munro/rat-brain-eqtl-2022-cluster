rule aFC_from_eQTL_model_main:
    input:
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        # bed = "data/tensorqtl/{region}/{region}.main.expression.bed.gz",
        # bedi = "data/tensorqtl/{region}/{region}.main.expression.bed.gz.tbi",
        bed = "data/expression/{region}.rsem_expected_count.bed.gz",
        bedi = "data/expression/{region}.rsem_expected_count.bed.gz.tbi",
        qtl = "data/tensorqtl/{region}/{region}.main.cis_qtl.txt.gz",
        covar = "data/tensorqtl/{region}/{region}.main.combined_covariates.txt"
    output:
        "data/afc/{region}.main.aFC.txt"
    conda:
        "envs/afc.yaml"
    shell:
        """
        python3 tools/aFC/aFC.py \
        --vcf {input.vcf} \
        --pheno {input.bed} \
        --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl}) \
        --cov <(sed 's/-N//g' {input.covar}) \
        --log_xform 1 \
        --output {output}
        """

rule aFC_from_eQTL_model_main2:
    input:
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        bed = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz",
        bedi = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz.tbi",
        qtl = "data/tensorqtl/{region}/{region}.main2.cis_qtl.txt.gz",
        covar = "data/tensorqtl/{region}/{region}.main2.combined_covariates.txt"
    output:
        "data/afc/{region}.main2.aFC.txt"
    conda:
        "envs/afc.yaml"
    shell:
        """
        python3 tools/aFC/aFC.py \
        --vcf {input.vcf} \
        --pheno {input.bed} \
        --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl}) \
        --cov <(sed 's/-N//g' {input.covar}) \
        --log_xform 1 \
        --output {output}
        """

rule aFC_from_eQTL_model_basic:
    input:
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        # bed = "data/tensorqtl/{region}/{region}.main.expression.bed.gz",
        # bedi = "data/tensorqtl/{region}/{region}.main.expression.bed.gz.tbi",
        bed = "data/expression/{region}.rsem_expected_count.bed.gz",
        bedi = "data/expression/{region}.rsem_expected_count.bed.gz.tbi",
        qtl = "data/tensorqtl/{region}/{region}.basic.cis_qtl.txt.gz"
    output:
        "data/afc/{region}.basic.aFC.txt"
    conda:
        "envs/afc.yaml"
    shell:
        """
        python3 tools/aFC/aFC.py \
        --vcf {input.vcf} \
        --pheno {input.bed} \
        --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl}) \
        --log_xform 1 \
        --output {output}
        """

rule aFC_from_eQTL_model_basic2:
    input:
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        bed = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz",
        bedi = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz.tbi",
        qtl = "data/tensorqtl/{region}/{region}.basic2.cis_qtl.txt.gz"
    output:
        "data/afc/{region}.basic2.aFC.txt"
    conda:
        "envs/afc.yaml"
    shell:
        """
        python3 tools/aFC/aFC.py \
        --vcf {input.vcf} \
        --pheno {input.bed} \
        --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl}) \
        --log_xform 1 \
        --output {output}
        """

rule aFC_from_eQTL_model_qtl2:
    input:
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        bed = "data/tensorqtl/{region}/{region}.main.expression.bed.gz",
        bedi = "data/tensorqtl/{region}/{region}.main.expression.bed.gz.tbi",
        qtl = "data/qtl2/{region}.gene_var_pval.tsv.gz",
        covar = "data/tensorqtl/other_covariates.tsv"
    output:
        "data/afc/{region}.qtl2.aFC.txt"
    conda:
        "envs/afc.yaml"    
    shell:
        """
        python3 tools/aFC/aFC.py \
        --vcf {input.vcf} \
        --pheno {input.bed} \
        --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl} 1 2) \
        --cov {input.covar} \
        --log_xform 1 \
        --output {output}
        """

###################
## ASE-based aFC ##
###################

rule phaser_sample_map:
    output:
        "data/afc/{region}.sample_map.txt"
    run:
        with open(output[0], "w") as out:
            out.write("vcf_sample\tbed_sample\n")
            for sample_id in sample_ids[wildcards.region]:
                rat_id = sample_id.split("_")[0]
                bam_base = "{}.Aligned.sortedByCoord.out".format(sample_id)
                out.write("{}\t{}\n".format(rat_id, bam_base))

rule phaser_gene_var_pairs:
    input:
        "data/afc/{region}.{method}.aFC.txt"
    output:
        "data/afc/{region}.{method}.pairs.txt"
    run:
        inlines = open(input[0], "r").read().splitlines()
        gene_var = [x.split("\t")[:2] for x in inlines[1:]]
        with open(output[0], "w") as out:
            out.write("gene_id\tvar_id\tvar_contig\tvar_pos\tvar_ref\tvar_alt\n")
            for x in gene_var:
                chrom, pos = tuple(x[1].split(":"))
                chrom = chrom.replace("chr", "")
                out.write("{}\t{}\t{}\t{}\t\t\n".format(x[0], x[1], chrom, pos))

rule phaser_expr_matrix:
    input:
        SFTP.remote("garibaldi.scripps.edu/gpfs/home/dmunro/data/phaser_pop_out/{region}.expr_matrix.gw_phased.bed.gz")
    output:
        "data/expression/{region}.expr_matrix.gw_phased.bed.gz"
    shell:
        "cp {input} {output}"

# rule bed_index:
#     input:
#         "{base}.bed.gz"
#     output:
#         "{base}.bed.gz.tbi"
#     shell:
#         "tabix {input}"

rule phaser_cis_var:
    input:
        bed = "data/expression/{region}.expr_matrix.gw_phased.bed.gz",
        # bedi = "data/expression/{region}.expr_matrix.gw_phased.bed.gz.tbi",
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        gene_var_pairs = "data/afc/{region}.{method}.pairs.txt",
        sample_map = "data/afc/{region}.sample_map.txt"
    output:
        "data/afc/{region}.{method}.ASE_aFC.txt"
    conda:
        "envs/phaser.yaml"
    threads: 8
    shell:
        """
        python tools/phaser/phaser_pop/phaser_cis_var.py \
        --bed {input.bed} \
        --vcf {input.vcf} \
        --pair {input.gene_var_pairs} \
        --map {input.sample_map} \
        --t 8 \
        --o {output}
        """

# rule aFC_from_eQTL_model_main:
#     input:
#         # vcf = "data/genotype/round8_unpruned.vcf.gz",
#         # vcfi = "data/genotype/round8_unpruned.vcf.gz.tbi",
#         vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
#         vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
#         bed = "data/fastqtl/main.expression.bed.gz",
#         bedi = "data/fastqtl/main.expression.bed.gz.tbi",
#         qtl = "data/fastqtl/main.genes.txt.gz",
#         covar = "data/fastqtl/main.combined_covariates.txt"
#     output:
#         "data/fastqtl/main.aFC.txt"
#     shell:
#         """
#         python3 tools/aFC/aFC.py \
#         --vcf {input.vcf} \
#         --pheno {input.bed} \
#         --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl}) \
#         --cov {input.covar} \
#         --log_xform 1 \
#         --output {output}
#         """

# rule aFC_from_eQTL_model_basic:
#     input:
#         # vcf = "data/genotype/round8_unpruned.vcf.gz",
#         # vcfi = "data/genotype/round8_unpruned.vcf.gz.tbi",
#         vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
#         vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
#         bed = "data/fastqtl/main.expression.bed.gz",
#         bedi = "data/fastqtl/main.expression.bed.gz.tbi",
#         qtl = "data/fastqtl/basic.genes.txt.gz"
#     output:
#         "data/fastqtl/basic.aFC.txt"
#     shell:
#         """
#         python3 tools/aFC/aFC.py \
#         --vcf {input.vcf} \
#         --pheno {input.bed} \
#         --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl}) \
#         --log_xform 1 \
#         --output {output}
#         """

# rule aFC_from_eQTL_model_qtl2:
#     input:
#         # vcf = "data/genotype/round8_unpruned.vcf.gz",
#         # vcfi = "data/genotype/round8_unpruned.vcf.gz.tbi",
#         vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
#         vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
#         bed = "data/fastqtl/main.expression.bed.gz",
#         bedi = "data/fastqtl/main.expression.bed.gz.tbi",
#         qtl = "data/qtl2/Acbc_gene_var_pval.tsv.gz"
#     output:
#         "data/qtl2/qtl2.aFC.txt"
#     shell:
#         """
#         python3 tools/aFC/aFC.py \
#         --vcf {input.vcf} \
#         --pheno {input.bed} \
#         --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl} 1 2) \
#         --log_xform 1 \
#         --output {output}
#         """
