localrules:
    phaser_expr_matrix,

rule phaser:
    input:
        # bam = bam_file,
        # bai = lambda w: bam_file(w) + ".bai",
        bam = "data/markdup_out/{tissue}/{sample_id}.bam",
        bai = "data/markdup_out/{tissue}/{sample_id}.bam.bai",
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi"
    output:
        "data/phaser_out/{tissue}/{sample_id}.vcf.gz",
        "data/phaser_out/{tissue}/{sample_id}.vcf.gz.tbi",
        "data/phaser_out/{tissue}/{sample_id}.haplotypic_counts.txt"
    params:
        rat_id = lambda w: w.sample_id.split("_")[0]
    conda:
        "../envs/phaser.yaml"
    # group:
    #     "phaser"
    resources:
        cpus = 16,
        walltime = 3
    shell:
        # "python tools/phaser/wrapper.py id {input.bam} {input.vcf} {input.gene_models} output --paired-end"
        """
        python ~/tools/phaser/phaser/phaser.py \
        --temp_dir $TMPDIR \
        --bam {input.bam} \
        --vcf {input.vcf} \
        --sample {params.rat_id} \
        --baseq 10 \
        --mapq 255 \
        --isize 1e6 \
        --paired_end 0 \
        --o data/phaser_out/{wildcards.tissue}/{wildcards.sample_id} \
        --include_indels 0 \
        --gw_phase_vcf 1 \
        --threads 16
        """

rule phaser_gene_ae:
    input:
        counts = "data/phaser_out/{tissue}/{sample_id}.haplotypic_counts.txt",
        gene_models = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.bed"
    output:
        "data/phaser_gene_ae_out/{tissue}/{sample_id}.gene_ae.txt"
    conda:
        "../envs/phaser.yaml"
    # group:
    #     "phaser"
    shell:
        # --min_haplo_maf 0.05 \
        """
        python ~/tools/phaser/phaser_gene_ae/phaser_gene_ae.py \
        --haplotypic_counts {input.counts} \
        --features {input.gene_models} \
        --o {output}
        """


def phaser_gene_ae_files(wildcards):
    return expand("data/phaser_gene_ae_out/{tissue}/{sample_id}.gene_ae.txt",
                  tissue=wildcards.tissue,
                  sample_id=sample_ids[wildcards.tissue])


rule phaser_expr_matrix:
    """Note: Due to phaser bug involving bed column order, internal tabix command doesn't work, but index isn't necessary.
    So I ignore the error and delete the bad index files.
    Also, it must be run from tissue-specific directories, since it uses './tmp/' and deletes all files within at the end.
    """
    input:
        gene_ae = phaser_gene_ae_files,
        gene_models = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.bed"
    output:
        bed = "data/phaser_pop_out/{tissue}/{tissue}.expr_matrix.bed.gz",
        bedgw = "data/phaser_pop_out/{tissue}/{tissue}.expr_matrix.gw_phased.bed.gz",
    params:
        work_dir = "data/phaser_pop_out/{tissue}",
        path_to_smk_dir = "../../..",
        gene_ae_dir = "data/phaser_gene_ae_out/{tissue}",
        prefix = "{tissue}.expr_matrix"
    conda:
        "../envs/phaser.yaml"
    shell:
        """
        cd {params.work_dir}
        python ~/tools/phaser/phaser_pop/phaser_expr_matrix.py \
            --gene_ae_dir {params.path_to_smk_dir}/{params.gene_ae_dir} \
            --features {params.path_to_smk_dir}/{input.gene_models} \
            --o {params.prefix} || true
        cd {params.path_to_smk_dir}
        rm {output.bed}.tbi {output.bedgw}.tbi
        """
