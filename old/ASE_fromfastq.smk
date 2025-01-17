rule phaser:
    input:
        bam = "data/star_out/{region}/{rat_id}_{region}.Aligned.sortedByCoord.out.bam",
        bai = "data/star_out/{region}/{rat_id}_{region}.Aligned.sortedByCoord.out.bam.bai",
        vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
        vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi"
    output:
        "data/phaser_out/{region}/{rat_id}_{region}.vcf.gz",
        "data/phaser_out/{region}/{rat_id}_{region}.vcf.gz.tbi",
        "data/phaser_out/{region}/{rat_id}_{region}.haplotypic_counts.txt"
    conda:
        "envs/phaser.yaml"
    group:
        "phaser"
    shell:
        # "python tools/phaser/wrapper.py id {input.bam} {input.vcf} {input.gene_models} output --paired-end"
        """
        python tools/phaser/phaser/phaser.py \
        --temp_dir $PBSTMPDIR \
        --bam {input.bam} \
        --vcf {input.vcf} \
        --sample {wildcards.rat_id} \
        --baseq 10 \
        --mapq 255 \
        --isize 1e6 \
        --paired_end 0 \
        --o data/phaser_out/{wildcards.region}/{wildcards.rat_id}_{wildcards.region} \
        --include_indels 0 \
        --gw_phase_vcf 1
        """

rule phaser_gene_ae:
    input:
        counts = "data/phaser_out/{region}/{rat_id}_{region}.haplotypic_counts.txt",
        gene_models = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.bed"
    output:
        "data/phaser_gene_ae_out/{region}/{rat_id}_{region}.gene_ae.txt"
    conda:
        "envs/phaser.yaml"
    group:
        "phaser"
    shell:
        # --min_haplo_maf 0.05 \
        """
        python tools/phaser/phaser_gene_ae/phaser_gene_ae.py \
        --haplotypic_counts {input.counts} \
        --features {input.gene_models} \
        --o {output}
        """


def phaser_gene_ae_files(wildcards):
    return expand("data/phaser_gene_ae_out/{region}/{sample_id}.gene_ae.txt",
                  region=wildcards.region,
                  sample_id=sample_ids[wildcards.region])


rule phaser_expr_matrix:
    """Note: Due to phaser bug involving bed column order, internal tabix command doesn't work, but index isn't necessary. So I ignore the error and delete the bad index files."""
    input:
        gene_ae = phaser_gene_ae_files,
        gene_models = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.bed"
    output:
        bed = "data/phaser_pop_out/{region}.expr_matrix.bed.gz",
        # bedi = "data/phaser_pop_out/{region}.expr_matrix.bed.gz.tbi",
        bedgw = "data/phaser_pop_out/{region}.expr_matrix.gw_phased.bed.gz",
        # bedgwi = "data/phaser_pop_out/{region}.expr_matrix.gw_phased.bed.gz.tbi"
    conda:
        "envs/phaser.yaml"
    shell:
        # tabix {output.bed}
        # tabix {output.bedgw}
        """
        python tools/phaser/phaser_pop/phaser_expr_matrix.py \
        --gene_ae_dir data/phaser_gene_ae_out/{wildcards.region} \
        --features {input.gene_models} \
        --o data/phaser_pop_out/{wildcards.region}.expr_matrix || true
        rm {output.bed}.tbi {output.bedgw}.tbi
        """
