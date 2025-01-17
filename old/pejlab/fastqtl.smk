rule sample_participant_lookup:
    input: "rat_ids.txt"
    output: "data/samples_participants.txt"
    run:
        with open(output[0], "w") as out:
            out.write("sample_id\tparticipant_id\n")
            for sample_id in sample_ids:
                rat_id = sample_id.split("_")[0]
                out.write("{0}\t{1}\n".format(sample_id, rat_id))

rule prepare_expression:
    input:
        # tpm_gct = "data/gene_tpm.gct.gz",
        # counts_gct = "data/gene_reads.gct.gz",
        tpm_gct = "data/rsem_TPM.gct.gz",
        counts_gct = "data/expression/rsem_expected_count.gct.gz",
        gtf = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf",
        samples = "data/samples_participants.txt",
        # chrlist = "data/genotype/round8_unpruned.chrlist.txt"
        chrlist = "data/genotype/phased/P50.rnaseq.88.unpruned.chrlist.txt"
    output:
        "data/fastqtl/main.expression.bed.gz",
        "data/fastqtl/main.expression.bed.gz.tbi"
    shell:
        """
        python3 src/eqtl_prepare_expression.py \
        {input.tpm_gct} {input.counts_gct} {input.gtf} \
        {input.samples} {input.chrlist} data/fastqtl/main \
        --tpm_threshold 0.1 \
        --count_threshold 6 \
        --sample_frac_threshold 0.2 \
        --normalization_method tmm
        """

rule expression_tsv_to_bed:
    input:
        tsv = "data/expression/ensemblGeneV2_Norm-counts-log2+1.txt",
        gtf = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf"
    output:
        bed = "data/expression/ensemblGeneV2_Norm-counts-log2+1.bed.gz",
        bedi = "data/expression/ensemblGeneV2_Norm-counts-log2+1.bed.gz.tbi"
    params:
        interm = "data/expression/ensemblGeneV2_Norm-counts-log2+1.bed",
        region = "Acbc"
    conda:
        "envs/bioinfo.yaml"
    shell:
        """
        python3 src/tsv_to_bed.py {input.tsv} {input.gtf} {params.interm} --region {params.region}
        bgzip {params.interm}
        tabix {output.bed}
        """

rule run_peer:
    input:
        "data/fastqtl/main.expression.bed.gz"
    output:
        "data/fastqtl/main.PEER_residuals.txt",
        "data/fastqtl/main.PEER_alpha.txt",
        "data/fastqtl/main.PEER_covariates.txt"
    params:
        num_peer = 15
    conda:
        "envs/peer.yaml"
    shell:
        "Rscript src/run_PEER.R {input} data/fastqtl/main {params.num_peer}"

rule other_covariates:
    input:
        "data/DemographicInfo_TN_RNASeq_88_LauraSaba.csv"
    output:
        "data/fastqtl/other_covariates.tsv"
    shell:
        "python3 src/fastqtl_covar.py {input} {output}"

rule combine_covariates:
    input:
        peer = "data/fastqtl/main.PEER_covariates.txt",
        # pcs = "data/genotype_pcs/genotype_pcs.tsv",
        other = "data/fastqtl/other_covariates.tsv"
    output:
        "data/fastqtl/main.combined_covariates.txt"
    shell:
    #   --genotype_pcs {input.pcs} \
        """
        python3 src/combine_covariates.py {input.peer} data/fastqtl/main \
        --add_covariates {input.other}
        """

rule run_fastqtl_nom:
    input:
        # vcf = "data/genotype/round8_unpruned.vcf.gz",
        vcf = "data/genotype/phased/P50.rnaseq.88.unpruned.vcf.gz",
        vcfi = "data/genotype/phased/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        bed = "data/fastqtl/main.expression.bed.gz",
        bedi = "data/fastqtl/main.expression.bed.gz.tbi",
        covar = "data/fastqtl/main.combined_covariates.txt"
    output:
        "data/fastqtl/main.allpairs.txt.gz",
    shell:
        """
        tools/fastqtl/python/run_FastQTL_threaded.py \
        {input.vcf} {input.bed} data/fastqtl/main \
        --covariates {input.covar}
        """

rule run_fastqtl_perm:
    input:
        # vcf = "data/genotype/round8_unpruned.vcf.gz",
        vcf = "data/genotype/phased/P50.rnaseq.88.unpruned.vcf.gz",
        vcfi = "data/genotype/phased/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        bed = "data/fastqtl/main.expression.bed.gz",
        bedi = "data/fastqtl/main.expression.bed.gz.tbi",
        covar = "data/fastqtl/main.combined_covariates.txt"
    output:
        "data/fastqtl/main.genes.txt.gz"
    shell:
        """
        tools/fastqtl/python/run_FastQTL_threaded.py \
        {input.vcf} {input.bed} data/fastqtl/main \
        --covariates {input.covar} \
        --permute 1000 10000
        """

rule run_fastqtl_nom_basic:
    input:
        # vcf = "data/genotype/round8_unpruned.vcf.gz",
        vcf = "data/genotype/phased/P50.rnaseq.88.unpruned.vcf.gz",
        vcfi = "data/genotype/phased/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        bed = "data/fastqtl/main.expression.bed.gz",
        bedi = "data/fastqtl/main.expression.bed.gz.tbi"
    output:
        "data/fastqtl/basic.allpairs.txt.gz",
    shell:
        """
        tools/fastqtl/python/run_FastQTL_threaded.py \
        {input.vcf} {input.bed} data/fastqtl/basic
        """

rule run_fastqtl_perm_basic:
    input:
        # vcf = "data/genotype/round8_unpruned.vcf.gz",
        vcf = "data/genotype/phased/P50.rnaseq.88.unpruned.vcf.gz",
        vcfi = "data/genotype/phased/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        bed = "data/fastqtl/main.expression.bed.gz",
        bedi = "data/fastqtl/main.expression.bed.gz.tbi"
    output:
        "data/fastqtl/basic.genes.txt.gz"
    shell:
        """
        tools/fastqtl/python/run_FastQTL_threaded.py \
        {input.vcf} {input.bed} data/fastqtl/basic \
        --permute 1000 10000
        """

rule run_fastqtl_perm_basic2:
    input:
        vcf = "data/genotype/phased/P50.rnaseq.88.unpruned.vcf.gz",
        vcfi = "data/genotype/phased/P50.rnaseq.88.unpruned.vcf.gz.tbi",
        # bed = "data/expression/ensemblGeneV2_Norm-counts-log2+1.bed.gz",
        # bedi = "data/expression/ensemblGeneV2_Norm-counts-log2+1.bed.gz.tbi"
        bed = "data/expression/expr_short.bed.gz",
        bedi = "data/expression/expr_short.bed.gz.tbi"
    output:
        "data/fastqtl/basic2.genes.txt.gz"
    conda:
        "envs/fastqtl.yaml"
    shell:
        """
        tools/fastqtl/python/run_FastQTL_threaded.py \
        {input.vcf} {input.bed} data/fastqtl/basic2 \
        --permute 1000 10000
        """
