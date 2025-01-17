rule retrieve_fastq:
    output:
        "data/fastq/{batch}/{file}"
    shell:
        "rsync -av pejlab:~/data/fastq/{wildcards.batch}/{wildcards.file} data/fastq/{wildcards.batch}/"

rule retrieve_vcf:
    output:
        "data/vcfs_for_star/{vcf}"
    shell:
        "rsync -av pejlab:~/{output} data/vcfs_for_star/"

rule retrieve_sj:
    output:
        "data/star_out_sj/{sj}"
    shell:
        "rsync -av pejlab:~/{output} data/star_out_sj/"

rule bam_to_fastq:
    input:
        lambda wildcards: input_bams[wildcards.sample_id]
    output:
        "data/fastq_for_star/{sample_id}.fastq.gz"
    shell:
        """
        mkfifo read_pipe_{wildcards.sample_id}
        gzip -1 -c < read_pipe_{wildcards.sample_id} > {output} &
        java -jar -Xmx8g tools/picard.jar SamToFastq \
        I={input} FASTQ=read_pipe_{wildcards.sample_id} \
        INCLUDE_NON_PF_READS=true VALIDATION_STRINGENCY=SILENT
        rm read_pipe_{wildcards.sample_id}
        """

rule star_load_genome:
    output:
        "star_genome_loaded"
    shell:
        """
        STAR --genomeDir data/star_index --genomeLoad LoadAndExit
        touch {output}
        """


def star_unload_genome():
    if os.path.exists("star_genome_loaded"):
        shell("STAR --genomeDir data/star_index --genomeLoad Remove")
        os.remove("star_genome_loaded")


onsuccess:
    star_unload_genome()
onerror:
    star_unload_genome()

rule star_align_1st_pass:
    # Get novel junctions to use for 2nd pass.
    input:
        fastq = fastq_paths,
        index = "data/star_index/SAindex",
        genome = "star_genome_loaded"
    params:
        fastq_list = lambda wildcards, input: ",".join(input.fastq)
        # unsorted = "data/star_out_sj/{sample_id}.Aligned.out.bam",
        # trans = "data/star_out_sj/{sample_id}.Aligned.toTranscriptome.out.bam",
        # chim = "data/star_out_sj/{sample_id}.Chimeric.out.junction"
    output:
        # 2nd pass deletes this file, so output 1st pass in different dir:
        "data/star_out_sj/{sample_id}.SJ.out.tab"
    # threads: 1000       # One STAR job should run alone due to memory.
    shell:
        """
        STAR --runMode alignReads \
        --runThreadN 8 \
        --genomeDir data/star_index \
        --genomeLoad LoadAndKeep \
        --readFilesIn {params.fastq_list} \
        --readFilesCommand zcat \
        --outSAMtype None \
        --outFileNamePrefix data/star_out_sj/{wildcards.sample_id}.
        """

rule star_align_2nd_pass:
    input:
        fastq = fastq_paths,
        vcf = vcf_path,
        # sj = expand("data/star_out_sj/{sample_id}.SJ.out.tab", sample_id=sample_ids),
        sj = "data/star_out_sj/{sample_id}.SJ.out.tab",
        index = "data/star_index/SAindex"
    params:
        fastq_list = lambda wildcards, input: ",".join(input.fastq),
        unsorted = "data/star_out/{sample_id}.Aligned.out.bam",
        sjnew = "data/star_out/{sample_id}.SJ.out.tab"
    output:
        coord = "data/star_out/{sample_id}.Aligned.sortedByCoord.out.bam",
        coordi = "data/star_out/{sample_id}.Aligned.sortedByCoord.out.bam.bai",
        trans = "data/star_out/{sample_id}.Aligned.toTranscriptome.out.bam",
        chim = "data/star_out/{sample_id}.Chimeric.out.junction"
    threads: 1000       # One STAR job should run alone due to memory.
    shell:
        """
        STAR --runMode alignReads \
        --runThreadN 8 \
        --genomeDir data/star_index \
        --readFilesIn {params.fastq_list} \
        --readFilesCommand zcat \
        --varVCFfile <(zcat {input.vcf}) \
        --waspOutputMode SAMtag \
        --sjdbFileChrStartEnd {input.sj} \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --chimOutType Junctions WithinBAM SoftClip \
        --chimMainSegmentMultNmax 1 \
        --quantMode TranscriptomeSAM \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --outFileNamePrefix data/star_out/{wildcards.sample_id}.
        tools/samtools sort --threads 8 -o {output.coord} {params.unsorted}
        tools/samtools index {output.coord}
        rm {params.unsorted} {params.sjnew}
        """

rule star_align_2nd_pass:
    input:
        fastq = fastq_paths,
        vcf = vcf_path,
        # sj = expand("data/star_out_sj/{sample_id}.SJ.out.tab", sample_id=sample_ids),
        sj = "data/star_out_sj/{sample_id}.SJ.out.tab",
        index = "data/star_index/SAindex"
    output:
        coord = "data/star_out/{sample_id}.Aligned.sortedByCoord.out.bam",
        trans = "data/star_out/{sample_id}.Aligned.toTranscriptome.out.bam",
        chim = "data/star_out/{sample_id}.Chimeric.out.junction"
    shell:
        """
        rsync -av garibaldi:~/data/star_out/{wildcards.sample_id}.* data/star_out/
        """

# Need to get bamsync binary.  Source code included in gtex pipeline.
rule bamsync:
    # sync BAMs (optional; copy QC flags and read group IDs)
    input:
        lambda wildcards: input_bams[wildcards.sample_id],
        "/secure/dmunro/star_out/{sample_id}.Aligned.sortedByCoord.out.bam"
    output:
        "/secure/dmunro/star_out/{sample_id}.Aligned.sortedByCoord.out.patched.bam"
    shell:
docker run --rm -v /secure/dmunro:/secure/dmunro -t broadinstitute/gtex_rnaseq:V8 \
    /bin/bash -c "/src/run_bamsync.sh \
        $input_bam \
        /secure/dmunro/star_out/${sample_id}.Aligned.sortedByCoord.out.bam \
        /secure/dmunro/star_out/${sample_id}"

rule rnaseqc:
    input:
        bam = "data/markdup_out/{sample_id}.Aligned.sortedByCoord.out.md.bam",
        gtf = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.gtf"
    params:
        interm1 = "data/rnaseqc_out/{sample_id}.exon_reads.gct",
        interm2 = "data/rnaseqc_out/{sample_id}.gene_reads.gct",
        interm3 = "data/rnaseqc_out/{sample_id}.gene_tpm.gct"
    output:
        "data/rnaseqc_out/{sample_id}.exon_reads.gct.gz",
        "data/rnaseqc_out/{sample_id}.gene_reads.gct.gz",
        "data/rnaseqc_out/{sample_id}.gene_tpm.gct.gz",
        "data/rnaseqc_out/{sample_id}.metrics.tsv"
    shell:
        """
        tools/rnaseqc {input.gtf} {input.bam} data/rnaseqc_out \
        -s {wildcards.sample_id} --unpaired
        # Seems identical to gene_reads file:
        rm data/rnaseqc_out/{wildcards.sample_id}.gene_fragments.gct
        gzip {params.interm1} {params.interm2} {params.interm3}
        """

rule combine_gcts:
    input:
        expand("data/rnaseqc_out/{sample_id}.{{quant_type}}.gct.gz",
               sample_id=sample_ids)
    output:
        "data/{quant_type,\w+}.gct.gz"
    shell:
        "python3 src/combine_GCTs.py {input} {wildcards.quant_type} -o data/"


# rule add_read_groups:
#     # Needed for GATK
#     input:
#         "data/star_out/{sample_id}.Aligned.sortedByCoord.out.bam"
#     output:
#         "data/with_rg/{sample_id}.Aligned.sortedByCoord.out.rg.bam"
#     group:
#         "asereadcounter"
#     shell:
#         """
#         picard AddOrReplaceReadGroups \
#         I={input} \
#         O={output} \
#         RGID=rg1 \
#         RGLB=lib1 \
#         RGPL=unknown \
#         RGPU=unknown \
#         RGSM={wildcards.sample_id}
#         """

rule ASEReadCounter:
    input:
        ref = "data/Rnor_6.0/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa",
        refi = "data/Rnor_6.0/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.fai",
        refd = "data/Rnor_6.0/Rattus_norvegicus.Rnor_6.0.dna.toplevel.dict",
        bam = "data/with_rg/{sample_id}.Aligned.sortedByCoord.out.rg.bam",
        vcf = vcf_path,
        vcfi = lambda wildcards: vcf_path(wildcards) + ".tbi"
    output:
        "data/asereadcounter/{sample_id}.readcounts.txt"
    group:
        "asereadcounter"
    conda:
        "envs/gatk.yaml"
    shell:
        """
        gatk ASEReadCounter \
        -R {input.ref} \
        -I {input.bam} \
        -V {input.vcf} \
        -O {output} \
        --min-depth-of-non-filtered-base 1 \
        --min-mapping-quality 250 \
        --min-base-quality 15 \
        -DF NotDuplicateReadFilter
        """

# def aFC_inputs(wildcards):
#     vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz"
#     bed = {
#         "basic": "data/expression/{region}.rsem_expected_count.bed.gz",
#         "basic2": "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz",
#         "basic3": "data/expression/{region}.rsem_TPM.bed.gz",
#         "basic4": "data/expression/{region}.rsem_expected_count.bed.gz",
#         "main": "data/expression/{region}.rsem_expected_count.bed.gz",
#         "main2": "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz",
#         "main3": "data/expression/{region}.rsem_TPM.bed.gz",
#         "main4": "data/expression/{region}.rsem_expected_count.bed.gz",
#         "qtl2": "data/expression/{region}.rsem_expected_count.bed.gz",
#     }
#     bed = bed[wildcards.method].format(region=wildcards.region)
#     covar = {
#         "main": "data/tensorqtl/{region}/{region}.main.combined_covariates.txt",
#         "main2": "data/tensorqtl/{region}/{region}.main2.combined_covariates.txt",
#         "main3": "data/tensorqtl/{region}/{region}.main3.combined_covariates.txt",
#         "main4": "data/tensorqtl/{region}/{region}.main4.combined_covariates.txt",
#         "qtl2": "data/tensorqtl/other_covariates.tsv",
#     }
#     if wildcards.method == "qtl2":
#         qtl = "data/qtl2/{region}.gene_var_pval.tsv.gz"
#     else:
#         qtl = "data/tensorqtl/{region}/{region}.{method}.cis_qtl.txt.gz"
#     qtl = expand(qtl, region=wildcards.region, method=wildcards.method)
#     inputs = dict(
#         vcf = vcf,
#         vcfi = vcf + ".tbi",
#         bed = bed,
#         bedi = bed + ".tbi",
#         qtl = qtl,
#     )
#     if wildcards.method not in {"basic", "basic2", "basic3", "basic4"}:
#         inputs["covar"] = covar[wildcards.method].format(region=wildcards.region)
#     return inputs


# rule aFC_from_eQTL_model_main:
#     input:
#         vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
#         vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
#         # bed = "data/tensorqtl/{region}/{region}.main.expression.bed.gz",
#         # bedi = "data/tensorqtl/{region}/{region}.main.expression.bed.gz.tbi",
#         bed = "data/expression/{region}.rsem_expected_count.bed.gz",
#         bedi = "data/expression/{region}.rsem_expected_count.bed.gz.tbi",
#         qtl = "data/tensorqtl/{region}/{region}.main.cis_qtl.txt.gz",
#         covar = "data/tensorqtl/{region}/{region}.main.combined_covariates.txt"
#     output:
#         "data/afc/{region}.main.aFC.txt"
#     conda:
#         "envs/afc.yaml"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python3 tools/aFC/aFC.py \
#         --vcf {input.vcf} \
#         --pheno {input.bed} \
#         --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl}) \
#         --cov <(sed 's/-N//g' {input.covar}) \
#         --log_xform 1 \
#         --output {output}
#         """

# rule aFC_from_eQTL_model_main2:
#     input:
#         vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
#         vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
#         bed = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz",
#         bedi = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz.tbi",
#         qtl = "data/tensorqtl/{region}/{region}.main2.cis_qtl.txt.gz",
#         covar = "data/tensorqtl/{region}/{region}.main2.combined_covariates.txt"
#     output:
#         "data/afc/{region}.main2.aFC.txt"
#     conda:
#         "envs/afc.yaml"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python3 tools/aFC/aFC.py \
#         --vcf {input.vcf} \
#         --pheno {input.bed} \
#         --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl}) \
#         --cov <(sed 's/-N//g' {input.covar}) \
#         --log_xform 1 \
#         --output {output}
#         """

# rule aFC_from_eQTL_model_main3:
#     input:
#         vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
#         vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
#         # bed = "data/tensorqtl/{region}/{region}.main.expression.bed.gz",
#         # bedi = "data/tensorqtl/{region}/{region}.main.expression.bed.gz.tbi",
#         bed = "data/expression/{region}.rsem_TPM.bed.gz",
#         bedi = "data/expression/{region}.rsem_TPM.bed.gz.tbi",
#         qtl = "data/tensorqtl/{region}/{region}.main3.cis_qtl.txt.gz",
#         covar = "data/tensorqtl/{region}/{region}.main3.combined_covariates.txt"
#     output:
#         "data/afc/{region}.main3.aFC.txt"
#     conda:
#         "envs/afc.yaml"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python3 tools/aFC/aFC.py \
#         --vcf {input.vcf} \
#         --pheno {input.bed} \
#         --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl}) \
#         --cov <(sed 's/-N//g' {input.covar}) \
#         --log_xform 1 \
#         --output {output}
#         """

# rule aFC_from_eQTL_model_basic:
#     input:
#         vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
#         vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
#         # bed = "data/tensorqtl/{region}/{region}.main.expression.bed.gz",
#         # bedi = "data/tensorqtl/{region}/{region}.main.expression.bed.gz.tbi",
#         bed = "data/expression/{region}.rsem_expected_count.bed.gz",
#         bedi = "data/expression/{region}.rsem_expected_count.bed.gz.tbi",
#         qtl = "data/tensorqtl/{region}/{region}.basic.cis_qtl.txt.gz"
#     output:
#         "data/afc/{region}.basic.aFC.txt"
#     conda:
#         "envs/afc.yaml"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python3 tools/aFC/aFC.py \
#         --vcf {input.vcf} \
#         --pheno {input.bed} \
#         --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl}) \
#         --log_xform 1 \
#         --output {output}
#         """

# rule aFC_from_eQTL_model_basic2:
#     input:
#         vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
#         vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
#         bed = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz",
#         bedi = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz.tbi",
#         qtl = "data/tensorqtl/{region}/{region}.basic2.cis_qtl.txt.gz"
#     output:
#         "data/afc/{region}.basic2.aFC.txt"
#     conda:
#         "envs/afc.yaml"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python3 tools/aFC/aFC.py \
#         --vcf {input.vcf} \
#         --pheno {input.bed} \
#         --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl}) \
#         --log_xform 1 \
#         --output {output}
#         """

# rule aFC_from_eQTL_model_basic3:
#     input:
#         vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
#         vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
#         # bed = "data/tensorqtl/{region}/{region}.main.expression.bed.gz",
#         # bedi = "data/tensorqtl/{region}/{region}.main.expression.bed.gz.tbi",
#         bed = "data/expression/{region}.rsem_TPM.bed.gz",
#         bedi = "data/expression/{region}.rsem_TPM.bed.gz.tbi",
#         qtl = "data/tensorqtl/{region}/{region}.basic3.cis_qtl.txt.gz"
#     output:
#         "data/afc/{region}.basic3.aFC.txt"
#     conda:
#         "envs/afc.yaml"
#     resources:
#         walltime = 12
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
#         vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
#         vcfi = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz.tbi",
#         bed = "data/tensorqtl/{region}/{region}.main.expression.bed.gz",
#         bedi = "data/tensorqtl/{region}/{region}.main.expression.bed.gz.tbi",
#         qtl = "data/qtl2/{region}.gene_var_pval.tsv.gz",
#         covar = "data/tensorqtl/other_covariates.tsv"
#     output:
#         "data/afc/{region}.qtl2.aFC.txt"
#     conda:
#         "envs/afc.yaml"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python3 tools/aFC/aFC.py \
#         --vcf {input.vcf} \
#         --pheno {input.bed} \
#         --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl} 1 2) \
#         --cov {input.covar} \
#         --log_xform 1 \
#         --output {output}
#         """

# rule phaser_expr_matrix:
#     input:
#         SFTP.remote("garibaldi.scripps.edu/gpfs/home/dmunro/data/phaser_pop_out/{region}.expr_matrix.gw_phased.bed.gz")
#     output:
#         "data/expression/{region}.expr_matrix.gw_phased.bed.gz"
#     shell:
#         "cp {input} {output}"

# rule bed_index:
#     input:
#         "{base}.bed.gz"
#     output:
#         "{base}.bed.gz.tbi"
#     shell:
#         "tabix {input}"

# rule phaser_sample_map:
#     output:
#         "data/phaser_pop_out/sample_map.txt"
#     run:
#         with open(output[0], "w") as out:
#             out.write("vcf_sample\tbed_sample\n")
#             for sample_id in sample_ids:
#                 rat_id = sample_id.split("_")[0]
#                 bam_base = "{}.Aligned.sortedByCoord.out.md".format(sample_id)
#                 out.write("{}\t{}\n".format(rat_id, bam_base))

# rule phaser_gene_var_pairs:
#     input:
#         "data/phaser_pop_out/{version}.aFC.txt"
#     output:
#         "data/phaser_pop_out/{version}.pairs.txt"
#     # shell:
#     #     """
#     #     echo "gene_id\tvar_contig" > {output}
#     #     tail -n+2 {input} | cut -f1,2 >> {output}
#     #     """
#     run:
#         inlines = open(input[0], "r").read().splitlines()
#         gene_var = [x.split("\t")[:2] for x in inlines[1:]]
#         with open(output[0], "w") as out:
#             out.write("gene_id\tvar_id\tvar_contig\tvar_pos\tvar_ref\tvar_alt\n")
#             for x in gene_var:
#                 chrom, pos = tuple(x[1].split(":"))
#                 chrom = chrom.replace("chr", "")
#                 out.write("{}\t{}\t{}\t{}\t\t\n".format(x[0], x[1], chrom, pos))

# rule phaser_cis_var:
#     input:
#         bed = "data/phaser_pop_out/expr_matrix.bed.gz",
#         vcf = "data/genotype/P50.rnaseq.88.unpruned.vcf.gz",
#         gene_var_pairs = "data/phaser_pop_out/{version}.pairs.txt",
#         sample_map = "data/phaser_pop_out/sample_map.txt"
#     output:
#         "data/phaser_pop_out/{version}.ASE_aFC.txt"
#     conda:
#         "envs/phaser.yaml"
#     threads: 16
#     shell:
#         """
#         python tools/phaser/phaser_pop/phaser_cis_var.py \
#         --bed {input.bed} \
#         --vcf {input.vcf} \
#         --pair {input.gene_var_pairs} \
#         --map {input.sample_map} \
#         --t 16 \
#         --o {output}
#         """

# python tools/phaser/phaser_gene_ae/phaser_gene_ae.py \
#     --haplotypic_counts data/phaser_out/00077E67B5_Acbc.haplotypic_counts.txt \
#     --features data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.bed \
#     --min_haplo_maf 0.05 \
#     --o data/phaser_out/00077E67B5_Acbc.gene_ae.txt

# python3 tools/phaser/phaser_pop/phaser_expr_matrix.py \
#     --gene_ae_dir data/phaser_gene_ae_out \
#     --features data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.bed \
#     --o data/phaser_pop_out/expr_matrix

# python tools/phaser/phaser_pop/phaser_cis_var.py \
#     --bed data/phaser_pop_out/expr_matrix.bed.gz \
#     --vcf data/genotype/P50.rnaseq.88.unpruned.vcf.gz \
#     --pair data/phaser_pop_out/main.pairs.txt \
#     --map data/phaser_pop_out/sample_map.txt \
#     --o data/phaser_pop_out/main.ASE_aFC.txt

# python tools/phaser/wrapper.py phase 00077E67B5 data/markdup_out/00077E67B5_Acbc.Aligned.sortedByCoord.out.md.bam data/genotype/P50.rnaseq.88.unpruned.vcf.gz data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.genes.bed data/phaser_out

# python tools/phaser/phaser/phaser.py \
#     --temp_dir $PBSTMPDIR \
#     --bam data/markdup_out/00077E67B5_Acbc.Aligned.sortedByCoord.out.md.bam \
#     --vcf data/genotype/P50.rnaseq.88.unpruned.vcf.gz \
#     --sample 00077E67B5 \
#     --baseq 10 \
#     --mapq 255 \
#     --isize 1e6 \
#     --paired_end 0 \
#     --o data/phaser_out/00077E67B5 \
#     --include_indels 0 \
#     --gw_phase_vcf 1

# Join chromosome VCFs if --chr

# bgzip VCF?

# tabix VCF?









# rule run_tensorqtl_perm_main:
#     input:
#         geno = multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam"),
#         bed = "data/expression/{region}.rsem_expected_count.bed.gz",
#         bedi = "data/expression/{region}.rsem_expected_count.bed.gz.tbi",
#         covar = "data/tensorqtl/{region}/{region}.main.combined_covariates.txt"
#     output:
#         "data/tensorqtl/{region}/{region}.main.cis_qtl.txt.gz"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python -m tensorqtl --mode cis \
#             data/genotype/P50.rnaseq.88.unpruned \
#             {input.bed} \
#             {wildcards.region}.main \
#             --covariates {input.covar} \
#             --output_dir data/tensorqtl/{wildcards.region}
#         """

# rule run_tensorqtl_perm_main2:
#     input:
#         geno = multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam"),
#         bed = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz",
#         bedi = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz.tbi",
#         covar = "data/tensorqtl/{region}/{region}.main2.combined_covariates.txt"
#     output:
#         "data/tensorqtl/{region}/{region}.main2.cis_qtl.txt.gz"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python -m tensorqtl --mode cis \
#             data/genotype/P50.rnaseq.88.unpruned \
#             {input.bed} \
#             {wildcards.region}.main2 \
#             --covariates {input.covar} \
#             --output_dir data/tensorqtl/{wildcards.region}
#         """

# rule run_tensorqtl_perm_main3:
#     input:
#         geno = multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam"),
#         bed = "data/expression/{region}.rsem_TPM.bed.gz",
#         bedi = "data/expression/{region}.rsem_TPM.bed.gz.tbi",
#         covar = "data/tensorqtl/{region}/{region}.main3.combined_covariates.txt"
#     output:
#         "data/tensorqtl/{region}/{region}.main3.cis_qtl.txt.gz"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python -m tensorqtl --mode cis \
#             data/genotype/P50.rnaseq.88.unpruned \
#             {input.bed} \
#             {wildcards.region}.main3 \
#             --covariates {input.covar} \
#             --output_dir data/tensorqtl/{wildcards.region}
#         """

# rule run_tensorqtl_perm_basic:
#     input:
#         geno = multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam"),
#         bed = "data/expression/{region}.rsem_expected_count.bed.gz",
#         bedi = "data/expression/{region}.rsem_expected_count.bed.gz.tbi",
#         covar = "data/tensorqtl/{region}/{region}.basic.covar_empty.txt"
#     output:
#         "data/tensorqtl/{region}/{region}.basic.cis_qtl.txt.gz"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python -m tensorqtl --mode cis \
#             data/genotype/P50.rnaseq.88.unpruned \
#             {input.bed} \
#             {wildcards.region}.basic \
#             --covariates {input.covar} \
#             --output_dir data/tensorqtl/{wildcards.region}
#         """

# rule run_tensorqtl_perm_basic2:
#     input:
#         geno = multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam"),
#         bed = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz",
#         bedi = "data/expression/{region}.ensemblGeneV2_Norm-counts-log2+1.bed.gz.tbi",
#         covar = "data/tensorqtl/{region}/{region}.basic2.covar_empty.txt"
#     output:
#         "data/tensorqtl/{region}/{region}.basic2.cis_qtl.txt.gz"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python -m tensorqtl --mode cis \
#             data/genotype/P50.rnaseq.88.unpruned \
#             {input.bed} \
#             {wildcards.region}.basic2 \
#             --covariates {input.covar} \
#             --output_dir data/tensorqtl/{wildcards.region}
#         """

# rule run_tensorqtl_perm_basic3:
#     input:
#         geno = multiext("data/genotype/P50.rnaseq.88.unpruned", ".bed", ".bim", ".fam"),
#         bed = "data/expression/{region}.rsem_TPM.bed.gz",
#         bedi = "data/expression/{region}.rsem_TPM.bed.gz.tbi",
#         covar = "data/tensorqtl/{region}/{region}.basic3.covar_empty.txt"
#     output:
#         "data/tensorqtl/{region}/{region}.basic3.cis_qtl.txt.gz"
#     resources:
#         walltime = 12
#     shell:
#         """
#         python -m tensorqtl --mode cis \
#             data/genotype/P50.rnaseq.88.unpruned \
#             {input.bed} \
#             {wildcards.region}.basic3 \
#             --covariates {input.covar} \
#             --output_dir data/tensorqtl/{wildcards.region}
#         """







# rule gzip_rsem_output:
#     input:
#         "data/rsem_out/{sample_id}.genes.results"
#     output:
#         "data/rsem_out/{sample_id}.genes.results.gz"
#     group:
#         "rsem"
#     shell:
#         "gzip {input}"


rule phaser_sample_map:
    output:
        "data/afc/{region}.sample_map.txt"
    run:
        df = pd.read_csv("data/bam/bam_files.txt", sep="\t")
        df = df[df["sample"].isin(sample_ids[wildcards.region])]
        df["vcf_sample"] = df["sample"].str.replace("_" + wildcards.region, "")
        df["bed_sample"] = df["file"].str.replace(".bam", "")
        df = df[["vcf_sample", "bed_sample"]]
        df.to_csv(output[0], sep="\t", index=False)
