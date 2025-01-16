def fastq_paths(wildcards):
    ## I've now renamed the fastq files so no need to swap.
    # # These two samples have swapped IDs:
    # if wildcards.sample_id == "000789FFF0_LHB":
    #     sample = "000789FFF9_LHB"
    # elif wildcards.sample_id == "000789FFF9_LHB":
    #     sample = "000789FFF0_LHB"
    # else:
    #     sample = wildcards.sample_id
    files = fastqs.loc[fastqs["library"] == wildcards.sample_id]
    files = files.to_dict(orient="records")
    return ["data/fastq/batch{}/{}".format(x["batch"], x["filename"]) for x in files]


rule star_index:
    input:
        fasta = "data/Rnor_6.0/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa",
        gtf = "data/Rnor_6.0_anno/Rattus_norvegicus.Rnor_6.0.99.gtf"
    output:
        # Among others:
        "/gpfs/group/pejman/dmunro/star_index/SAindex"
    params:
        outdir = "/gpfs/group/pejman/dmunro/star_index"
    # threads: 8
    resources:
        cpus = 8
    shell:
        """
        STAR --runMode genomeGenerate \
        --genomeDir {params.outdir} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 99 \
        --runThreadN 8
        """


def vcf_path(wildcards):
    rat_id = wildcards.sample_id.split("_")[0]
    return "data/genotype/individual/{}.vcf.gz".format(rat_id)


def read_groups(wildcards, input):
    rgs = expand("ID:{fq} SM:{sam}", fq=input.fastq, sam=wildcards.sample_id)
    return " , ".join(rgs)


rule star_align:
    input:
        fastq = fastq_paths,
        vcf = vcf_path,
        index = "data/star_index/SAindex"
    output:
        "data/star_out/{tissue}/{sample_id}.Aligned.sortedByCoord.out.bam"
    params:
        fastq_list = lambda wildcards, input: ",".join(input.fastq),
        index_dir = "/gpfs/group/pejman/dmunro/star_index",
        prefix = "data/star_out/{tissue}/{sample_id}.",
        read_groups = read_groups,
    resources:
        mem_mb = 60000,
        cpus = 16,
        walltime = 6
    shell:
        """
        STAR --runMode alignReads \
        --runThreadN 16 \
        --genomeDir {params.index_dir} \
        --readFilesIn {params.fastq_list} \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --varVCFfile <(zcat {input.vcf}) \
        --waspOutputMode SAMtag \
        --outSAMstrandField intronMotif \
        --outSAMattrRGline {params.read_groups} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outFileNamePrefix {params.prefix}
        """


rule filter_wasp:
    input:
        "data/star_out/{tissue}/{sample_id}.Aligned.sortedByCoord.out.bam"
    output:
        "data/star_out/{tissue}/{sample_id}.Aligned.sortedByCoord.out.wasp.bam"
    resources:
        cpus = 16,
        walltime = 6
    shell:
        """
        samtools view -H {input} > {output}.sam
        samtools view {input} | grep -v "vW:i:[2-7]" >> {output}.sam
        samtools view -1 {output}.sam --threads 16 > {output}
        rm {output}.sam
        """


rule mark_duplicates:
    input:
        # "data/star_out/{tissue}/{sample_id}.Aligned.sortedByCoord.out.bam"
        "data/star_out/{tissue}/{sample_id}.Aligned.sortedByCoord.out.wasp.bam"
    output:
        bam = "data/markdup_out/{tissue}/{sample_id}.bam",
        metrics = "data/markdup_out/{tissue}/{sample_id}.marked_dup_metrics.txt"
    shell:
        # MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500
        """
        picard MarkDuplicates \
        I={input} \
        O={output.bam} \
        M={output.metrics} \
        ASSUME_SORT_ORDER=coordinate \
        PROGRAM_RECORD_ID=null \
        TMP_DIR=$TMPDIR \
        MAX_RECORDS_IN_RAM=2000000
        """
