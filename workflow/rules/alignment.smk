reference_fasta = Path(config['ref_fna']).name
ref_base = Path(reference_fasta).with_suffix('')

"""
I'm unsure what this is sorting, but I think it's the indexes for the reads
"""
rule gen_sorted_reads:
    input:
        reads = "{main_dir}/{SRR}/gen_duplicates_marked_reads/reads.bam",
    output:
        sorted_reads = "{main_dir}/{SRR}/gen_sorted_reads/reads.bam"
    log: 
        stderr = "{main_dir}/{SRR}/gen_sorted_reads/stderr",
        stdout = "{main_dir}/{SRR}/gen_sorted_reads/stdout"
    shell:
        """
        java {config[java_args]} \
        -jar {config[gatk_jar]} \
        SortSam \
        INPUT={input.reads} \
        OUTPUT={output.sorted_reads} \
        SORT_ORDER="coordinate" \
        CREATE_INDEX=true \
        MAX_RECORDS_IN_RAM=300000 2> {log.stderr} > {log.stdout}
        """

#TODO: consider playing with the optical_duplicate_pixel_distance
rule gen_duplicates_marked_reads:
    input:
        reads = "{main_dir}/{SRR}/gen_merged_aligned_reads/reads.bam"
    output:
        duplicates_marked_reads = "{main_dir}/{SRR}/gen_duplicates_marked_reads/reads.bam",
        metric_file = "{main_dir}/{SRR}/gen_duplicates_marked_reads/metrics.txt"
    log: 
        stderr = "{main_dir}/{SRR}/gen_duplicates_marked_reads/stderr",
        stdout = "{main_dir}/{SRR}/gen_duplicates_marked_reads/stdout"
    shell:
        """
        java {config[java_args]} \
        -jar {config[gatk_jar]} \
        MarkDuplicates \
        INPUT={input.reads} \
        OUTPUT={output.duplicates_marked_reads} \
        METRICS_FILE={output.metric_file} \
        VALIDATION_STRINGENCY=SILENT \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        ASSUME_SORT_ORDER="queryname" \
        CLEAR_DT="false" \
        ADD_PG_TAG_TO_READS=false 2> {log.stderr} > {log.stdout}
        """

rule gen_merged_aligned_reads: 
    input:
        aligned_reads = "{main_dir}/{SRR}/gen_aligned_reads/reads.bam",
        unmapped_reads = "{main_dir}/{SRR}/gen_unmapped_bam/reads.unmapped.bam",
        ref_fa = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}",
        ref_dict = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.dict"
    output:
        m_b = "{main_dir}/{SRR}/gen_merged_aligned_reads/reads.bam"
    params:
        bwa_command = "bwa mem -K 100000000 -v 3 -t 2 -Y BW25113_reference_fasta {SRR}_1.fastq {SRR}_2.fastq -o {main_dir}/{SRR}.bam"
    log: 
        stdout = "{main_dir}/{SRR}/gen_merged_aligned_reads/stdout",
        stderr = "{main_dir}/{SRR}/gen_merged_aligned_reads/stderr"
    shell:
        """
        java {config[java_args]} \
        -jar {config[gatk_jar]} \
        MergeBamAlignment \
        VALIDATION_STRINGENCY=SILENT \
        EXPECTED_ORIENTATIONS=FR \
        ATTRIBUTES_TO_RETAIN=X0 \
        ATTRIBUTES_TO_REMOVE=NM \
        ATTRIBUTES_TO_REMOVE=MD \
        ALIGNED_BAM={input.aligned_reads} \
        UNMAPPED_BAM={input.unmapped_reads} \
        OUTPUT={output.m_b} \
        REFERENCE_SEQUENCE={input.ref_fa} \
        PAIRED_RUN=true \
        SORT_ORDER="unsorted" \
        IS_BISULFITE_SEQUENCE=false \
        ALIGNED_READS_ONLY=false \
        CLIP_ADAPTERS=false \
        MAX_RECORDS_IN_RAM=2000000 \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
        ALIGNER_PROPER_PAIR_FLAGS=true \
        UNMAP_CONTAMINANT_READS=true \
        ADD_PG_TAG_TO_READS=false 2> {log.stderr} > {log.stdout}
        """ 

rule gen_unmapped_bam:
    input:
        fq_1 = "{main_dir}/{SRR}/fetch_reads/{SRR}_1.fastq",
        fq_2 = "{main_dir}/{SRR}/fetch_reads/{SRR}_2.fastq" 
    output:
        unmapped_reads = "{main_dir}/{SRR}/gen_unmapped_bam/reads.unmapped.bam"
    log: 
        stderr = "{main_dir}/{SRR}/gen_unmapped_bam/stderr",
        stdout = "{main_dir}/{SRR}/gen_unmapped_bam/stdout"
    shell: 
        """
        java {config[java_args]} \
        -jar {config[gatk_jar]} \
        FastqToSam \
        --FASTQ {input.fq_1} \
        --FASTQ2 {input.fq_2} \
        --OUTPUT {output.unmapped_reads} \
        --SAMPLE_NAME {wildcards.SRR} \
        --LIBRARY_NAME {config[read_group_library]} \
        --PLATFORM_UNIT {config[read_group_platform_unit]} \
        --RUN_DATE 0000-00-00 \
        --PLATFORM {config[read_group_platform]} \
        --SEQUENCING_CENTER {config[read_group_sequencing_center]} 2> {log.stderr} > {log.stdout}
        """

rule gen_aligned_reads:
    input:
        fq_1 = "{main_dir}/{SRR}/fetch_reads/{SRR}_1.fastq",
        fq_2 = "{main_dir}/{SRR}/fetch_reads/{SRR}_2.fastq",
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}",
        ref_amb = f"{main_dir}/{ref_base}/ref_processing/{ref_base}{ref_ext}.amb",
        ref_bwt = f"{main_dir}/{ref_base}/ref_processing/{ref_base}{ref_ext}.bwt",
        ref_ann = f"{main_dir}/{ref_base}/ref_processing/{ref_base}{ref_ext}.ann",
        ref_pac = f"{main_dir}/{ref_base}/ref_processing/{ref_base}{ref_ext}.pac",
        ref_sa = f"{main_dir}/{ref_base}/ref_processing/{ref_base}{ref_ext}.sa",
    output:
        aligned_reads = "{main_dir}/{SRR}/gen_aligned_reads/reads.bam"
    params:
        ref_processing_dir = f"{main_dir}/{ref_base}/ref_processing"
    log: 
        stderr = "{main_dir}/{SRR}/gen_aligned_reads/stderr"
    shell:
        """
        WORKDIR=$(pwd); \
        RULEDIR="$(dirname -- "$(realpath -- "{output.aligned_reads}")")"; \
        
        cd "$RULEDIR"; \

        bwa mem \
        -K 100000000 \
        -v 3 \
        -t 2 \
        -Y \
        "$WORKDIR"/{input.ref} "$WORKDIR"/{input.fq_1} "$WORKDIR"/{input.fq_2} \
        2> "$WORKDIR"/{log.stderr} \
        > reads.bam
        """

rule fetch_reads:
    output:
        fq_1 = "{main_dir}/{SRR}/fetch_reads/{SRR}_1.fastq",
        fq_2 = "{main_dir}/{SRR}/fetch_reads/{SRR}_2.fastq"
    log: 
        stderr = "{main_dir}/{SRR}/fetch_reads/{SRR}_fetch_reads_stderr",
        stdout = "{main_dir}/{SRR}/fetch_reads/{SRR}_fetch_reads_stdout"
    shell:
        """
        RULEDIR="$(dirname -- "$(realpath -- "{output.fq_1}")")"; \

        fasterq-dump {wildcards.SRR} \
        -O "$RULEDIR"
        """