reference_fasta = Path(config['ref_fna']).name
ref_base = Path(reference_fasta).with_suffix('')

use rule gen_sorted_reads as gen_sorted_reads_shifted with:
    input:
        reads = "{main_dir}/{SRR}/gen_duplicates_marked_reads_shifted/reads.bam",
    output:
        sorted_reads = "{main_dir}/{SRR}/gen_sorted_reads_shifted/reads.bam"
    log: 
        stderr = "{main_dir}/{SRR}/gen_sorted_reads_shifted/stderr",
        stdout = "{main_dir}/{SRR}/gen_sorted_reads_shifted/stdout"
    container: "docker://us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"

use rule gen_duplicates_marked_reads as gen_duplicates_marked_reads_shifted with:
    input:
        reads = "{main_dir}/{SRR}/gen_merged_aligned_reads_shifted/reads.bam"
    output:
        duplicates_marked_reads = "{main_dir}/{SRR}/gen_duplicates_marked_reads_shifted/reads.bam",
        metric_file = "{main_dir}/{SRR}/gen_duplicates_marked_reads_shifted/metrics.txt"
    log: 
        stderr = "{main_dir}/{SRR}/gen_duplicates_marked_reads_shifted/stderr",
        stdout = "{main_dir}/{SRR}/gen_duplicates_marked_reads_shifted/stdout"
    container:
        "docker://us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"

use rule gen_merged_aligned_reads as gen_merged_aligned_reads_shifted with: 
    input:
        aligned_reads = "{main_dir}/{SRR}/gen_aligned_reads_shifted/reads.bam",
        unmapped_reads = "{main_dir}/{SRR}/gen_unmapped_bam/reads.unmapped.bam",
        ref_fa = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted{ref_ext}",
        ref_dict = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted.dict"
    output:
        m_b = "{main_dir}/{SRR}/gen_merged_aligned_reads_shifted/reads.bam"
    params:
        bwa_command = "bwa mem -K 100000000 -v 3 -t 2 -Y BW25113_reference_fasta_shifted {SRR}_1.fastq {SRR}_2.fastq -o {main_dir}/{SRR}.bam"
    log: 
        stdout = "{main_dir}/{SRR}/gen_merged_aligned_reads_shifted/stdout",
        stderr = "{main_dir}/{SRR}/gen_merged_aligned_reads_shifted/stderr"
    container:
        "docker://us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"

use rule gen_aligned_reads as gen_aligned_reads_shifted with:
    input:
        fq_1 = "{main_dir}/{SRR}/fetch_reads/{SRR}_1.fastq",
        fq_2 = "{main_dir}/{SRR}/fetch_reads/{SRR}_2.fastq",
        ref = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted{ref_ext}",
        ref_amb = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted{ref_ext}.amb",
        ref_ann = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted{ref_ext}.ann",
        ref_bwt = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted{ref_ext}.bwt",
        ref_pac = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted{ref_ext}.pac",
        ref_sa = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted{ref_ext}.sa"
    output:
        aligned_reads = "{main_dir}/{SRR}/gen_aligned_reads_shifted/reads.bam"
    params:
        ref_processing_dir = f"{main_dir}/{ref_base}/ref_processing_shifted"
    log: 
        stderr = "{main_dir}/{SRR}/gen_aligned_reads_shifted/stderr",
        stdout = "{main_dir}/{SRR}/gen_aligned_reads_shifted/stdout"
    container: 
        "docker://us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
