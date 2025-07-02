use rule gen_ref_bwa_indexes as gen_ref_bwa_indexes_shifted with: 
    input:
        ref = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted{ref_ext}"
    params:
        ref_fna = f"{ref_base}.shifted{ref_ext}",
        ref_base = f"{ref_base}.shifted",
        rule_dir = f"{main_dir}/{ref_base}/ref_processing_shifted"
    output:
        ref_amb = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted{ref_ext}.amb",
        ref_ann = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted{ref_ext}.ann",
        ref_bwt = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted{ref_ext}.bwt",
        ref_pac = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted{ref_ext}.pac",
        ref_sa = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted{ref_ext}.sa"
    log:
        stderr = f"{main_dir}/{ref_base}/ref_processing_shifted/logs/{ref_base}_gen_ref_bwa_indexes_shifted_stderr",
        stdout = f"{main_dir}/{ref_base}/ref_processing_shifted/logs/{ref_base}_gen_ref_bwa_indexes_shifted_stdout"
    container:
        "docker://us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"

rule shift_fasta:
    input:
        ref = f"{main_dir}/{ref_base}/ref_processing_shifted/{reference_fasta}",
        ref_index = f"{main_dir}/{ref_base}/ref_processing/{ref_base}{ref_ext}.fai",
        ref_dict = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.dict"
    output:
        shifted_ref_fna = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted{ref_ext}",
        shifted_ref_fna_index = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted{ref_ext}.fai",
        shiftback_chain = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shiftback.chain",
        shifted_intervals = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted.intervals",
        unshifted_intervals = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.intervals",
        shifted_dict = f"{main_dir}/{ref_base}/ref_processing_shifted/{ref_base}.shifted.dict"
    params:
        rule_dir = f"{main_dir}/{ref_base}/ref_processing_shifted",
        ref_base = ref_base
    log: 
        stderr = f"{main_dir}/{ref_base}/ref_processing_shifted/logs/{ref_base}_shift_fasta_stderr",
        stdout = f"{main_dir}/{ref_base}/ref_processing_shifted/logs/{ref_base}_shift_fasta_stdout"
    container:
        "docker://broadinstitute/gatk"
    shell:
        """
        RULEDIR="$(dirname -- "$(realpath -- "{output.shifted_ref_fna}")")"; \

        cp {input.ref_index} "$RULEDIR"; \
        cp {input.ref_dict} "$RULEDIR"; \

        gatk ShiftFasta \
        -R {input.ref} \
        -O {output.shifted_ref_fna} \
        --interval-file-name "$RULEDIR"/{ref_base} \
        --shift-back-output {output.shiftback_chain} 2> {log.stderr} > {log.stdout}; \

        rm "$RULEDIR"/{params.ref_base}{ref_ext}.fai; \
        rm "$RULEDIR"/{params.ref_base}.dict
        """

use rule move_ref as move_ref_shifted with:
    output:
        ref = f"{main_dir}/{ref_base}/ref_processing_shifted/{reference_fasta}"