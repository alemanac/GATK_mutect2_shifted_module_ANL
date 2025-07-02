#generates indexes for ref

rule gen_ref_bwa_indexes:
    input:
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}"
    output:
        ref_amb = f"{main_dir}/{ref_base}/ref_processing/{ref_base}{ref_ext}.amb",
        ref_bwt = f"{main_dir}/{ref_base}/ref_processing/{ref_base}{ref_ext}.bwt",
        ref_ann = f"{main_dir}/{ref_base}/ref_processing/{ref_base}{ref_ext}.ann",
        ref_pac = f"{main_dir}/{ref_base}/ref_processing/{ref_base}{ref_ext}.pac",
        ref_sa = f"{main_dir}/{ref_base}/ref_processing/{ref_base}{ref_ext}.sa"
    params:
        ref_fna = f"{reference_fasta}",
        ref_base = f"{ref_base}",
        rule_dir = f"{main_dir}/{ref_base}/ref_processing"
    log: 
        stderr = f"{main_dir}/{ref_base}/ref_processing/logs/{ref_base}_gen_ref_bwa_indexes_stderr",
        stdout = f"{main_dir}/{ref_base}/ref_processing/logs/{ref_base}_gen_ref_bwa_indexes_stdout"
    shell:
        """
        
        CWD=$(pwd); \
        RULEDIR="$(dirname -- "$(realpath -- "{output.ref_amb}")")"; \
        cd "$RULEDIR"; \
        
        bwa index {params.ref_fna} 2> "$CWD"/{log.stderr} > "$CWD"/{log.stdout}
        """

rule gen_ref_faidx:
    input:
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}"
    output:
        fai = f"{main_dir}/{ref_base}/ref_processing/{ref_base}{ref_ext}.fai"
    log:
        stderr = f"{main_dir}/{ref_base}/ref_processing/logs/{ref_base}_gen_ref_faidx_stderr",
        stdout = f"{main_dir}/{ref_base}/ref_processing/logs/{ref_base}_gen_ref_faidx_stdout"
    shell:
        """
        samtools faidx {input.ref} \
        -o {output.fai} 2> {log.stderr} > {log.stdout}
        """

rule gen_ref_dict:
    input:
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}"
    output:
        ref_dict = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.dict"
    log: 
        stderr = f"{main_dir}/{ref_base}/ref_processing/logs/{ref_base}_gen_ref_dict_stderr",
        stdout = f"{main_dir}/{ref_base}/ref_processing/logs/{ref_base}_gen_ref_dict_stdout"
    shell:
        """
        java {config[java_args]} \
        -jar {config[gatk_jar]} \
        CreateSequenceDictionary \
        -R {config[ref_fna]} \
        -O {output.ref_dict} 2> {log.stderr} > {log.stdout}
        """    

rule move_ref:
    output:
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}"
    shell:
        """
        cp {config[ref_fna]} {output.ref}
        """

rule move_ref_python_pipeline:
    output:
        ref = f"{main_dir}/run_gatk_python_pipeline/{reference_fasta}"