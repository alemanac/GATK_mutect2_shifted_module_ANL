rule remove_duplicate_rows:
    input:
        variants = "{main_dir}/{SRR}/lifted_over_and_combined_vcfs/variants.final.vcf",
    output:
        variants = "{main_dir}/{SRR}/remove_duplicate_rows/variants.vcf"
    log:
        stderr = "{main_dir}/{SRR}/remove_duplicate_rows/stderr",
        stdout = "{main_dir}/{SRR}/remove_duplicate_rows/stdout"
    container:
        "docker://broadinstitute/gatk"
    shell:
        """
        bcftools norm -d none {input.variants} -o {output.variants}
        """


rule lifted_over_and_combined_vcfs:
    input:
        variants = "{main_dir}/{SRR}/gen_filtered_vcfs/variants_filtered.vcf",
        shifted_variants = "{main_dir}/{SRR}/gen_filtered_vcfs_shifted/variants_filtered.shifted.vcf",
        shiftback_chain = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shiftback.chain",
        ref_dict = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.dict",
        ref_fa = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}"
    output:
        shifted_back_variants = "{main_dir}/{SRR}/lifted_over_and_combined_vcfs/variants.shifted_back.vcf",
        rejected_vcf = "{main_dir}/{SRR}/lifted_over_and_combined_vcfs/variants.rejected.vcf",
        final_variants = "{main_dir}/{SRR}/lifted_over_and_combined_vcfs/variants.final.vcf",
        final_variants_index = "{main_dir}/{SRR}/lifted_over_and_combined_vcfs/variants.final.vcf.idx"
    log:
        stderr_lift_over = "{main_dir}/{SRR}/lifted_over_and_combined_vcfs/stderr",
        stdout_lift_over = "{main_dir}/{SRR}/lifted_over_and_combined_vcfs/stdout",
        stderr_merge = "{main_dir}/{SRR}/lifted_over_and_combined_vcfs/stderr",
        stdout_merge = "{main_dir}/{SRR}/lifted_over_and_combined_vcfs/stdout"
    container:
        "docker://broadinstitute/gatk"
    shell:
        """
        gatk LiftoverVcf \
        -I {input.shifted_variants} \
        -O {main_dir}/{wildcards.SRR}/lifted_over_and_combined_vcfs/variants.shifted_back.vcf \
        -R {input.ref_fa} \
        --CHAIN {input.shiftback_chain} \
        --REJECT {output.rejected_vcf} 2> {log.stderr_lift_over} > {log.stdout_lift_over}; \

        gatk MergeVcfs \
        -I {main_dir}/{wildcards.SRR}/lifted_over_and_combined_vcfs/variants.shifted_back.vcf \
        -I {input.variants} \
        -O {output.final_variants} 2> {log.stderr_merge} > {log.stdout_merge}
        """

"""
filters output .vcf file from GATK variation call
"""
rule gen_filtered_vcfs:
    input:
        variants = "{main_dir}/{SRR}/gen_mutect2_vcfs/variants.vcf",
        stats = "{main_dir}/{SRR}/gen_mutect2_vcfs/variants.vcf.stats",
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}"
    output:
        variants_filtered = "{main_dir}/{SRR}/gen_filtered_vcfs/variants_filtered.vcf"
    params:
        allele_fraction_thres = "0.0",
        microbial_mode = "--microbial-mode"
    log: 
        stderr = "{main_dir}/{SRR}/gen_filtered_vcfs/stderr",
        stdout = "{main_dir}/{SRR}/gen_filtered_vcfs/stdout"
    container:
        "docker://broadinstitute/gatk"
    shell:
        """
        gatk --java-options "-Xmx2500m" FilterMutectCalls \
        -V {input.variants} \
        -R {input.ref} \
        -O {output.variants_filtered} \
        --stats {input.stats} \
        {params.microbial_mode}\
        --min-allele-fraction {params.allele_fraction_thres} 2> {log.stderr} > {log.stdout}
        """	

"""
filters lifted & combined vcf
"""
use rule gen_filtered_vcfs as gen_filtered_vcfs_shifted with:
    input:
        variants ="{main_dir}/{SRR}/gen_mutect2_vcfs_shifted/variants.vcf",
        stats = "{main_dir}/{SRR}/gen_mutect2_vcfs_shifted/variants.vcf.stats",
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}"
    output:
        variants_filtered = "{main_dir}/{SRR}/gen_filtered_vcfs_shifted/variants_filtered.shifted.vcf"
    params:
        allele_fraction_thres = "0.0",
        microbial_mode = "--microbial-mode "
    log:
        stderr = "{main_dir}/{SRR}/gen_filtered_vcfs_shifted/stderr",
        stdout = "{main_dir}/{SRR}/gen_filtered_vcfs_shifted/stdout"
    container:
        "docker://broadinstitute/gatk"

#TODO: check whether to specify the sequence dictionary command
"""
performs a variation call via GATK
"""
rule gen_mutect2_vcfs: 
    input:
        reads = "{main_dir}/{SRR}/gen_mq_filtered_reads/reads.bam",
        ref = f"{main_dir}/{ref_base}/ref_processing/{reference_fasta}",
        fai = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.fna.fai" #although not explictly use in the shell, the command relies on it
    output:
        variants = "{main_dir}/{SRR}/gen_mutect2_vcfs/variants.vcf",
        stats = "{main_dir}/{SRR}/gen_mutect2_vcfs/variants.vcf.stats",
        local_assemblies = "{main_dir}/{SRR}/gen_mutect2_vcfs/local_assemblies.bam",
        assembly_regions = "{main_dir}/{SRR}/gen_mutect2_vcfs/assembly_regions.tsv",
        graphs = "{main_dir}/{SRR}/gen_mutect2_vcfs/graph"
    threads: 8
    log:
        stderr = "{main_dir}/{SRR}/gen_mutect2_vcfs/stderr",
        stdout = "{main_dir}/{SRR}/gen_mutect2_vcfs/stdout"
    container:
        "docker://broadinstitute/gatk"
    shell: 
        """
        gatk --java-options "-Xmx3000m" Mutect2 \
        -R {input.ref} \
        -I {input.reads} \
        -O {output.variants} \
        --bam-output {output.local_assemblies} \
        --annotation StrandBiasBySample \
        --num-matching-bases-in-dangling-end-to-recover 1 \
        --max-reads-per-alignment-start 75 \
        --native-pair-hmm-threads 8 \
        --assembly-region-out {output.assembly_regions} \
        --graph-output {output.graphs} 2> {log.stderr} > {log.stdout}
        """

"""
performs a variation call with GATK using basically all the shifted files
"""
use rule gen_mutect2_vcfs as gen_mutect2_vcfs_shifted with:
    input:
        reads = "{main_dir}/{SRR}/gen_mq_filtered_reads_shifted/reads.bam",
        ref = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna",
        fai = f"{main_dir}/{ref_base}/ref_processing/{ref_base}.shifted.fna.fai",
    output:
        variants = "{main_dir}/{SRR}/gen_mutect2_vcfs_shifted/variants.vcf",
        stats = "{main_dir}/{SRR}/gen_mutect2_vcfs_shifted/variants.vcf.stats",
        local_assemblies = "{main_dir}/{SRR}/gen_mutect2_vcfs_shifted/local_assemblies.shifted.bam",
        assembly_regions = "{main_dir}/{SRR}/gen_mutect2_vcfs_shifted/assembly_regions",
        graphs = "{main_dir}/{SRR}/gen_mutect2_vcfs_shifted/graphs"
    threads: 8
    log:
        stderr = "{main_dir}/{SRR}/gen_mutect2_vcfs_shifted/stderr",
        stdout = "{main_dir}/{SRR}/gen_mutect2_vcfs_shifted/stdout"
    container:
        "docker://broadinstitute/gatk"