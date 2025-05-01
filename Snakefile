#!/usr/bin/env python3
import peppy

#############
# FUNCTIONS #
#############

def get_reads(wildcards):
    input_keys = ['l1r1', 'l2r1', 'l3r1', 'l1r2', 'l2r2', 'l3r2']
    my_pep = pep.get_sample(wildcards.sample).to_dict()
    return {k: my_pep[k] for k in input_keys}

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_samples = pep.sample_table['sample_name']

## containers ##
fastqc_container = 'docker://staphb/fastqc:0.12.1'
multiqc_container = 'docker://staphb/multiqc:1.22.3'
trimgalore_container = 'docker://quay.io/biocontainers/trim-galore:0.6.9--hdfd78af_0'
salmon_container = 'docker://combinelab/salmon:1.10.3'
bioconductor_container = 'library://sinwood/bioconductor/bioconductor_3.20:0.0.1'

#########
# RULES #
#########

rule target:
    input:
        ## trimming and QC report
        'output/01_read_prep/multiqc/multiqc_report.html',
        ## salmon report
        'output/02_salmon/multiqc/multiqc_report.html',
        ## deseq2
        expand('output/03_deseq/PCA/PCA_timepoint_treatment_{seq_batch}.pdf', seq_batch=["both_batches", "first_batch"]),
        expand('output/03_deseq/time_treatment_pairwise_both_batches/{dds_file}/sig_annots_12hr.csv', dds_file=["both_batches", "both_batches_filtered"]),
        expand('output/03_deseq/time_treatment_pairwise_first_batches/{dds_file}/sig_annots_12hr.csv', dds_file=["first_batch", "first_batch_filtered"])

##########################
## 03 - deseq2 analysis ##
##########################

## far less
rule deseq2_time_treatment_pairwise_first_batch:
    input:
        dds_file = 'output/03_deseq/dds_{dds_file}.rds',
        gtf_file = 'data/misc_dmel-r6.63_files/dmel-all-r6.63.gtf'
    output:
        dds = 'output/03_deseq/time_treatment_pairwise_first_batches/{dds_file}/dds.rds',
        res_12hr = 'output/03_deseq/time_treatment_pairwise_first_batches/{dds_file}/sig_annots_12hr.csv',
        res_24hr = 'output/03_deseq/time_treatment_pairwise_first_batches/{dds_file}/sig_annots_24hr.csv',
        res_48hr = 'output/03_deseq/time_treatment_pairwise_first_batches/{dds_file}/sig_annots_48hr.csv',
        res_120hr = 'output/03_deseq/time_treatment_pairwise_first_batches/{dds_file}/sig_annots_120hr.csv'
    log:
        'output/logs/deseq2_time_treatment_pairwise_{dds_file}.log'
    singularity:
        bioconductor_container
    script:
        'src/deseq2_time_treatment_pairwise_first_batch.R'

## filtered or not pretty similar except for 120 hrs
rule deseq2_time_treatment_pairwise:
    input:
        dds_file = 'output/03_deseq/dds_{dds_file}.rds',
        gtf_file = 'data/misc_dmel-r6.63_files/dmel-all-r6.63.gtf'
    output:
        dds = 'output/03_deseq/time_treatment_pairwise_both_batches/{dds_file}/dds.rds',
        res_12hr = 'output/03_deseq/time_treatment_pairwise_both_batches/{dds_file}/sig_annots_12hr.csv',
        res_24hr = 'output/03_deseq/time_treatment_pairwise_both_batches/{dds_file}/sig_annots_24hr.csv',
        res_48hr = 'output/03_deseq/time_treatment_pairwise_both_batches/{dds_file}/sig_annots_48hr.csv',
        res_120hr = 'output/03_deseq/time_treatment_pairwise_both_batches/{dds_file}/sig_annots_120hr.csv'
    log:
        'output/logs/deseq2_time_treatment_pairwise_{dds_file}.log'
    singularity:
        bioconductor_container
    script:
        'src/deseq2_time_treatment_pairwise.R'

rule deseq2_treatment_pairwise:
    input:
        dds_file = 'output/03_deseq/dds_{dds_file}.rds',
        gtf_file = 'data/misc_dmel-r6.63_files/dmel-all-r6.63.gtf'
    output:
        dds = 'output/03_deseq/treatment_pairwise/both_batches/{dds_file}/{dds_file}_dds.rds',
        res_12hr = 'output/03_deseq/treatment_pairwise/both_batches/{dds_file}/{dds_file}_sig_annots_12hr.csv',
        res_24hr = 'output/03_deseq/treatment_pairwise/both_batches/{dds_file}/{dds_file}_sig_annots_24hr.csv',
        res_48hr = 'output/03_deseq/treatment_pairwise/both_batches/{dds_file}/{dds_file}_sig_annots_48hr.csv',
        res_120hr = 'output/03_deseq/treatment_pairwise/both_batches/{dds_file}/{dds_file}_sig_annots_120hr.csv'
    log:
        'output/logs/deseq2_time_treatment_pairwise_{dds_file}.log'
    singularity:
        bioconductor_container
    script:
        'src/deseq2_time_treatment_pairwise.R'

rule deseq2_pca:
    input:
        dds_file = 'output/03_deseq/dds_{seq_batch}.rds'
    output:
        PCA_time = 'output/03_deseq/PCA/PCA_time_{seq_batch}.pdf',
        PCA_treatment = 'output/03_deseq/PCA/PCA_treatment_{seq_batch}.pdf',
        PCA_timepoint_treatment = 'output/03_deseq/PCA/PCA_timepoint_treatment_{seq_batch}.pdf',
    log:
        'output/logs/deseq2_pca_{seq_batch}.log'
    singularity:
        bioconductor_container
    script:
        'src/deseq2_pca.R'

rule make_deseq2_object:
    input:
        salmon_quant = expand('output/02_salmon/quant_{sample}/quant.sf', sample=all_samples),
        gtf_file = 'data/misc_dmel-r6.63_files/dmel-all-r6.63.gtf',
        sample_data_file = 'data/sample_table.csv'
    output:
        dds_file = 'output/03_deseq/dds_both_batches.rds',
        dds_filtered_file = 'output/03_deseq/dds_both_batches_filtered.rds',
        dds_first_batch_file = 'output/03_deseq/dds_first_batch.rds',
        dds_first_batch_filtered_file = 'output/03_deseq/dds_first_batch_filtered.rds'
    log:
        'output/logs/make_deseq2_object.log'
    singularity:
        bioconductor_container
    script:
        'src/make_deseq2_object.R'

## gff has lines commented out with ## that cause reading into R to be difficult

####################################
## 02 - quasi-mapping with salmon ##
####################################

rule salmon_multiqc:
    input:
        expand('output/02_salmon/quant_{sample}/quant.sf', sample=all_samples)
    output:
        'output/02_salmon/multiqc/multiqc_report.html'
    params:
        analysis_dir = 'output/02_salmon',
        outdir = 'output/02_salmon/multiqc'
    log:
        'output/logs/salmon_multiqc.log'
    singularity:
        multiqc_container
    shell:
        'multiqc -f ' ##force to write over old output if it exists
        '{params.analysis_dir} '
        '-o {params.outdir} '
        '&> {log}'

## check auto library type
rule salmon_quant:
    input:
        r1_trimmed = 'output/01_read_prep/trimgalore/{sample}_r1.fq.gz',
        r2_trimmed = 'output/01_read_prep/trimgalore/{sample}_r2.fq.gz',
        index = 'output/02_salmon/index'
    output:
        'output/02_salmon/quant_{sample}/quant.sf'
    params:
        outdir = 'output/02_salmon/quant_{sample}'
    log:
        'output/logs/salmon/quant_{sample}.log'
    threads:
        50
    singularity:
        salmon_container
    shell:
        'salmon quant '
        '--index {input.index} '
        '--libType A ' ## salmon will infer library type itself - should be ISR based on Tom's code
        '-1 {input.r1_trimmed} '
        '-2 {input.r2_trimmed} '
        '-o {params.outdir} '
        '-p {threads} '
        '&> {log}'

## downloaded from current (v6.63) @: https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/index.html
rule salmon_index:
    input:
        'data/dmel-all-transcript-r6.63.fasta'
    output:
        'output/02_salmon/index'
    log:
        'output/logs/salmon/index.log'
    threads:
        50
    singularity:
        salmon_container
    shell:
        'salmon index '
        '--transcripts {input} '
        '--index {output} '
        '--threads {threads} '
        '&> {log}'

####################
## 01 - read prep ##
####################

rule read_prep_multiqc:
    input:
        untrimmed_fastqc = expand('output/01_read_prep/fastqc/untrimmed/{sample}_r{n}_fastqc.html', sample=all_samples, n=[1, 2]),
        trimmed_fastqc = expand('output/01_read_prep/fastqc/trimmed/{sample}_r{n}_fastqc.html', sample=all_samples, n=[1, 2]),
        trimming = expand('output/logs/trimgalore/{sample}.log', sample=all_samples)
    output:
        'output/01_read_prep/multiqc/multiqc_report.html'
    params:
        analysis_dir = 'output/01_read_prep',
        outdir = 'output/01_read_prep/multiqc'
    log:
        'output/logs/read_prep_multiqc.log'
    singularity:
        multiqc_container
    shell:
        'multiqc -f ' ##force to write over old output if it exists
        '{params.analysis_dir} '
        '-o {params.outdir} '
        '&> {log}'

#    reason: Missing output files: output/01_read_prep/fastqc/trimmed/4872etoh-2_r1_fastqc.html; Input files updated by another job: output/01_read_prep/trimgalore/4872etoh-2_r1.fq.gz

rule fastqc_posttrimming:
    input:
        'output/01_read_prep/trimgalore/{sample}_r{n}.fq.gz'
    output:
        'output/01_read_prep/fastqc/trimmed/{sample}_r{n}_fastqc.html'
    params:
        wd = 'output/01_read_prep/fastqc/trimmed'
    singularity:
        fastqc_container
    shell:
        'fastqc --outdir {params.wd} {input}'

rule rename_trimgalore_output:
    input:
        'output/01_read_prep/trimgalore/{sample}_r{n}_val_{n}.fq.gz'
    output:
        'output/01_read_prep/trimgalore/{sample}_r{n}.fq.gz'
    shell:
        'mv {input} {output}'

rule trimgalore:
    input:
        r1 = 'output/01_read_prep/joined/{sample}_r1.fq.gz',
        r2 = 'output/01_read_prep/joined/{sample}_r2.fq.gz'
    output:
        r1_trimmed = 'output/01_read_prep/trimgalore/{sample}_r1_val_1.fq.gz',
        r2_trimmed = 'output/01_read_prep/trimgalore/{sample}_r2_val_2.fq.gz'
    params:
        wd = 'output/01_read_prep/trimgalore',
        fastqc_wd = 'output/01_read_prep/trimgalore/fastqc'
    log:
        'output/logs/trimgalore/{sample}.log'
    threads:
        8
    singularity:
        trimgalore_container
    shell:
        'trim_galore '
        '--paired '
        '--cores {threads} '
        '--output_dir {params.wd} '
        '{input.r1} {input.r2} '
        '&> {log}'

rule fastqc_pretrimming:
    input:
        'output/01_read_prep/joined/{sample}_r{n}.fq.gz'
    output:
        'output/01_read_prep/fastqc/untrimmed/{sample}_r{n}_fastqc.html'
    params:
        wd = 'output/01_read_prep/fastqc/untrimmed'
    singularity:
        fastqc_container
    shell:
        'fastqc --outdir {params.wd} {input}'

rule join_reads:
    input:
        unpack(get_reads)
    output:
        r1 = 'output/01_read_prep/joined/{sample}_r1.fq.gz',
        r2 = 'output/01_read_prep/joined/{sample}_r2.fq.gz',
    shell:
        'cat {input.l1r1} {input.l2r1} > {output.r1} & '
        'cat {input.l1r2} {input.l2r2} > {output.r2} & '
        'wait'