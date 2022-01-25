import os

configfile: "config.yaml"

SAMPLES = [f.split('.fastq.gz')[0] for f in os.listdir('data')]
SAMPLE_BASES = [s.split('_R')[0] for s in SAMPLES]

rule all:
    input:
        "results/raw_qc/multiqc_report.html",
        expand("results/trimmed/{sample}_sickle_cutadapt_trimmed.fastq.gz", sample=SAMPLES),
        "results/trimmed_qc/multiqc_report.html",
        expand("results/merged/{sb}.fasta", sb=SAMPLE_BASES)

rule raw_fastqc:
    input:
        "data/{sample}.fastq.gz"
    output:
        "results/raw_qc/{sample}_fastqc.html",
        "results/raw_qc/{sample}_fastqc.zip"
    threads: 16
    shell:
        "fastqc --noextract -o=results/raw_qc -t={threads} {input}"

rule raw_multiqc:
    input:
        expand("results/raw_qc/{sample}_fastqc.zip", sample=SAMPLES)
    output:
        "results/raw_qc/multiqc_report.html"
    threads: 16
    shell:
        "multiqc results/raw_qc/ -o results/raw_qc/"

rule cutadapt:
    input:
        "data/{sample}.fastq.gz",
        config['adapters']
    output:
        temp("results/trimmed/{sample}_cutadapt_trimmed.fastq.gz")
    threads: 16
    script:
        "scripts/adapter_trim.py"

rule sickle:
    input: 
        "results/trimmed/{s}_R1_001_cutadapt_trimmed.fastq.gz",
        "results/trimmed/{s}_R2_001_cutadapt_trimmed.fastq.gz"
    output:
        "results/trimmed/{s}_R1_001_sickle_cutadapt_trimmed.fastq.gz",
        "results/trimmed/{s}_R2_001_sickle_cutadapt_trimmed.fastq.gz"
    threads: 1
    shell:
        "sickle pe -t sanger -l 50 -q 20 -g -x -f {input[0]} -o {output[0]} -r {input[1]} -p {output[1]} -s /dev/null"

rule trimmed_fastqc:
    input:
        "results/trimmed/{sample}_sickle_cutadapt_trimmed.fastq.gz"
    output:
        "results/trimmed_qc/{sample}_sickle_cutadapt_trimmed_fastqc.html",
        "results/trimmed_qc/{sample}_sickle_cutadapt_trimmed_fastqc.zip"
    threads: 16
    shell:
        "fastqc --noextract -o=results/trimmed_qc -t={threads} {input}"

rule trimmed_multiqc:
    input:
        expand("results/trimmed_qc/{sample}_sickle_cutadapt_trimmed_fastqc.zip", sample=SAMPLES)
    output:
        "results/trimmed_qc/multiqc_report.html"
    threads: 16
    shell:
        "multiqc results/trimmed_qc/ -o results/trimmed_qc/"

rule pandaseq:
    input: 
        "results/trimmed/{m}_R1_001_sickle_cutadapt_trimmed.fastq.gz",
        "results/trimmed/{m}_R2_001_sickle_cutadapt_trimmed.fastq.gz"
    output:
        "results/merged/{m}.fasta
    threads: 16
    shell:
        "pandaseq -f {input[0]} -r {input[1]} -A simple_bayesian -d rbfkms -T {threads} -w {wildcards.m}.fasta"