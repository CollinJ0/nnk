import os

configfile: "config.yaml"

SAMPLES = [f.split('.fastq.gz')[0] for f in os.listdir('data')]

rule all:
    input:
        "results/raw_qc/multiqc_report.html",
        expand("results/trimmed/{sample}_cutadapt_trimmed.fastq.gz", sample=SAMPLES)

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
        "results/trimmed/{sample}_cutadapt_trimmed.fastq.gz"
    threads: 16
    script:
        "scripts/adapter_trim.py"

