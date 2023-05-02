include: "rules/QC.smk"
# include: "rules/count.smk"

rule all:
    input:
        "results/raw_qc/multiqc_report.html",
        "results/trimmed_qc/multiqc_report.html",
        expand("results/summary_figures/{sb}_Read_Counts.pdf", sb=SAMPLE_BASES)