import os

import ast
import gzip
import zipfile

import matplotlib.pyplot as plt
import seaborn as sns

from Bio import SeqIO

def read_fastq(fastq_file):
    with gzip.open(fastq_file, "rt") as f:
        sequences = [Sequence(s) for s in SeqIO.parse(f, 'fastq')]
    return sequences

def read_fasta(fasta_file):
    with open(fasta_file) as f:
        sequences = [Sequence(s) for s in SeqIO.parse(f, 'fasta')]
    return sequences

def fastqc_counts(zipped_file):
    archive = zipfile.ZipFile(zipped_file, 'r')
    imgdata = archive.read(f"{os.path.basename(zipped_file).split('.')[0]}/fastqc_data.txt")
    counts = int(str(imgdata).split('Total Sequences')[1].split('\\t')[1].split('\\n')[0])
    return counts


read_counts=[]

# 1. Raw Fastq R1 and R2 counts
raw_r1=fastqc_counts(snakemake.input[0])
raw_r2=fastqc_counts(snakemake.input[1])
read_counts.append({'Sample': snakemake.wildcards.s, 'Read': 'R1', 'Step': 'Raw Fastq', 'Read Count': raw_r1})
read_counts.append({'Sample': snakemake.wildcards.s, 'Read': 'R2', 'Step': 'Raw Fastq', 'Read Count': raw_r2})

# 2. Adapter/Qual trimmed R1 and R2 Counts
trimmed_r1=fastqc_counts(snakemake.input[2])
trimmed_r2=fastqc_counts(snakemake.input[3])
read_counts.append({'Sample': snakemake.wildcards.s, 'Read': 'R1', 'Step': 'Adapter/Qual Trimmed', 'Read Count': trimmed_r1})
read_counts.append({'Sample': snakemake.wildcards.s, 'Read': 'R2', 'Step': 'Adapter/Qual Trimmed', 'Read Count': trimmed_r2})
    
# 3. Merged reads counts
merged = read_fasta(snakemake.input[4])
merged_count=len(merged)
read_counts.append({'Sample': snakemake.wildcards.s, 'Read': 'Merged All', 'Step': 'Merged', 'Read Count': merged_count})

# 4. Reads with correct length Counts
# correct_counts = len([s for s in merged if len(s)==392])
# read_counts.append({'Sample': sample, 'Read': 'Merged Correct Length', 'Step': 'Merged', 'Read Count': correct_counts})

sample_df = pd.DataFrame(read_counts)

plt.figure()
sns.barplot(data=sample_df, x='Step', y='Read Count', hue='Read')
plt.title({snakemake.wildcards.s})
plt.tight_layout()
plt.savefig(f"results/summary_figures/{snakemake.wildcards.s}_Read_Counts.pdf")