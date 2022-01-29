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

for i in range(1,19):
    sample=f'211225-Mati_{str(i)}'
    
    # 1. Raw Fastq R1 and R2 counts
    raw_r1=fastqc_counts(f'/data/preprocess/results/raw_qc/{sample}_R1_001_fastqc.zip')
    raw_r2=fastqc_counts(f'/data/preprocess/results/raw_qc/{sample}_R2_001_fastqc.zip')
    read_counts.append({'Sample': sample, 'Read': 'R1', 'Step': 'Raw Fastq', 'Read Count': raw_r1})
    read_counts.append({'Sample': sample, 'Read': 'R2', 'Step': 'Raw Fastq', 'Read Count': raw_r2})
    
    # 2. Adapter/Qual trimmed R1 and R2 Counts
    trimmed_r1=fastqc_counts(f'/data/preprocess/results/trimmed_qc/{sample}_R1_001_sickle_cutadapt_trimmed_fastqc.zip')
    trimmed_r2=fastqc_counts(f'/data/preprocess/results/trimmed_qc/{sample}_R2_001_sickle_cutadapt_trimmed_fastqc.zip')
    read_counts.append({'Sample': sample, 'Read': 'R1', 'Step': 'Adapter/Qual Trimmed', 'Read Count': trimmed_r1})
    read_counts.append({'Sample': sample, 'Read': 'R2', 'Step': 'Adapter/Qual Trimmed', 'Read Count': trimmed_r2})
    
    # 3. Merged reads counts
    merged = read_fasta(f'/data/preprocess/results/merged/{sample}.fasta')
    merged_count=len(merged)
    read_counts.append({'Sample': sample, 'Read': 'Merged All', 'Step': 'Merged', 'Read Count': merged_count})
    
    # 4. Reads with correct length Counts
    correct_counts = len([s for s in merged if len(s)==392])
    read_counts.append({'Sample': sample, 'Read': 'Merged Correct Length', 'Step': 'Merged', 'Read Count': correct_counts})

for i in range(1,19):
    sample=f'211225-Mati_{str(i)}'
    sample_df = df.loc[df['Sample']==sample]
    
    plt.figure()
    sns.barplot(data=sample_df, x='Step', y='Read Count', hue='Read')
    plt.title(sample)
    plt.tight_layout()
    plt.savefig(f'./Read_Count_Figures/{sample}_Read_Counts.pdf')

for i in range(1,19):
    sample=f'211225-Mati_{str(i)}'
    with open(f'./20210521_Mati_N49P96VH-FW3/1_Count_Mutations_Output/Library_{sample}.txt', 'r') as file:
        text = file.read()
        d = ast.literal_eval(text.split('Number of barcodes per sequence: ')[1].split('\n#\n')[0])
        
    bc_count_data=[]
    for k in d:
        bc_count_data.append({'Sample': sample, 'Number of Barcodes': k, 'Count': d[k]})
        
    dff = pd.DataFrame(bc_count_data)
    
    plt.figure()
    sns.barplot(data=dff, x='Number of Barcodes', y='Count')
    plt.title(sample)
    plt.tight_layout()
    plt.savefig(f'./Number_of_BCs_per_Sequence/{sample}_BCs-per-seq.pdf')