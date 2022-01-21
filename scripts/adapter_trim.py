from Bio import SeqIO
from subprocess import Popen, PIPE

adapters = [str(s.seq) for s in SeqIO.parse(open(snakemake.input[1], 'r'), 'fasta')]
adapters = '-b ' + ' -b '.join(adapters)

cutadapt_cmd = f"cutadapt -o {snakemake.output[0]} {adapters} {snakemake.input[0]}"

p = Popen(cutadapt_cmd, stdout=PIPE, stderr=PIPE, shell=True)
stdout, stderr = p.communicate()