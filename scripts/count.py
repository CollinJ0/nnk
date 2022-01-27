import os, sys

from Bio import SeqIO
from Bio import pairwise2
import time
import itertools
import multiprocessing as mp
from abutils.utils.alignment import global_alignment, muscle
from abutils.core.sequence import Sequence, read_fasta

from collections import Counter

import pandas as pd
import numpy as np

import multiprocessing as mp
import subprocess as sp
import tempfile

from abutils.utils.jobs import monitor_mp_jobs
from abutils.utils.pipeline import list_files, make_dir
from abutils.utils.progbar import progress_bar

from natsort import natsorted
import re

full_ref=''
ref=''
excluded_positions=[]

#place fasta files in preferred directory and change code below
merged_fasta_files_dir='/data/preprocess/results/merged/'

#list of files
libraries = [os.path.join(merged_fasta_files_dir, l) for l in os.listdir(merged_fasta_files_dir)]

#get the actual variants from the Twist DNA variant file
variant_df = pd.read_csv('./N49P96_FW3_Variant_Sheet.csv')

# #make truncated ref start and end
truncated_offset=254
variant_df['Truncated Ref Start'] = variant_df['Start']-truncated_offset
variant_df['Truncated Ref End'] = variant_df['End']-truncated_offset

pos_start_row=0
pos_end_row=2488

output_dir = './1_Count_Mutations_Output/'


with open('./N49P7_FW3_NNK_Variants.txt', 'w') as vf:
    vf.write('\n'.join(variant_df['Variant'].tolist()))

# Now test with actual variants
def write_potential_variants_to_file(new_sequence, 
                                     variant_df, 
                                     temp_file, 
                                     pos_start_row=pos_start_row, 
                                     pos_end_row=pos_end_row):
    #length of variant barcodes
    variant_length = 21

    #length from beginning of NNK to end of variant
    downstream_of_nnk_length = 12

    #length from beginning to variant
    uptream_of_nnk_length = 9

#     #get first truncated variants
#     first_three_range = np.arange(0, 9, 3)

    #get full variants
    full_variant_range = np.unique(variant_df.loc[pos_start_row:pos_end_row]['Truncated Ref Start']).tolist()

#     #get last truncated variants
#     last_three_range = np.arange(435, 443, 3)
    
#     #write first three truncated to temp_file
#     for first_three in first_three_range:
#         temp_file.write(bytes(new_sequence[0:first_three+downstream_of_nnk_length] + '\n', encoding = 'utf-8'))
    
    #write full variants to file
    for full_variant_pos in full_variant_range:
        temp_file.write(bytes(new_sequence[full_variant_pos:full_variant_pos+variant_length] + '\n', encoding = 'utf-8'))
    
#     #write last three truncated variants to file
#     for i in last_three_range:
#         temp_file.write(bytes(new_sequence[i-uptream_of_nnk_length:] + '\n', encoding = 'utf-8'))

    temp_file.flush()
    return temp_file.name
    
    
#Find present > 1
def find_multiples(file, variant_file='./N49P7_FW3_NNK_Variants.txt', temp_dir='./tempdir/'):
    outs=[]
    cmd = "cat {} {} | sort -T {} | uniq -c | grep -v '^      1'".format(file, variant_file, temp_dir)
    p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, encoding='utf8')
    stdout, stderr = p.communicate()
    
    for line in stdout.split('\n'):
        if not line.strip():
            continue
        else:
            outs.append(line[6:])
#         if not line.strip():
#             continue
#         count = int(line.split()[0])
#         if count not in counts:
#             counts[count] = 1
#         else:
#             counts[count] += 1
    return outs
    
    
def analyze_seq(seq, variant_df=variant_df):
    with tempfile.NamedTemporaryFile() as temp:
        new_file = write_potential_variants_to_file(seq, variant_df, temp)
        _multis=find_multiples(new_file)
    return _multis

def print_splash(txt_file):
    # will print splash from text file made by the command
    #'$ figlet 'Sample Text' > Sample_text_file.txt 
    with open(txt_file, 'r') as f:
        lines = f.readlines()
    for l in lines[:1]:
        print(' ' + l, end=' ')
    for l in lines[1:]:
        print(l, end=' ')
        
def print_loading():
    for i in range(25):
        print('Loading Program{}'.format('.'*i), end='\r')
        time.sleep(0.2)
        
    print(' Loading Program{}Done!'.format('.'*25))


# In[9]:


if __name__ == "__main__":
    #Print the Banner
    for f in sorted([d for d in os.listdir('.') if 'banner' in d]):
        print_splash(f)
        time.sleep(0.3)
    
    #load the program
    print(' ')
    print_loading()
    print(' ')
    
    print('#######################')
    print('Starting Pipeline')
    print(' ')
    print('{} Fasta Files Selected:'.format(len(libraries)))
    for _lib in natsorted(libraries):
        print("  - Library_{}".format(_lib.split('/')[-1].split('.')[0]))
    print(' ')
    print('#######################')
    
    #go through each library
    for lib in natsorted(libraries):
        #get lib number
        lib_number = "Library_{}".format(lib.split('/')[-1].split('.')[0])
        print('-'*30)
        print(' ')
        print('Running {}'.format(lib_number))
        print(' ')
        print('-'*30)
        print(' ')
        print("Reading Sequences....", end=' ')
        #read in sequences and DO NOT get reverse comps
        sequences = read_fasta(lib)
        #sequences = [s.reverse_complement for s in sequences if len(s.sequence)==len(ref)]
        sequences = [s for s in sequences if len(s.sequence)==len(ref)]
        print("\rReading Sequences....Done! Number of Sequences: {}".format(str(len(sequences))))
        
        #initialize multiprocessing pool
        p = mp.Pool(maxtasksperchild=500)

        #initialize list to capture asynchronized results
        async_results = []
        
        print(' ')
        print('Analyzing Sequences...', end=' ')
        #analyze each sequence
        for s in sequences:
            async_results.append(p.apply_async(analyze_seq, args=(s, )))
        monitor_mp_jobs(async_results)
        
        #capture results
        results = [ar.get() for ar in async_results]
        
        #count number of times each BC was seen per sequence
        num_times_each_bc = str(dict(Counter(re.findall(r'\d+', ''.join([''.join(r) for r in results])))))
        
        #Count number of barcodes per sequence
        num_bcs_per_seq = str(dict(Counter([len(r) for r in results])))
        
        print('\rAnalyzing Sequences...Done!')
        print(' ')
        print('********************')
        print('Summary of Barcodes:')
        print('********************')
        print('Number of times seen each Barcode: {}\n'.format(num_times_each_bc))
        print('Number of barcodes per sequence: {}\n'.format(num_bcs_per_seq))
        print('********************')
        
        print(' ')
        print(' ')
        print('Writing Results to File...', end='\r')
        #write the results to a file
        with open(os.path.join(output_dir, lib_number+'.txt'), 'w') as out_file:
            out_file.write('#\n')
            out_file.write('Number of times seen each Barcode: {}\n'.format(num_times_each_bc))
            out_file.write('Number of barcodes per sequence: {}\n'.format(num_bcs_per_seq))
            out_file.write('#\n')
            out_file.write('\n'.join([r[0][1:] for r in results if len(r)==1]))
            
        print('Writing Results to File...Done!')
        print(' ')
        
        results=''
        
        print('-'*30)
        print(' ')
        
        #close multiprocessing
        p.terminate()
        p.join()