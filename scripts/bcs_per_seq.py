#open the Count output
with open(f'./1_Count_Mutations_Output/Library_{sample}.txt', 'r') as file:
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