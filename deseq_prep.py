import pandas as pd
import pickle

samples = pickle.load(open('heat_stroke/data.p','rb'))
samples = samples[samples.columns[[2,3,4]]]

idx = [i for i in samples.index if i not in [93,94,86,92,67,63]]
samples=samples.loc[idx]

for i in range(1,49):
    genes = pd.read_csv('heat_stroke/mapped/'+str(i)+'.genes.results', sep='\t', usecols=[0,4], index_col=0)
    genes.columns = [i+48]
    
    genes = pd.DataFrame([int(i[0]) for i in genes.values],index=genes.index,columns=genes.columns)
    
    if i == 1:
        rna_genes = genes
    else:
        rna_genes = rna_genes.join(genes)
        
del genes

rna_genes.index = [gene_name(i) for i in rna_genes.index.values]
rna_genes = rna_genes[~rna_genes.index.duplicated(keep='first')]
rna_genes = rna_genes[idx]

with open('heat_stroke/counts.csv','w') as f:
    f.write(rna_genes.to_csv())
    
samples['Condition'] = [samples['Condition'][i] for i in samples.index]
samples['Injection'] = [i if i == 'Saline' else 'PolyIC' for i in samples['Injection']]
samples['Sac'] = [i if i == 'TcMax' else '1Day' for i in samples['Sac']]
with open('heat_stroke/counts_coldata.csv','w') as f:
    f.write(samples.to_csv())
