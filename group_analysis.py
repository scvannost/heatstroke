from helpers import *

## import the necessary data
samples = pickle.load(open('heat_stroke/data.p','rb'))
samples = samples[['Mouse','Time','Condition','Sac','Injection','[D-dimer]','PLT', 'GranzymeB', 'AvgDayTc',
            'AvgNightTc','CorrBW(g)','Dehydration(%)','Hem_WBC(10^9/L)','Hem_RBC(10^12/L)', 'Hem_HGB(g/dL)',
            'Hem_HCT(%)', 'Hem_MCV(fL)', 'Hem_MCH(pg)', 'Hem_MCHC(g/dL)', 'Hem_RDWc(%)', 'Hem_PLT(10^9/L)',
            'Hem_PCT(%)', 'Hem_MPV(fL)', 'Hem_PDWc(%)']]

rna_genes = pickle.load(open('heat_stroke/mapped/genes.p', 'rb'))
rna_genes = rna_genes[~rna_genes.index.duplicated(keep='first')]
rna_trans = pickle.load(open('heat_stroke/mapped/transcripts.p', 'rb'))

groups = samples['Condition'] + ' ' + samples['Sac']
groups_dict = dict(zip(list(set(groups)),range(8)))
samples['Group'] = [groups_dict[i] for i in groups]
print({'rgbm'[v]: k for k,v in groups_dict.items()})

# Getting rid of bad data
# 93, 94 had low read counts
# rest have oviAri3 contamination
idx = [i for i in rna_genes.columns if i not in [93,94,86,92,67,63]]
rna_genes = rna_genes[idx]
rna_trans = rna_trans[idx]
samples = samples.loc[idx]

# E = log2(TPM + 1)
rna_genes = np.log2(rna_genes + 1)
rna_trans = np.log2(rna_trans + 1)

# load the count file
counts = pd.read_csv('heat_stroke/counts.csv',header=0,index_col=0)
counts.columns = samples.index


### Begin group analysis
# make the files for DESeq
idx = samples.index[samples['Condition'] == 'Control'].tolist()
temp = counts[idx]
with open('heat_stroke/deseq/sac_control.csv','w') as f:
    f.write(temp.to_csv())
    
sac = samples['Sac'][idx]
sac = pd.Series(['TcMax' if i == 'TcMax' else '1Day' for i in sac], index=sac.index)

with open('heat_stroke/deseq/sac_control_coldata.csv','w') as f:
    f.writelines([
        ',"Sac"\n',
        sac.to_csv()
    ])

with open('heat_stroke/deseq/sac_control.cls','w') as f:
    f.writelines([
        str(len(idx)) +' 2 1\n',
        '# ' + ' '.join(list(set(sac))) + '\n',
        ' '.join(sac.tolist())
    ])
    
    
# Run deseq_groups.R here

ranked = pd.read_csv('heat_stroke/DEcondition.csv',index_col=0,usecols=[0,2,6])
padj = np.array(ranked['padj'].fillna(1))

check = [10**-((i+13)/10) for i in list(range(30)) + [40,50,60,70,80,90,100]]
plt.plot(check,[sum(padj < i) for i in check])
plt.axvline(0.05,c='r')
plt.annotate(str(sum(padj < 0.05)),xy=(0.045, sum(padj < 0.05)))
plt.axvline(0.001,c='r',linestyle=':')
plt.annotate(str(sum(padj < 0.001)), xy=(0.001, sum(padj < 0.001)))
plt.axhline(200,c='k',linestyle='--')
plt.title('Condition for all')

# make file for gsea
ranked.index = [i.upper() for i in ranked.index.values]
idx = padj < 0.05
with open('heat_stroke/deseq/condition_all_ranked.rnk','w') as f:
    f.write(ranked['log2FoldChange'][idx].fillna(0).to_csv(header=False,sep='\t'))
