from helpers import *

# load all the pickled data
samples = pickle.load(open('heat_stroke/data.p','rb'))
samples = samples[['Mouse','Time','Condition','Sac','Injection','[D-dimer]','PLT', 'GranzymeB', 'AvgDayTc',
            'AvgNightTc','CorrBW(g)','Dehydration(%)','Hem_WBC(10^9/L)','Hem_RBC(10^12/L)', 'Hem_HGB(g/dL)',
            'Hem_HCT(%)', 'Hem_MCV(fL)', 'Hem_MCH(pg)', 'Hem_MCHC(g/dL)', 'Hem_RDWc(%)', 'Hem_PLT(10^9/L)',
            'Hem_PCT(%)', 'Hem_MPV(fL)', 'Hem_PDWc(%)']]

activity = pickle.load(open('heat_stroke/activity.p','rb'))
temperature = pickle.load(open('heat_stroke/temperature.p','rb'))
hist = pickle.load(open('heat_stroke/hist.p','rb'))
rna_genes = pickle.load(open('heat_stroke/mapped/genes.p', 'rb'))
rna_genes = rna_genes[~rna_genes.index.duplicated(keep='first')]
rna_trans = pickle.load(open('heat_stroke/mapped/transcripts.p', 'rb'))

groups = samples['Condition'] + ' ' + samples['Sac']
groups_dict = dict(zip(list(set(groups)),range(4)))
samples['Group'] = [groups_dict[i] for i in groups]
print({'rgbm'[v]: k for k,v in groups_dict.items()})

# Getting rid of bad data
# 93, 94 had low read counts
# rest have oviAri3 contamination

idx = [i for i in rna_genes.columns if i not in [93,94,86,92,67,63]]
rna_genes = rna_genes[idx]
rna_trans = rna_trans[idx]
samples = samples.loc[idx]

## Log scale
# E = log2(TPM + 1)
rna_genes = np.log2(rna_genes + 1)
rna_trans = np.log2(rna_trans + 1)

# drop constant rows - std == 0 breaks zscore
rna_genes = rna_genes.loc[rna_genes.index[[min(rna_genes.loc[i]) != max(rna_genes.loc[i]) for i in rna_genes.index]]]
rna_trans = rna_trans.loc[rna_trans.index[[min(rna_trans.loc[i]) != max(rna_trans.loc[i]) for i in rna_trans.index]]]

## Zscore
rna_genes = pd.DataFrame(sp.stats.zscore(rna_genes,axis=1),index=rna_genes.index, columns=rna_genes.columns)
rna_trans = pd.DataFrame(sp.stats.zscore(rna_trans,axis=1),index=rna_trans.index, columns=rna_trans.columns)


print('All further analyses follow as in log2scaled.py')
