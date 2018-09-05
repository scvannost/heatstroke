import json
import numpy as np #redundant from DataFrameN
import pandas as pd #redundant from DataFrameN
import pickle
from DataFrameN import *

for i in range(1,49):
    genes = pd.read_csv('heat_stroke/mapped/'+str(i)+'.genes.results', sep='\t', usecols=[0,5], index_col=0)
    genes.columns = [i+48]
    transcripts = pd.read_csv('heat_stroke/mapped/'+str(i)+'.isoforms.results', sep='\t', usecols=[0,5], index_col=0)
    transcripts.columns = [i+48]
    
    if i == 1:
        rna_genes = genes
        rna_trans = transcripts
    else:
        rna_genes = rna_genes.join(genes)
        rna_trans = rna_trans.join(transcripts)
        
del genes, transcripts

rna_genes.index = [gene_name(i) for i in rna_genes.index.values]
rna_trans.index = [isoform_name(i) for i in rna_trans.index.values]

pickle.dump(rna_genes, open('heat_stroke/mapped/genes.p', 'wb'))
pickle.dump(rna_trans, open('heat_stroke/mapped/transcripts.p', 'wb'))

samples = pd.read_csv('heat_stroke/sample_key.tsv',sep='\t',index_col=0,usecols=range(3))
samples = samples[samples['Organ'] == 'Kidney']
del samples['Organ']

with open('heat_stroke/mice.json','r') as f:
    mice = json.load(f)

data = np.zeros((2,2,3,3,10),dtype='object')
rows, cols, aisles, halls, = [], [], [], []
for m,(a,b) in enumerate(mice.items()):
    rows += [a]
    for n,(c,d) in enumerate(b.items()):
        if m == 0: cols += [c]
        for o,(e,f) in enumerate(d.items()):
            if m == 0 and n == 0: aisles += [e]
            for p,(g,h) in enumerate(f.items()):
                if m == 0 and n == 0 and o == 0: halls += [g]
                data[m,n,o,p,:] = h
mice = DataFrameN([rows,cols,aisles,halls],data) # see DataFrameN in github.com/scvannost/randompy
del data

labels = np.array([mice.labels(mice.index(i)) for i in samples['Mouse']]).T.tolist()
samples['Time']      = labels[0] #[i.split()[0] for i in labels[0]]
samples['Condition'] = labels[1]
samples['Sac']       = labels[2]
samples['Injection'] = labels[3]
del labels, mice

coag = pd.read_excel('heat_stroke/1408A Plasma Coagulation Markers (10-23-2017).xlsx',
                     sheet_name=1,index_col=0,skiprows=[0,2])
samples = pd.concat([
    samples,
    pd.DataFrame(
        [coag.loc[i].tolist()[2:8] if i in coag.index else [np.nan]*6 for i in samples['Mouse']]
        ,index=samples.index,columns=coag.columns[2:8])
    ],axis=1)
del coag

granz = pd.read_excel('heat_stroke/1408A_Collated Liver Granzyme B Data (10-23-2017).xlsx',
              index_col=0,skiprows=[0,1],usecols=[0,3])
samples['GranzymeB'] = [granz.loc[i]['Granzyme B (pg/ml)'] if i in granz.index else np.nan for i in samples['Mouse']]
del granz

heme = pd.read_excel('heat_stroke/1408A_Heat and Hematology (7-6-16).xlsx',
             sheet_name=2, index_col=0,usecols=[1]+list(range(7,39)))
heme.index = [str(i)[0:2] + '-' + str(i)[2:] for i in heme.index]
temp = heme.columns.values
heme = np.array([heme.loc[i].tolist() for i in samples['Mouse']]).T.tolist()

for n,i in enumerate(temp):
    samples[i] = heme[n]
del temp, heme

pickle.dump(samples, open('heat_stroke/data.p','wb'))

# load hourly activity and temperature data
activity = pd.read_excel('heat_stroke/1408A_Act Hour Avg Post Heat (6-29-2016).xlsx',
              sheet_name=1, usecols=range(350), skiprows=[1,2,3]+list(range(78,172)))
activity.columns = [i if type(i) is str else str(i)[0:2]+'-'+str(i)[2:] for i in activity.columns.values]
activity = activity[activity.columns[[i in samples['Mouse'].tolist() for i in activity.columns.values]]]

temperature = pd.read_excel('heat_stroke/1408A_Tc Hour Avg Post Heat (6-29-2016).xlsx',
             sheet_name=1,usecols=range(351),skiprows=[1,2,3]+list(range(78,172)))
temperature.columns = [i if type(i) is str else str(i)[0:2] + '-' + str(i)[2:] for i in temperature.columns.values]
temperature = temperature[temperature.columns[[i in samples['Mouse'].tolist() for i in temperature.columns.values]]]

pickle.dump(activity, open('heat_stroke/activity.p','wb'))
pickle.dump(temperature, open('heat_stroke/temperature.p','wb'))

hist = pd.read_excel('heat_stroke/organHistopathology_US Army 70464 jhh 4-26-17.xlsx', sheet_name=0,
             skiprows=[0,1,2,4,5]+list(range(18,24))+list(range(36,42))+list(range(54,60))+list(range(72,79))+list(range(91,96)),
             index_col=0,
             usecols=list(range(9))+list(range(14,20))+list(range(25,33))+list(range(36,46))+list(range(51,59))+list(range(64,72))+list(range(77,85))+list(range(90,98)))
hist.index = pd.MultiIndex.from_tuples(list(zip(['Liver']*12 + ['Kidney']*12 + ['Spleen']*12 + ['Lung']*12 + ['Duodenum']*12,hist.index)))
hist = hist[hist.columns[[i in samples['Mouse'].tolist() for i in hist.columns]]]

pickle.dump(hist, open('heat_stroke/hist.p','wb'))
