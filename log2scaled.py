from helpers import *

## Load all the pickled data
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

# add a groups column for easy identification
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

temperature = temperature[[i for i in temperature.columns.values if i in samples['Mouse'].tolist()]]
hist = hist[[i for i in hist.columns.values if i in samples['Mouse'].tolist()]]



## Log2 scaling
# E = log2(TPM + 1)
rna_genes = np.log2(rna_genes + 1)
rna_trans = np.log2(rna_trans + 1)


## Look at all measurements for simple, by eye, judgement
fig, ax = plt.subplots(int(len(samples.columns[5:-1])/3)+1,3,figsize=(15,40))
for n,i in enumerate(samples.columns[5:-1]):
    a = (n//3, n%3)
    ax[a].scatter(range(len(samples.index)),samples[i])
    ax[a].set_title('$'+i.replace(r'%',r'\%').replace(r'_',r'\_')+'$')
    ax[a].set_xlabel('Samples')
    
    if i == 'GranzymeB':
        ax[a].axhline(3000,c='r')
    elif i == 'AvgDayTc':
        ax[a].axhline(36.55,c='r')
    elif i == 'Dehydration(%)':
        ax[a].axhline(5.5,c='r')
plt.show()

## Show only the 3 I judged were separable by eye
fig, ax = plt.subplots(3,1,figsize=(5,15), sharex=True)
for n,i in enumerate(['GranzymeB','AvgDayTc','Dehydration(%)']):
    ax[n].scatter(range(len(samples.index)),samples[i])
    
ax[0].set_title('GranzymeB - 3000')
ax[1].set_title('AvgDayTc - 36.55')
ax[2].set_title('Dehydration(\%) - 5.5')

ax[0].axhline(3000,c='r')
ax[1].axhline(36.55,c='r')
ax[2].axhline(5.5,c='r')

ax[2].set_xlabel('Sample')
plt.show()




## Plot temperature data
fig, ax = plt.subplots(figsize=(15,7))

print({'rgbm'[v]:k for k,v in groups_dict.items()})
print('Dotted is Saline, solid is PolyI:C')

for i in temperature.columns.values:
    if i in samples['Mouse'].tolist():
        l = samples.iloc[samples['Mouse'].tolist().index(i)]
        ax.plot(temperature.index/60, temperature[i], c='rbgm'[l['Group']], linestyle=':' if l['Injection'] == 'Saline' else '-')

ax.set_xlabel('Time (hrs)')
ax.set_ylabel(r'Core Temperature $T_c$ (\textdegree C)')
ax.set_title('Body Temperature')
plt.show()




## Plot activity data
fig, ax = plt.subplots(figsize=(15,7))

print({'rgbm'[v]:k for k,v in groups_dict.items()})
print('Dotted is PolyI:C, solid is Saline')
for i in activity.columns.values:
    if i in samples['Mouse'].tolist():
        l = samples.iloc[samples['Mouse'].tolist().index(i)]
        ax.plot(activity.index/60, activity[i], c='rgbm'[l['Group']], linestyle='-' if l['Injection'] == 'Saline' else ':')

ax.set_xlabel('Time (hrs)')
ax.set_ylabel('Activity')
ax.set_title('Activity')
plt.show()



## Clustermap histology data
print({'rgbm'[v]:k for k,v in groups_dict.items()})
h = sns.clustermap(
    pd.DataFrame(np.array(hist),index=hist.index,
                 columns=[samples.iloc[samples['Mouse'].tolist().index(i)]['Injection'] for i in hist.columns]),
    row_cluster=False,
    col_colors=['rbgm'[samples.iloc[samples['Mouse'].tolist().index(i)]['Group']] for i in hist.columns])

plt.show()




## PCA for Transcripts
from sklearn.decomposition import PCA

gene_trans = PCA()
pc_trans = gene_trans.fit_transform(rna_trans.T)
fig, ax = plt.subplots(2,2,figsize=(15,15))
fig.suptitle('Log2-scaled Transcripts')

ax[0,0].scatter(pc_trans[:,0], pc_trans[:,1], c=['rgbm'[i] for i in samples['Group']])
ax[0,0].set_title(str({'rgbm'[v]: k for k,v in groups_dict.items()}))

ax[0,1].scatter(pc_trans[:,0], pc_trans[:,1], c=['b' if samples['Sac'][i] == '1 Day' else 'r' for i in samples.index])
ax[0,1].set_title('Sac - \'r\': TcMax, \'b\': Day 1')

ax[1,0].scatter(pc_trans[:,0], pc_trans[:,1], c=['b' if samples['Condition'][i] == 'Control' else 'r' for i in samples.index])
ax[1,0].set_title('Condition - \'r\': Heat, \'b\': Control')
ax[1,1].scatter(pc_trans[:,0], pc_trans[:,1], c=['b' if samples['Injection'][i] == 'Saline' else 'r' for i in samples.index])
ax[1,1].set_title('Injection - \'r\': Poly I:C, \'b\': Saline')

for i in [(0,0),(0,1),(1,0),(1,1)]:
    ax[i].set_xlabel('Principal Component 1 - '+str(np.round(100*gene_trans.explained_variance_ratio_[0],1))+r'\%')
    ax[i].set_ylabel('Principal Component 2 - '+str(np.round(100*gene_trans.explained_variance_ratio_[1],1))+r'\%')

for i in range(len(pc_trans[:,0])):
    ax[0,0].annotate(samples.index.values[i],xy=(pc_trans[i,0]+1,pc_trans[i,1]+1))
    ax[0,1].annotate(samples.index.values[i],xy=(pc_trans[i,0]+1,pc_trans[i,1]+1))
    ax[1,0].annotate(samples.index.values[i],xy=(pc_trans[i,0]+1,pc_trans[i,1]+1))
    ax[1,1].annotate(samples.index.values[i],xy=(pc_trans[i,0]+1,pc_trans[i,1]+1))
        
fig, ax = plt.subplots()
ax.bar(range(42),gene_trans.explained_variance_ratio_.tolist())
# ax.plot(range(42),[sum(gene_trans.explained_variance_ratio_[:i+1]) for i in range(42)],c='r')
ax.set_title('Log2-scaled Transcripts')
ax.set_xlabel('Principal Component')
ax.set_ylabel(r'Variance Explained (\%)')

plt.show()



## PCA for Genes
from sklearn.decomposition import PCA

gene_pca = PCA()
pc_genes = gene_pca.fit_transform(rna_genes.T)
fig, ax = plt.subplots(2,2,figsize=(15,15))
fig.suptitle('Log2-scaled Genes')

ax[0,0].scatter(pc_genes[:,0], pc_genes[:,1], c=['rgbm'[i] for i in samples['Group']])
ax[0,0].set_title(str({'rgbm'[v]: k for k,v in groups_dict.items()}))

ax[0,1].scatter(pc_genes[:,0], pc_genes[:,1], c=['b' if samples['Sac'][i] == '1 Day' else 'r' for i in samples.index])
ax[0,1].set_title('Sac - \'r\': TcMax, \'b\': Day 1')

ax[1,0].scatter(pc_genes[:,0], pc_genes[:,1], c=['b' if samples['Condition'][i] == 'Control' else 'r' for i in samples.index])
ax[1,0].set_title('Condition - \'r\': Heat, \'b\': Control')
ax[1,1].scatter(pc_genes[:,0], pc_genes[:,1], c=['b' if samples['Injection'][i] == 'Saline' else 'r' for i in samples.index])
ax[1,1].set_title('Injection - \'r\': Poly I:C, \'b\': Saline')

for i in [(0,0),(0,1),(1,0),(1,1)]:
    ax[i].set_xlabel('Principal Component 1 - '+str(np.round(100*gene_pca.explained_variance_ratio_[0],1))+r'\%')
    ax[i].set_ylabel('Principal Component 2 - '+str(np.round(100*gene_pca.explained_variance_ratio_[1],1))+r'\%')

for i in range(len(pc_genes[:,0])):
    ax[0,0].annotate(samples.index.values[i],xy=(pc_genes[i,0]+1,pc_genes[i,1]+1))
    ax[0,1].annotate(samples.index.values[i],xy=(pc_genes[i,0]+1,pc_genes[i,1]+1))
    ax[1,0].annotate(samples.index.values[i],xy=(pc_genes[i,0]+1,pc_genes[i,1]+1))
    ax[1,1].annotate(samples.index.values[i],xy=(pc_genes[i,0]+1,pc_genes[i,1]+1))

fig, ax = plt.subplots()
ax.bar(range(42),gene_pca.explained_variance_ratio_.tolist())
# ax.plot(range(42),[sum(gene_pca.explained_variance_ratio_[:i+1]) for i in range(42)],c='r')
ax.set_title('Log2-scaled Genes')
ax.set_xlabel('Principal Component')
ax.set_ylabel(r'Variance Explained (\%)')

plt.show()



## Statistics for what is separated by each PC from genes
num=3
fig, ax = plt.subplots(figsize=(12,1.25*num))
for i in range(num):
    print('PC'+str(i+1)+':')
    ax.scatter(pc_genes[:,i],[i+1]*pc_genes.shape[0],
               c=['b' if samples['Condition'][i] == 'Control' else 'r' for i in samples.index])
    p = sp.stats.fisher_exact(
        [[sum([pc_genes[n,i] > 0 for n,j in enumerate(samples.index) if samples['Condition'][j] == 'Control']),
          sum([pc_genes[n,i] < 0 for n,j in enumerate(samples.index) if samples['Condition'][j] == 'Control'])],         
        [sum([pc_genes[n,i] > 0 for n,j in enumerate(samples.index) if samples['Condition'][j] == 'Heat']),
         sum([pc_genes[n,i] < 0 for n,j in enumerate(samples.index) if samples['Condition'][j] == 'Heat'])]]
    )[1]
    print('Condition p =',p,bold = p<0.05/4)
    
    ax.scatter(pc_genes[:,i],[i+1.15]*pc_genes.shape[0],
               c=['k' if samples['Sac'][i] == '1 Day' else 'y' for i in samples.index])
    p = sp.stats.fisher_exact(
        [[sum([pc_genes[n,i] > 0 for n,j in enumerate(samples.index) if samples['Sac'][j] == '1 Day']),
          sum([pc_genes[n,i] < 0 for n,j in enumerate(samples.index) if samples['Sac'][j] == '1 Day'])],         
        [sum([pc_genes[n,i] > 0 for n,j in enumerate(samples.index) if samples['Sac'][j] == 'TcMax']),
         sum([pc_genes[n,i] < 0 for n,j in enumerate(samples.index) if samples['Sac'][j] == 'TcMax'])]]
    )[1]
    print('Sac p =',p,bold = p < 0.05/4)
    
    ax.scatter(pc_genes[:,i],[i+1.3]*pc_genes.shape[0],
               c=['c' if samples['Injection'][i] == 'Saline' else 'm' for i in samples.index])
    p = sp.stats.fisher_exact(
        [[sum([pc_genes[n,i] > 0 for n,j in enumerate(samples.index) if samples['Injection'][j] == 'Saline']),
          sum([pc_genes[n,i] < 0 for n,j in enumerate(samples.index) if samples['Injection'][j] == 'Saline'])],         
        [sum([pc_genes[n,i] > 0 for n,j in enumerate(samples.index) if samples['Injection'][j] != 'Saline']),
         sum([pc_genes[n,i] < 0 for n,j in enumerate(samples.index) if samples['Injection'][j] != 'Saline'])]]
    )[1]
    print('Injection p =',p,bold = p<0.05/4)
    
    p = sp.stats.fisher_exact(
        [[sum([pc_genes[n,i] > 0 for n,j in enumerate(samples.index) if samples['Group'][j] == 3]),
            sum(pc_genes[n,i] < 0 for n,j in enumerate(samples.index) if samples['Group'][j] == 3)],
         [sum([pc_genes[n,i] > 0 for n,j in enumerate(samples.index) if samples['Group'][j] != 3]),
            sum(pc_genes[n,i] < 0 for n,j in enumerate(samples.index) if samples['Group'][j] != 3)]
        ]
    )[1]
    print('Heat TcMax p =',p,'\n',bold=p<0.05/4)
        
ax.set_title('Log2-scaled: Red vs Blue = Heated vs Control;\nYellow vs Black = TcMax  vs 1 Day; Magenta vs Cyan = polyI:C vs Saline')
ax.set_xlabel('Value')
ax.axvline(0,c='k',linestyle='--')

fig.gca().invert_yaxis()
ax.set_ylabel('Principal Component')
ax.set_yticks(list(range(1,num+1)))

plt.show()




## Clustermapping of Genes
g = sns.clustermap(rna_genes.T, col_cluster=False,
              row_colors=['rgbm'[i] for i in samples['Group']])
g.fig.suptitle('Log2-scaled Genes')
plt.show()


## Kmeans clustering of Genes
# kmeans clustering scaled data
from sklearn.cluster import KMeans

silo = []
for i in range(2,21):
    kmeans = KMeans(n_clusters=i,max_iter=1000)
    lbls = kmeans.fit_predict(rna_genes.T)
    print(i)
    silo += [sk.metrics.silhouette_score(rna_genes.T,lbls)]
        
    if i == 2 or i == 3:
        fig, ax = plt.subplots()
        ax.scatter(pc_genes[:,0], pc_genes[:,1], c=['rgb'[i] for i in lbls])
        
        ax.set_xlabel('Principal Component 1 - '+str(np.round(100*gene_pca.explained_variance_ratio_[0],1))+r'\%')
        ax.set_ylabel('Principal Component 2 - '+str(np.round(100*gene_pca.explained_variance_ratio_[1],1))+r'\%')
        ax.set_title('Log2-scaled Genes - K-means clustering')
        fig.savefig('heat_stroke/log2genes_kmeans_pca'+str(i)+'.png')
        
    if i == 3:
        g = sns.clustermap(rna_genes.T,col_cluster=False,row_colors=['rgb'[i] for i in lbls])
        g.fig.suptitle('Log2-scaled Genes')
        g.savefig('heat_stroke/log2genes_kmeans_clustermap.png')
        
        

fig, ax = plt.subplots()
ax.plot(range(2,21),silo)
ax.set_xticks(list(range(2,21)))
ax.set_xlabel('k')
ax.set_ylabel('Silhouette Score')
ax.set_title('Log2-scaled Genes')
plt.show()



### LassoCV for Condition
## plotting alphas
from sklearn.linear_model import LassoCV
genes_tcmax = rna_genes[[i for i in samples.index if samples['Sac'][i] == 'TcMax']]
genes_oneday = rna_genes[[i for i in samples.index if samples['Sac'][i] == '1 Day']]

# vary these pairwise to get all possible combinations
base = genes_oneday
test = genes_tcmax

tcmax_condition = LassoCV(n_jobs=-1,cv=3)
tcmax_condition.fit(base.T,[1 if samples['Condition'][i] == 'Heat' else 0 for i in samples.index if i in base.columns])
print('Condition:',uline=True)
for i in range(10):
    plt.plot(tcmax_condition.alphas_,tcmax_condition.mse_path_[:,i])
plt.xlim([max(tcmax_condition.alphas_), min(tcmax_condition.alphas_)])
plt.axvline(tcmax_condition.alpha_,c='r')
plt.xlabel('Alpha')
plt.ylabel('CVE')
plt.title('TcMax Condition')
plt.show()

## Check Heat vs Control
p = np.array(list(zip(['Heat' if i > np.mean(tcmax_condition.predict(test.T)) else 'Control' for i in tcmax_condition.predict(test.T)], [samples['Condition'][i] for i in samples.index if i in test.columns])))
print(p)
p = sum([i[0] == i[1] for i in p])

## Check empirical p-value for Heat vs Control
labels = [1 if samples['Condition'][i] == 'Heat' else 0 for i in samples.index if i in base.columns]
dist = []
for i in range(50):
    np.random.shuffle(labels) # in place
    guess = LassoCV(n_jobs=-1,cv=10).fit(base.T,labels).predict(test.T)
    guess = [1 if i >= np.mean(guess) else 0 for i in guess]
    dist += [sum(labels == np.array(guess))]
    if i%10 == 9:
        print(str(i+1)+' done.')
    
fig,ax = plt.subplots()
ax.hist(dist)
plt.axvline(p,c='r')
for n,i in enumerate(sorted(dist, reverse=True)):
    if i <= p:
        break
        
p = n/len(dist)
print('p =',p)


## Check PolyIC vs Saline, with alphas graph
tcmax_injection = LassoCV(n_jobs=-1,cv=2)
tcmax_injection.fit(base.T, [1 if samples['Injection'][i] != 'Saline' else 0 for i in samples.index if i in base.columns])
print('Injection:',uline=True)
for i in range(2):
    plt.plot(tcmax_injection.alphas_,tcmax_injection.mse_path_[:,i])
plt.xlim([max(tcmax_injection.alphas_), min(tcmax_injection.alphas_)])
plt.axvline(tcmax_injection.alpha_,c='r')
plt.xlabel('Alpha')
plt.ylabel('CVE')
plt.title('TcMax Injection')
print(np.array(list(zip(['100 ug Poly I:C' if i > np.median(tcmax_injection.predict(test.T)) else 'Saline' for i in tcmax_injection.predict(test.T)], [samples['Injection'][i] for i in samples.index if i in test.columns]))))
plt.show()




## DESeq overview
de_sac = pd.read_csv('heat_stroke/DEsac.csv',index_col=0).fillna(
    {'baseMean':0,'log2FoldChange':0,'lfcSE':np.nan,'stat':np.inf,
    'pvalue':1,'padj':1}
)
de_injection = pd.read_csv('heat_stroke/DEinjection.csv',index_col=0).fillna(
    {'baseMean':0,'log2FoldChange':0,'lfcSE':np.nan,'stat':np.inf,
    'pvalue':1,'padj':1}
)
de_condition = pd.read_csv('heat_stroke/DEcondition.csv',index_col=0).fillna(
    {'baseMean':0,'log2FoldChange':0,'lfcSE':np.nan,'stat':np.inf,
    'pvalue':1,'padj':1}
)

sac_genes = de_sac.index[de_sac['padj'] < 0.05]
injection_genes = de_injection.index[de_injection['padj'] < 0.05]
condition_genes = de_condition.index[de_condition['padj'] < 0.05]

print('Sac:',sum(de_sac['padj'] < 0.05), str(np.round(100*sum(de_sac['padj'] < 0.05) / len(de_sac.index),4)) + '%')
print('Injection:',sum(de_injection['padj'] < 0.05), str(np.round(100*sum(de_injection['padj'] < 0.05) / len(de_injection.index),4)) + '%')
print('Condition:',sum(de_condition['padj'] < 0.05), str(np.round(100*sum(de_condition['padj'] < 0.05) / len(de_condition.index),4)) + '%')

# Venn diagram
from matplotlib_venn import venn3
abc = sum([i in sac_genes and i in condition_genes for i in injection_genes])
ab, ac = sum([i in sac_genes for i in injection_genes])-abc,sum([i in condition_genes for i in injection_genes])-abc
bc = sum([i in sac_genes for i in condition_genes])-abc

fig = plt.figure(figsize=(10,10))
venn3(subsets=(len(injection_genes)-ab-ac-abc, len(sac_genes)-ab-bc-abc, ab, len(condition_genes)-ac-bc-abc, ac, bc, abc),
     set_labels=('Injection','Sac','Condition'))
plt.show()



## DESeq clustermapping
rna_sac = rna_genes.loc[sac_genes]
rna_inj = rna_genes.loc[injection_genes]
rna_con = rna_genes.loc[condition_genes]

print({'rgbm'[v]:k for k,v in groups_dict.items()})
for n,i in enumerate([rna_sac, rna_inj, rna_con]):
    h = sns.clustermap(pd.DataFrame(np.array(i),
                                columns=[samples[['Sac', 'Injection', 'Condition'][n]][j] for j in i.columns.values],
                                index=i.index
                               ),
                   col_colors=['rgbm'[samples['Group'][j]] for j in i.columns.values]
                  )
    h.fig.suptitle(['Sac', 'Injection', 'Condition'][n])
plt.show()



## Setting up for Lasso
# use second line to make 3-way instead of 4-way

genes_tcmax = rna_genes[[i for i in samples.index if samples['Sac'][i] == 'TcMax']]
# genes_tcmax = genes_tcmax[genes_tcmax.columns[[samples.loc[i]['Injection'] == 'Saline' or samples.loc[i]['Condition'] == 'Heat' for i in genes_tcmax.columns]]]
genes_tcmax.loc['Group'] = [['Heat PolyI:C', 'Heat Saline', 'Control Saline', 'Control PolyI:C'].index(samples.loc[i]['Condition'] + ' ' + ('PolyI:C' if samples.loc[i]['Injection'] != 'Saline' else 'Saline')) for i in genes_tcmax.columns]

genes_oneday = rna_genes[[i for i in samples.index if samples['Sac'][i] == '1 Day']]
# genes_oneday = genes_oneday[genes_oneday.columns[[samples.loc[i]['Injection'] == 'Saline' or samples.loc[i]['Condition'] == 'Heat' for i in genes_oneday.columns]]]
genes_oneday.loc['Group'] = [['Heat PolyI:C', 'Heat Saline', 'Control Saline', 'Control PolyI:C'].index(samples.loc[i]['Condition']  + ' ' + ('PolyI:C' if samples.loc[i]['Injection'] != 'Saline' else 'Saline')) for i in genes_oneday.columns]


## run 3/4-way discriminating LassoCV
from sklearn.linear_model import LassoCV
tcmax_3way = LassoCV(n_jobs=-1,cv=5)
tcmax_x = genes_tcmax.loc[genes_tcmax.index[[not i is 'Group' for i in genes_tcmax.index]]]
tcmax_3way.fit(tcmax_x.T, genes_tcmax.loc['Group'])
print('TcMax:',uline=True)
for i in list(zip(np.round(tcmax_3way.predict(tcmax_x.T)), genes_tcmax.loc['Group'])):
    print(int(i[0]), int(i[1]), i[0] == i[1],ignore=True)
print()
    
day1_3way = LassoCV(n_jobs=-1,cv=3)
oneday_x = genes_oneday.loc[genes_oneday.index[[not i is 'Group' for i in genes_oneday.index]]]
day1_3way.fit(oneday_x.T, genes_oneday.loc['Group'])
print('One Day:',uline=True)
for i in list(zip(np.floor(tcmax_3way.predict(oneday_x.T)), genes_oneday.loc['Group'])):
    print(int(i[0]), int(i[1]), i[0] == i[1],ignore=True)
    
## Make empirical p-value for LogReg and plot
labels = np.array(genes_oneday.loc['Group'].tolist())
dist = []
for i in range(50):
    np.random.shuffle(labels) # in place
    lg = LogisticRegressionCV(n_jobs=-1,cv=3).fit(oneday_x.T, genes_oneday.loc['Group'])
    guess = lg.predict(oneday_x.T)
    
    dist += [sum(labels == guess)]
    
fig, ax = plt.subplots()
ax.hist(dist)
ax.axvline(len(labels),c='r') # H_a gets them all perfectly, so value = len(labels)
for n,i in enumerate(sorted(dist, reverse=True)):
    if i <= len(labels):
        break
        
print('p =',str(n/len(dist)))
plt.show()
