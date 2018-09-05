library("DESeq2")
cnts = read.csv('~/Documents/DAL/heat_stroke/counts.csv', header=TRUE)
rownames(cnts) <- cnts$X
cnts = cnts[2:43]

coldata = read.csv('~/Documents/DAL/heat_stroke/counts_coldata.csv', header=TRUE)
coldata = coldata[2:4]

dds <- DESeqDataSetFromMatrix(countData=cnts, colData=coldata, design=~Condition+Injection+Sac)

dds <- DESeq(dds)

write.csv(results(dds, c('Condition','Heat','Control')),file='~/Documents/DAL/heat_stroke/DEcondition.csv')
write.csv(results(dds, c('Sac','TcMax','1Day')), file='~/Documents/DAL/heat_stroke/DEsac.csv')
write.csv(results(dds, c('Injection','PolyIC','Saline')), file='~/Documents/DAL/heat_stroke/DEinjection.csv')
