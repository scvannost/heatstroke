library("DESeq2")
cnts = read.csv('/path/to/deseq/sac_control.csv', header=TRUE)
rownames(cnts) <- cnts$X
cnts = cnts[!(names(cnts) %in% c("X"))]

coldata = read.csv('/path/to/deseq/sac_control_coldata.csv', header=TRUE)

dds <- DESeqDataSetFromMatrix(countData=cnts, colData=coldata, design=~Sac)

dds <- DESeq(dds)
resultsNames(dds)

write.csv(results(dds, c('Sac','TcMax','1Day')),file='/path/to/deseq/sac_control_deseq.csv')
