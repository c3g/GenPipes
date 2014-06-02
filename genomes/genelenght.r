library(gqSeqUtils)
le = calculateGeneLengthsFromGtf("BosTau7.gtf")
le = data.frame(le,row.names=names(le))
write.table(le,file='gene_lengths.tsv',sep='\t',col.names=FALSE,row.names=TRUE,quote=FALSE)