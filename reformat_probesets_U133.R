normalized.table = read.table("custom_RMA/rma-sketch.summary.txt", sep="\t", header=T)
probesetID = normalized.table$probeset_id

rmaIDs = gsub(".CEL$","",names(normalized.table)[-1])
rmaIDs = gsub("^GSM\\d+_","",rmaIDs)
rmaIDs = gsub("_",".",rmaIDs)
rma.mat = normalized.table[,2:ncol(normalized.table)]
colnames(rma.mat) = rmaIDs

annotated.expression = data.frame(Gene=probesetID, rma.mat)
write.table(annotated.expression,"APT_Gene_RMA_expression.txt", sep="\t", row.names=F)