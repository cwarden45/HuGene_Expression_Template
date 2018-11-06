probe.annotations = read.csv("/Volumes/user_data/Seq/cwarden/Array_Annotation_Files/HG-U133_Plus_2.na35.annot.csv", comment.char="#", header=T)
probesetID = as.character(probe.annotations$Probe.Set.ID)
gene = as.character(probe.annotations$Gene.Symbol)
probeCount = rep(11,length(probesetID))

gene.symbols = as.character(levels(as.factor(as.character(gene))))
probesets.per.gene = tapply(probeCount, gene, sum)
probeset.text = tapply(probesetID, gene, paste, collapse=" ")

#remove multi-mapped probes
probesets.per.gene = probesets.per.gene[-grep(" /// ", gene.symbols)]
probeset.text = probeset.text[-grep(" /// ", gene.symbols)]
gene.symbols = gene.symbols[-grep(" /// ", gene.symbols)]

#remove unannotated probes
probesets.per.gene = probesets.per.gene[-grep("---", gene.symbols)]
probeset.text = probeset.text[-grep("---", gene.symbols)]
gene.symbols = gene.symbols[-grep("---", gene.symbols)]

mps.table = data.frame(probeset_id=gene.symbols, transcript_cluster_id=gene.symbols,
						probeset_list = probeset.text,	probe_count=probesets.per.gene)
write.table(mps.table,"probeset_gene.mps",quote=F, row.names=F, sep="\t")