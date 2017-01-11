parse.gene.info <- function(char.value)
{
	return.arr = c(NA, NA, NA)
	if (char.value != "---"){
		per.gene = unlist(strsplit(char.value, split=" /// "))
		temp.info = unlist(strsplit(per.gene[1], split=" // "))
		
		temp.refseq = temp.info[1]
		temp.symbol = temp.info[2]
		
		genes = c()
		
		temp.num.genes = 1
		if(length(per.gene) > 1){
			for (i in 1:length(per.gene)){
				temp.info = unlist(strsplit(per.gene[i], split=" // "))
				temp.symbol = temp.info[2]
				#print(temp.symbol)
				
				if (!(temp.symbol %in% genes)){
					genes = c(genes, temp.symbol)
				} 
			}#end if(length(per.gene) > 1)
			temp.num.genes = length(genes)
			temp.symbol = paste(genes, collapse=" // ")
		}#end if(length(per.gene) > 1)
		
		return.arr = c(temp.refseq, temp.symbol, temp.num.genes)
	}#end if (char.value != "---")
	return(return.arr)
}#end def parse.gene.info

param.table = read.table("parameters.txt", header=T, sep="\t")
output.folder = as.character(param.table$Value[param.table$Parameter == "Raw_Code_PC"])
annotation.file = as.character(param.table$Value[param.table$Parameter == "probeset_gene_annotations"])
mps.file = as.character(param.table$Value[param.table$Parameter == "MPS_Gene_File"])

setwd(output.folder)

normalized.table = read.csv(annotation.file,comment.char="#", header=T)
probesetID = as.character(normalized.table$probeset_id)
gene = as.character(normalized.table$gene_assignment)
probeCount = normalized.table$probe_count

control.probeID = probesetID[gene == "---"]
control.probe.count = probeCount[gene == "---"]
control.probeID = control.probeID[control.probe.count > 4]
control.probe.count = control.probe.count[control.probe.count > 4]

gene.probeID = probesetID[gene != "---"]
gene.info = gene[gene != "---"]
gene.probe.count = probeCount[gene != "---"]

temp.gene.info = t(sapply(gene.info, parse.gene.info))

symbol = as.character(temp.gene.info[,2])

gene.symbols = as.character(levels(as.factor(as.character(symbol))))
probesets.per.gene = tapply(gene.probe.count, symbol, sum)
probeset.text = tapply(gene.probeID, symbol, paste, collapse=" ")

#just remove multi-mapped probes
probesets.per.gene = probesets.per.gene[-grep(" // ", gene.symbols)]
probeset.text = probeset.text[-grep(" // ", gene.symbols)]
gene.symbols = gene.symbols[-grep(" // ", gene.symbols)]

gene.symbols = c(gene.symbols, control.probeID)
probeset.text = c(probeset.text, control.probeID)
probesets.per.gene = c(probesets.per.gene, control.probe.count)

mps.table = data.frame(probeset_id=gene.symbols, transcript_cluster_id=gene.symbols,
						probeset_list = probeset.text,	probe_count=probesets.per.gene)
write.table(mps.table,mps.file,quote=F, row.names=F, sep="\t")