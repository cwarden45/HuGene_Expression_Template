parse.gene.info <- function(char.value)
{
	return.arr = c(NA, NA, NA, NA, NA, NA)
	if (char.value != "---"){
		per.gene = unlist(strsplit(char.value, split=" /// "))
		temp.info = unlist(strsplit(per.gene[1], split=" // "))
		
		temp.refseq = temp.info[1]
		temp.gene = temp.info[2]
		temp.description = temp.info[3]
		temp.cytoband = temp.info[4]
		temp.geneID = temp.info[5]
		
		genes = c()
		
		temp.num.genes = 1
		if(length(per.gene) > 1){
			for (i in 1:length(per.gene)){
				temp.info = unlist(strsplit(per.gene[i], split=" // "))
				temp.gene = temp.info[2]
				#print(temp.gene)
				
				if (!(temp.gene %in% genes)){
					genes = c(genes, temp.gene)
				} 
			}
			temp.num.genes = length(genes)
			temp.gene = paste(genes, collapse=" // ")
		}#end if(length(per.gene) > 1)
		
		return.arr = c(temp.refseq, temp.gene, temp.description, temp.cytoband, temp.geneID, temp.num.genes)
	}#end if (char.value != "---")
	return(return.arr)
}#end def parse.gene.info 

param.table = read.table("parameters.txt", header=T, sep="\t")
output.folder = as.character(param.table$Value[param.table$Parameter == "Raw_Code_PC"])
summary.type = as.character(param.table$Value[param.table$Parameter == "summary_method"])
transcript.annotation.file = as.character(param.table$Value[param.table$Parameter == "transcript_cluster_annotations"])
rma.folder = as.character(param.table$Value[param.table$Parameter == "RMA_Cluster_Folder"])
gene.mps.file = as.character(param.table$Value[param.table$Parameter == "MPS_Gene_File"])
sample.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
rma.table = as.character(param.table$Value[param.table$Parameter == "RMA_expression_file"])

setwd(output.folder)

normalized.table = read.table(paste(rma.folder, "rma.summary.txt",sep="/"), sep="\t", header=T)
probesetID = normalized.table$probeset_id
print(length(probesetID))

sample.table = read.table(sample.file, sep="\t", header=T)
mod.celFile = as.character(sample.table$Sample)
mod.celFile = paste("X",mod.celFile,sep="")
mod.celFile = gsub("-",".",mod.celFile)
mod.celFile = gsub("\\(",".",mod.celFile)
mod.celFile = gsub("\\)",".",mod.celFile)
sampleID = as.character(sample.table$ShortID)
rma.mat = normalized.table[,match(mod.celFile, names(normalized.table))]
colnames(rma.mat) = sampleID
print(dim(rma.mat))

if (summary.type == "gene_symbol"){
	probe.count.table = read.table(gene.mps.file, sep="\t", header=T)
	probe.count = probe.count.table$probe_count[match(probesetID, probe.count.table$probeset_id)]

	gene.per.probe.count = rep(1,length(probesetID))
	gene.per.probe.count[grep("^\\d",gene.per.probe.count,perl=T)]=NA
	
	annotated.expression = data.frame(transcript.cluster=probesetID, Symbol=probesetID,
										probe.count=probe.count,Num.Genes=gene.per.probe.count,
										rma.mat)
	write.table(annotated.expression, file = rma.table, sep="\t", row.names=F, quote=T)

}else if (summary.type == "transcript_cluster"){
	transcript.annotations = read.csv(transcript.annotation.file,comment.char="#", header=T)
	print(dim(transcript.annotations))
	transcript.annotations = transcript.annotations[match(probesetID, transcript.annotations$transcript_cluster_id),]

	target.chr = transcript.annotations$seqname
	target.start = transcript.annotations$start
	target.stop = transcript.annotations$stop
	target.strand = transcript.annotations$strand
	gene.text = as.character(transcript.annotations$gene_assignment)
	probe.count = transcript.annotations$total_probes
	probe.type = transcript.annotations$category

	gene.info = t(sapply(gene.text, parse.gene.info))
	#gene.info = matrix(unlist(sapply(gene.text, parse.gene.info)),ncol=6, byrow = T)
	colnames(gene.info) = c("Accession", "Symbol", "Name", "Cytoband", "geneID", "Num.Genes")

	annotated.expression = data.frame(transcript.cluster=probesetID, probe.type=probe.type,probe.count=probe.count,
										target.chr=target.chr, target.start=target.start, target.stop=target.stop, target.strand=target.strand,
										gene.info, rma.mat)
	write.table(annotated.expression, file = rma.table, sep="\t", row.names=F, quote=T)
}else{
	stop("'summary_method' must be 'gene_symbol' or 'transcript_cluster'")
}#end else
