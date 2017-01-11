percent.called = function(sample.dabg, cutoff){
	return(paste(round(100 * length(sample.dabg[sample.dabg < cutoff]) / length(sample.dabg), digits=2),"%",sep=""))
}#end def percent.called

param.table = read.table("parameters.txt", header=T, sep="\t")
output.folder = as.character(param.table$Value[param.table$Parameter == "Raw_Code_PC"])
dabg.folder = as.character(param.table$Value[param.table$Parameter == "DABG_Probeset_Folder"])
sample.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
call.rate.output = as.character(param.table$Value[param.table$Parameter == "Call_Rate_Table"])
pval.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "DABG_Pvalue_Cutoff"]))

normalized.table = read.table(paste(dabg.folder,"dabg.summary.txt",sep="/"), sep="\t", header=T)
probesetID = normalized.table$probeset_id
print(length(probesetID))

sample.table = read.table(sample.file, sep="\t", header=T)
mod.celFile = as.character(sample.table$Sample)
mod.celFile = paste("X",mod.celFile,sep="")
mod.celFile = gsub("-",".",mod.celFile)
mod.celFile = gsub("\\(",".",mod.celFile)
mod.celFile = gsub("\\)",".",mod.celFile)
sampleID = as.character(sample.table$ShortID)
dabg.mat = normalized.table[,match(mod.celFile, names(normalized.table))]
colnames(dabg.mat) = sampleID
print(dim(dabg.mat))

control.mat = dabg.mat[grep("^\\d+",probesetID,perl=T),]
control.call.rate = apply( dabg.mat, 2, percent.called, cutoff=pval.cutoff)
write.table(data.frame(sample=names(control.call.rate), call.rate = control.call.rate),call.rate.output, sep="\t", row.names=F, quote=F)