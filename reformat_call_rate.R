percent.called = function(sample.dabg, cutoff){
	return(paste(round(100 * length(sample.dabg[sample.dabg < cutoff]) / length(sample.dabg), digits=1),"%",sep=""))
}#end def percent.called

param.table = read.table("parameters.txt", header=T, sep="\t")
output.folder = as.character(param.table$Value[param.table$Parameter == "Raw_Code_PC"])
dabg.folder = as.character(param.table$Value[param.table$Parameter == "DABG_Probeset_Folder"])
sample.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])

setwd(output.folder)

stat.table = read.table(paste(dabg.folder,"dabg.report.txt",sep="/"), sep="\t", header=T)

sample.table = read.table(sample.file, sep="\t", header=T)
celFile = as.character(sample.table$Sample)
sampleID = as.character(sample.table$ShortID)
stat.table = stat.table[match(celFile,stat.table$cel_files),]

percent.called = paste(round(100*stat.table$all_probeset_percent_called, digits=2),"%",sep="")
write.table(data.frame(sample=sampleID, call.rate=percent.called), "probeset_call_rate.txt",sep="\t", quote=F, row.names=F)