normalizeTotalExpression = function (geneExpr, totalReads) {
	return(geneExpr / totalReads)
}#end def normalizeTotalExpression

avgGroupExpression = function (geneExpr, groups) {
	avg.expr = tapply(geneExpr, groups, mean)
	return(avg.expr)
}#end def avgGroupExpression

standardize.arr = function(arr)
{
	center.arr = as.numeric(arr) - mean(as.numeric(arr), na.rm=T)
	norm.arr = center.arr / sd(center.arr, na.rm=T)
	return(norm.arr)
}#end def standardize.arr

count.defined.values = function(arr, expr.cutoff)
{
	return(length(arr[arr > expr.cutoff]))
}#end def count.values

count.na.values = function(arr)
{
	return(length(arr[is.na(arr)]))
}#end def count.values

ratio2fc = function(value)
{
	if(value >= 0){
		return(2^value)
	} else {
		return (-2^(-value))
	}
}#end def ratio2fc

gene.lm = function(arr, var1, var2=c(), var3=c())
{	
	if (length(var2) == 0){
		fit = lm(as.numeric(arr) ~ var1)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
	} else if (length(var3) == 0){
		fit = lm(as.numeric(arr) ~ var1 + var2)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
	} else {
		fit = lm(as.numeric(arr) ~ var2*var3 + var2 + var3)
		result = summary(fit)
		pvalue = result$coefficients[4,4]
	}
	return(pvalue)
}#end def gene.lm

gene.aov = function(arr, var1, var2=c(), var3=c())
{	
	if (length(var2) == 0){
		fit = aov(as.numeric(arr) ~ var1)
		result = summary(fit)
		aov.pvalue = result[[1]][['Pr(>F)']][1]
	} else if (length(var3) == 0){
		fit = aov(as.numeric(arr) ~ var1 + var2)
		result = summary(fit)
		aov.pvalue = result[[1]][['Pr(>F)']][1]
	} else {
		fit = aov(as.numeric(arr) ~ var2*var3 + var2 + var3)
		result = summary(fit)
		aov.pvalue = result[[1]][['Pr(>F)']][3]
	}
	return(aov.pvalue)
}#end def gene.aov

calc.gene.cor = function(arr, indep.var)
{	
	na.count = length(arr[!is.na(arr)])
	if((na.count >= 3) & (sd(arr) != 0)){
		gene.cor.coef = cor(arr,indep.var)
	} else {
		gene.cor.coef = NA
	}
	return(gene.cor.coef)
}#end def calc.gene.cor

library(gplots)
fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black","darkgreen","thistle1","tan","orchid1",colors())

param.table = read.table("parameters.txt", header=T, sep="\t")
comp.name=as.character(param.table$Value[param.table$Parameter == "comp_name"])
detection.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "RMA_expression_cutoff"]))
min.fraction.expressed = as.numeric(as.character(param.table$Value[param.table$Parameter == "minimum_fraction_expressed"]))
fc.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "fold_change_cutoff"]))
cor.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "cor_cutoff"]))
cor.cutoff2 = as.numeric(as.character(param.table$Value[param.table$Parameter == "sec_cor_cutoff"]))
pvalue.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "pvalue_cutoff"]))
fdr.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "fdr_cutoff"]))
fc.cutoff2 = as.numeric(as.character(param.table$Value[param.table$Parameter == "sec_fold_change_cutoff"]))
pvalue.cutoff2 = as.numeric(as.character(param.table$Value[param.table$Parameter == "sec_pvalue_cutoff"]))
fdr.cutoff2 = as.numeric(as.character(param.table$Value[param.table$Parameter == "sec_fdr_cutoff"]))
deg.groups = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "deg_groups"]), split=","))
plot.groups = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "plot_groups"]), split=","))
trt.group = as.character(param.table$Value[param.table$Parameter == "treatment_group"])
trt.group2 = as.character(param.table$Value[param.table$Parameter == "secondary_trt"])
interaction.flag = as.character(param.table$Value[param.table$Parameter == "interaction"])
interaction.flag[interaction.flag == "none"]="no"
pvalue.method = as.character(param.table$Value[param.table$Parameter == "pvalue_method"])
fdr.method = as.character(param.table$Value[param.table$Parameter == "fdr_method"])
output.folder = as.character(param.table$Value[param.table$Parameter == "Raw_Code_PC"])
user.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
rma.file = as.character(param.table$Value[param.table$Parameter == "RMA_expression_file"])
gostat.type = as.character(param.table$Value[param.table$Parameter == "R_GO_type"])
run.gostat = as.character(param.table$Value[param.table$Parameter == "R_GO"])

setwd(output.folder)

sample.description.table = read.table(sample.description.file, sep="\t", header=T)
sample.label = sample.description.table$ShortID

deg.group.table = sample.description.table[,deg.groups]
if (length(deg.groups) == 1){
	deg.meta = sample.description.table[!is.na(deg.group.table),]
} else {
	deg.grp.na.counts = apply(deg.group.table, 1, count.na.values)
	deg.meta = sample.description.table[deg.grp.na.counts == 0,]
}

rma.table = read.table(rma.file, head=T, sep="\t")
rma.mat= rma.table[,match(sample.label,names(rma.table))]

print(dim(rma.mat))
expressed.sample.count = apply(rma.mat, 1, count.defined.values, expr.cutoff = detection.cutoff)
rma.table = rma.table[expressed.sample.count >= round(min.fraction.expressed * ncol(rma.mat)),]
rma.mat = rma.mat[expressed.sample.count >= round(min.fraction.expressed * ncol(rma.mat)),]
print(dim(rma.mat))

transcript.cluster = rma.table$transcript.cluster
genes = rma.table$Symbol
probes.per.cluster= rma.table$probe.count
genes.per.probe = rma.table$Num.Genes
colnames(rma.mat) = as.character(sample.label)
rownames(rma.mat) = as.character(transcript.cluster)

if(length(plot.groups) == 1){
	print("Averaging Expression for One Variable (for plot.groups)")
	grp = sample.description.table[,plot.groups]
} else if ((length(plot.groups) == 2)&(interaction.flag == "no")){
	print("Averaging Expression for First Variable (for plot.groups)")
	grp = sample.description.table[,plot.groups[1]]
} else if (length(plot.groups) == 2){
	print("Averaging Expression for Interaction Variable (for plot.groups)")
	grp = paste(sample.description.table[,plot.groups[1]],sample.description.table[,plot.groups[2]],sep=":")
} else {
	stop("Code only compatible with 2 variables (with or without a 3rd interaction variable")
}

groupIDs = as.character(levels(as.factor(grp)))
average.rma = data.frame(t(apply(rma.mat, 1, avgGroupExpression, groups = grp)))
if(length(groupIDs) == 1){
	average.rma = t(average.rma)
} else {
	average.rma = average.rma
}
colnames(average.rma) = paste("avg.rma", sub("-",".",groupIDs), sep=".")

#remove undefined group IDs (so, you can visualize samples that aren't in your comparison)
if(length(deg.groups) == 1){
	var1 = sample.description.table[,deg.groups]
	deg.RMA = rma.mat[,!is.na(var1)]
	print(dim(deg.RMA))
	var1 = var1[!is.na(var1)]
	if (trt.group != "continuous"){
		var1 = as.factor(as.character(var1[!is.na(var1)]))
	}
} else if (length(deg.groups) == 2){
	if(interaction.flag == "filter-overlap"){
		var1 = sample.description.table[,deg.groups[1]]
		var2 = sample.description.table[,deg.groups[2]]
		deg.samples = !is.na(var1)&!is.na(var2)
		deg.RMA = rma.mat[,deg.samples]
		var1 = var1[deg.samples]
		var2 = var2[deg.samples]
			
		prim.deg.grp = var1[var2 == trt.group2]
		prim.RMA = rma.mat[,var2 == trt.group2]
		print(dim(prim.RMA))

		sec.deg.grp = var1[var2 != trt.group2]
		sec.RMA = rma.mat[,var2 != trt.group2]
		print(dim(sec.RMA))

	} else{
		var1 = sample.description.table[,deg.groups[1]]
		var2 = sample.description.table[,deg.groups[2]]
		deg.samples = !is.na(var1)&!is.na(var2)
		deg.RMA = rma.mat[,deg.samples]
		print(dim(deg.RMA))
		var1 = var1[deg.samples]
		if (trt.group != "continuous"){
			var1 = as.factor(as.character(var1[!is.na(var1)]))
		}
		var2 = var2[deg.samples]
		if (trt.group2 != "continuous"){
			var2 = as.factor(as.character(var2[!is.na(var2)]))
		}
	}
} else {
	stop("Code currently doesn't support more than 2 group model for DEG (with or without interaction)")
}

if(length(deg.groups) == 1){
	print("Averaging Expression for One Variable (for deg.groups)")
	contrast.grp = var1
} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
	print("Averaging Expression for First Variable (for deg.groups)")
	contrast.grp = var1
} else if (length(deg.groups) == 2){
	print("Averaging Expression for Interaction Variable (for deg.groups)")
	contrast.grp = paste(var1,var2,sep=":")
} else {
	stop("Code only compatible with 2 variables (with or without a 3rd interaction variable")
}

if (trt.group == "continuous"){
	contrast.grp = as.numeric(contrast.grp)
	
	gene.cor = apply(deg.RMA, 1, calc.gene.cor, indep.var=contrast.grp)
	
	fc.table = data.frame(cor=gene.cor)
} else {
	groupIDs = as.character(levels(as.factor(contrast.grp)))
	contrast.rma = data.frame(t(apply(deg.RMA, 1, avgGroupExpression, groups = contrast.grp)))
	colnames(contrast.rma) = paste("avg.log2.rma", sub("-",".",groupIDs), sep=".")
}#end else

if((interaction.flag == "no") & (trt.group != "continuous")){
	print("Calculating fold-change for primary variable")
	trt.expr = contrast.rma[,paste("avg.log2.rma", sub("-",".",trt.group), sep=".")]
	cntl.expr = contrast.rma[,paste("avg.log2.rma", sub("-",".",groupIDs[groupIDs != trt.group]), sep=".")]

	log2ratio = round(trt.expr - cntl.expr, digits = 2)
	fc = round(sapply(log2ratio, ratio2fc), digits = 2)
	fc.table = data.frame(log2ratio=log2ratio, fold.change=fc)
} else if ((interaction.flag == "model")|(interaction.flag == "filter-overlap")){
	if ((trt.group == "continuous")&(trt.group2 == "continuous")){
		print("Calculating  correlation for secondary variable")
		sec.contrast.grp = as.numeric(var2)
		
		gene.cor2 = apply(deg.RMA, 1, calc.gene.cor, indep.var=sec.contrast.grp)
		
		fc.table = data.frame(prim.cor=gene.cor, sec.cor = gene.cor2)
	} else if (trt.group == "continuous"){
		print("Fold-change / correlation cutoff not used for mixed variable analysis")
		print("NOTE: 'Up-Regulated' R output refers to genes that vary with FDR and p-value cutoffs")
		print("However, fold-change / correlation values for each separate variable are still provided")

		sec.groupIDs = var2
		sec.groups = as.character(levels(as.factor(sec.groupIDs)))
		sec.contrast.rma = data.frame(t(apply(deg.RMA, 1, avgGroupExpression, groups = sec.groupIDs)))
		colnames(sec.contrast.rma) = paste("avg.log2.rma", sub("-",".",groupIDs), sep=".")
		sec.trt.expr = sec.contrast.rma[,paste("avg.log2.rma", sub("-",".",trt.group2), sep=".")]
		sec.cntl.expr = sec.contrast.rma[,paste("avg.log2.rma", sub("-",".",sec.groups[sec.groups != trt.group2]), sep=".")]

		sec.log2ratio = round(sec.trt.expr - sec.cntl.expr, digits = 2)
		sec.fc = round(sapply(sec.log2ratio, ratio2fc), digits = 2)
		
		fc.table = data.frame(prim.cor=gene.cor, sec.fc = sec.fc)
	} else if (trt.group2 == "continuous"){	
		print("Fold-change / correlation cutoff not used for mixed variable analysis")
		print("NOTE: 'Up-Regulated' R output refers to genes that vary with FDR and p-value cutoffs")
		print("However, fold-change / correlation values for each separate variable are still provided")

		prim.groupIDs = var1
		prim.groups = as.character(levels(as.factor(prim.groupIDs)))
		prim.contrast.rma = data.frame(t(apply(deg.RMA, 1, avgGroupExpression, groups = prim.groupIDs)))
		colnames(prim.contrast.rma) = paste("avg.log2.rma", sub("-",".",prim.groups), sep=".")
		prim.trt = trt.group
		prim.cntl = prim.groups[prim.groups != trt.group]
		prim.trt.expr = prim.contrast.rma[,paste("avg.log2.rma", sub("-",".",prim.trt), sep=".")]
		prim.cntl.expr = prim.contrast.rma[,paste("avg.log2.rma", sub("-",".",prim.cntl), sep=".")]

		prim.log2ratio = round(prim.trt.expr - prim.cntl.expr, digits = 2)
		prim.fc = round(sapply(prim.log2ratio, ratio2fc), digits = 2)
		
		sec.contrast.grp = as.numeric(var2)
		gene.cor2 = apply(deg.RMA, 1, calc.gene.cor, indep.var=sec.contrast.grp)
		
		fc.table = data.frame(prim.fc=prim.fc, sec.cor = gene.cor2)
	} else {
		print("Calculating fold-change table for primary variables (within subsets of secondary variable)")
		prim.groups = paste(var1,var2,sep=":")
		prim.trt = paste(trt.group,trt.group2,sep=":")
		prim.cntl = paste(prim.groups[prim.groups != trt.group],trt.group2,sep=":")
		prim.trt.expr = contrast.rma[,paste("avg.log2.rma", sub("-",".",prim.trt), sep=".")]
		prim.cntl.expr = contrast.rma[,paste("avg.log2.rma", sub("-",".",prim.cntl), sep=".")]

		prim.log2ratio = round(prim.trt.expr - prim.cntl.expr, digits = 2)
		prim.fc = round(sapply(prim.log2ratio, ratio2fc), digits = 2)

		sec.groups = as.character(levels(as.factor(sample.description.table[,deg.groups[2]])))
		sec.trt = paste(trt.group, sec.groups[sec.groups != trt.group2], sep=":")
		sec.cntl = paste(prim.groups[prim.groups != trt.group], sec.groups[sec.groups != trt.group2], sep=":")
		sec.trt.expr = contrast.rma[,paste("avg.log2.rma", sub("-",".",sec.trt), sep=".")]
		sec.cntl.expr = contrast.rma[,paste("avg.log2.rma", sub("-",".",sec.cntl), sep=".")]

		sec.log2ratio = round(sec.trt.expr - sec.cntl.expr, digits = 2)
		sec.fc = round(sapply(sec.log2ratio, ratio2fc), digits = 2)

		overall.log2ratio = prim.log2ratio - sec.log2ratio
		overall.fc = round(sapply(overall.log2ratio, ratio2fc), digits = 2)

		fc.table = data.frame(fc1 = prim.fc, fc2=sec.fc, fc3=overall.fc)
		colnames(fc.table) = c(paste("fold.change",trt.group,":",trt.group2,sep="."),
								paste("fold.change",trt.group,":",sec.groups[sec.groups != trt.group2], sep="."),
								"overall.fold.change")
	}#end else
}else if(trt.group == "continuous"){
	print("Skipping fold-change calculation for continuous variable")
}else{
		stop("interaction must be \"no\", \"model\", or \"filter-overlap\"")
}#end else

rep.check = 1
for (i in 1:length(deg.groups)){
	deg.group = deg.groups[i]
	
	if((i == 1) & (trt.group != "continuous")){
		deg.group.values = as.factor(as.character(deg.meta[,deg.group]))
		min.reps = min(table(deg.group.values))
		if (min.reps < 2){
			rep.check=0
			print("There are not at least 2 samples per-group in order to calculate p-value.")
			print("In the future, please make sure you at least have duplicate samples.")
		}#end if (min.reps < 2)
	} else if ((i == 2) & (trt.group2 != "continuous")){
		deg.group.values = as.factor(as.character(deg.meta[,deg.group]))
		min.reps = min(table(deg.group.values))
		if (min.reps < 2){
			rep.check=0
			print("There are not at least 2 samples per-group in order to calculate p-value.")
			print("In the future, please make sure you at least have duplicate samples.")
		}#end if (min.reps < 2)
	} else if (i > 2){
		stop("Workflow currently doesn't support use of more than 2 variables")
	}
}#end for (deg.group in deg.groups)

if(rep.check == 1){
	#start p-value calculation
	if (pvalue.method == "limma"){
		library(limma)
		
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("limma, Two-Step Analysis")

			if (trt.group == "continuous"){
				prim.deg.grp = as.numeric(prim.deg.grp)
			}
			
			design <- model.matrix(~prim.deg.grp)
			fit <- lmFit(prim.RMA,design)
			fit <- eBayes(fit)
			pvalue.mat = data.frame(fit$p.value)
			prim.pvalue = pvalue.mat[,2]
			

			if (trt.group2 == "continuous"){
				sec.deg.grp = as.numeric(sec.deg.grp)
			}
			
			design <- model.matrix(~sec.deg.grp)
			fit <- lmFit(sec.RMA,design)
			fit <- eBayes(fit)
			pvalue.mat = data.frame(fit$p.value)
			sec.pvalue = pvalue.mat[,2]
			
		} else {
			if (length(deg.groups) == 1){
				print("limma with 1 variable")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}
				design <- model.matrix(~var1)
				fit <- lmFit(deg.RMA,design)
				fit <- eBayes(fit)
				pvalue.mat = data.frame(fit$p.value)
				test.pvalue = pvalue.mat[,2]
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("limma-voom with 2 variables")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				} else{
					var1 = as.factor(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				} else{
					var2 = as.factor(var2)
				}
				design <- model.matrix(~var1 + var2)
				fit <- lmFit(deg.RMA,design)
				fit <- eBayes(fit)
				pvalue.mat = data.frame(fit$p.value)
				test.pvalue = pvalue.mat[,2]
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("limma with 2 variables plus interaction")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}
				design <- model.matrix(~var1*var2 + var1 + var2)
				fit <- lmFit(deg.RMA,design)
				fit <- eBayes(fit)
				pvalue.mat = data.frame(fit$p.value)
				test.pvalue = pvalue.mat[,4]
			}
		}#end else
	} else if (pvalue.method == "lm"){
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("RMA linear regression, Two-Step Analysis")

			if (trt.group == "continuous"){
				prim.deg.grp = as.numeric(prim.deg.grp)
			}
			
			prim.pvalue = apply(prim.RMA, 1, gene.lm, var1=prim.deg.grp)
			

			if (trt.group2 == "continuous"){
				sec.deg.grp = as.numeric(sec.deg.grp)
			}
			
			sec.pvalue = apply(sec.RMA, 1, gene.lm, var1=sec.deg.grp)
		} else {
			if (length(deg.groups) == 1){
				print("RMA linear regression with 1 variable")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}
				test.pvalue = apply(deg.RMA, 1, gene.lm, var1=var1)
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("RMA linear regression with 2 variables")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}
				test.pvalue = apply(deg.RMA, 1, gene.lm, var1=var1, var2=var2)
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("RMA linear regression with 2 variables plus interaction")
				var3 = as.factor(paste(var1,var2,sep=":"))

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group == "continuous"){
					var2 = as.numeric(var2)
				}
				test.pvalue = apply(deg.RMA, 1, gene.lm, var1=var3, var2=var1, var3=var2)
			}
		}#end else
	} else if (pvalue.method == "ANOVA"){
		if ((length(deg.groups) == 2)&(interaction.flag == "filter-overlap")){
			print("RMA ANOVA, Two-Step Analysis")

			if (trt.group == "continuous"){
				prim.deg.grp = as.numeric(prim.deg.grp)
			}
			
			prim.pvalue = apply(prim.RMA, 1, gene.aov, var1=prim.deg.grp)
			

			if (trt.group2 == "continuous"){
				sec.deg.grp = as.numeric(sec.deg.grp)
			}
			
			sec.pvalue = apply(sec.RMA, 1, gene.aov, var1=sec.deg.grp)
		} else {
			if (length(deg.groups) == 1){
				print("RMA ANOVA with 1 variable")
				
				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}
				test.pvalue = apply(deg.RMA, 1, gene.aov, var1=var1)
			} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
				print("RMA ANOVA with 2 variables")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}
				test.pvalue = apply(deg.RMA, 1, gene.aov, var1=var1, var2=var2)
			} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
				print("RMA ANOVA with 2 variables plus interaction")
				var3 = as.factor(paste(var1,var2,sep=":"))

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}
				test.pvalue = apply(deg.RMA, 1, gene.aov, var1=var3, var2=var1, var3=var2)
			}
		}#end else
	} else{
		stop("pvalue_method must be \"limma\", \"lm\", or \"ANOVA\"")
	}
} else{
	test.pvalue = rep(1,times=length(genes))
	prim.pvalue = rep(1,times=length(genes))
	sec.pvalue = rep(1,times=length(genes))
}#end else

if (trt.group == "continuous"){
	upID = "Increased Expression"
	downID = "Decreased Expression"
} else {
	upID = paste(trt.group," Up",sep="")
	downID = paste(trt.group," Down",sep="")	
}


if (interaction.flag == "no"){
	if (fdr.method == "BH"){
		fdr = p.adjust(test.pvalue, "fdr")
	} else if (fdr.method == "q-value"){
		library(qvalue)
		qobj <- qvalue(p = test.pvalue)
		fdr = qobj$qvalue
		png(paste(pvalue.method,"_qvalue_plot.png",sep=""))
		qHist = hist(qobj)
		print(qHist)
		dev.off()
	} else if (fdr.method == "q-lfdr"){
		library(qvalue)
		qobj <- qvalue(p = test.pvalue)
		fdr = qobj$lfdr
		png(paste(pvalue.method,"_qvalue_plot.png",sep=""))
		qHist = hist(qobj)
		print(qHist)
		dev.off()
	} else {
		stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
	}
	status = rep("No Change", times=length(fdr))
	if (trt.group == "continuous"){
		status[(gene.cor >= cor.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
		status[(gene.cor <= -cor.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
		pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
	} else{
		status[(fc >= fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
		status[(fc <= -fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
		pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
	}#end else
} else{
	trt.group = prim.trt
	if(interaction.flag == "model"){
		if (fdr.method == "BH"){
			fdr = p.adjust(test.pvalue, "fdr")
		} else if (fdr.method == "q-value"){
			library(qvalue)
			qobj <- qvalue(p = test.pvalue)
			fdr = qobj$qvalue
			png(paste(pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else if (fdr.method == "q-lfdr"){
			library(qvalue)
			qobj <- qvalue(p = test.pvalue)
			fdr = qobj$lfdr
			png(paste(pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else {
			stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
		}
		status = rep("No Change", times=length(fdr))
		if ((trt.group == "continuous")&(trt.group2 == "continuous")){
			status[(gene.cor.int >= cor.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
			status[(gene.cor.int <= -cor.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
			pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
		} else if ((trt.group != "continuous")&(trt.group2 != "continuous")){
			status[(overall.fc >= fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
			status[(overall.fc <= -fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
			pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
		} else {
			upID = "Variable Expression"
			status[(test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
			pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
		}#end else
	} else if (interaction.flag == "filter-overlap"){
		if (fdr.method == "BH"){
			fdr = p.adjust(prim.pvalue, "fdr")
		} else if (fdr.method == "q-value"){
			library(qvalue)
			qobj <- qvalue(p = prim.pvalue)
			fdr = qobj$qvalue
			png(paste(pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else if (fdr.method == "q-lfdr"){
			library(qvalue)
			qobj <- qvalue(p = prim.pvalue)
			fdr = qobj$lfdr
			png(paste(pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else {
			stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
		}
		pass1.status = rep("No Change", times=length(fdr))
		if (trt.group == "continuous"){
			pass1.status[(gene.cor >= cor.cutoff) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Up",sep="")
			pass1.status[(gene.cor <= -cor.cutoff) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Down",sep="")
		} else{
			pass1.status[(prim.fc >= fc.cutoff) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Up",sep="")
			pass1.status[(prim.fc <= -fc.cutoff) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Down",sep="")
		}#end else

		print(paste("Primary Up-Regulated: ",length(pass1.status[pass1.status == paste(trt.group," Up",sep="")]),sep=""))
		print(paste("Primary Down-Regulated: ",length(pass1.status[pass1.status == paste(trt.group," Down",sep="")]),sep=""))

		if (fdr.method == "BH"){
			sec.fdr = p.adjust(sec.pvalue, "fdr")
		} else if (fdr.method == "q-value"){
			library(qvalue)
			qobj <- qvalue(p = sec.pvalue)
			sec.fdr = qobj$qvalue
			png(paste(pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else if (fdr.method == "q-lfdr"){
			library(qvalue)
			qobj <- qvalue(p = sec.pvalue)
			sec.fdr = qobj$lfdr
			png(paste(pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else {
			stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
		}		

		pass2.status = rep("No Change", times=length(fdr))
		if (trt.group2 == "continuous"){
			pass2.status[(gene.cor2 >= cor.cutoff2) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Up",sep="")
			pass2.status[(gene.cor2 <= -cor.cutoff2) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Down",sep="")
		} else{
			pass2.status[(sec.fc >= fc.cutoff2) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Up",sep="")
			pass2.status[(sec.fc <= -fc.cutoff2) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Down",sep="")
		}#end else

		print(paste("Secondary Up-Regulated: ",length(pass2.status[pass2.status == paste(trt.group," Up",sep="")]),sep=""))
		print(paste("Secondary Down-Regulated: ",length(pass2.status[pass2.status == paste(trt.group," Down",sep="")]),sep=""))
			
		pvalue.table = data.frame(prim.pvalue = prim.pvalue, prim.FDR = fdr,
										sec.pvalue=sec.pvalue, sec.fdr=sec.fdr)
			
		status = rep("No Change", times=length(fdr))
		status[(pass1.status == paste(trt.group," Up",sep="")) & (pass2.status == "No Change")] = upID
		status[(pass1.status == paste(trt.group," Down",sep="")) & (pass2.status == "No Change")] = downID
	} else{
		stop("interaction must be \"no\", \"model\", or \"filter-overlap\"")
	}#end else
}#end else

print(paste("Up-Regulated: ",length(status[status == upID]),sep=""))
print(paste("Down-Regulated: ",length(status[status == downID]),sep=""))

if (interaction.flag == "filter-overlap"){
	pvalue.method = paste(pvalue.method,"two-step_filtered",sep="_")
}

if(rep.check == 1){
	deg.table = data.frame(transcript.cluster=transcript.cluster, symbol = genes,
							probes.per.cluster=probes.per.cluster, genes.per.probe=genes.per.probe,
							average.rma, fc.table,
							pvalue.table, status = status)
} else {
	deg.table = data.frame(transcript.cluster=transcript.cluster, symbol = genes,
							probes.per.cluster=probes.per.cluster, genes.per.probe=genes.per.probe,
							average.rma, fc.table, status = status)	
}#end else

deg.file = paste(comp.name,"_",pvalue.method,"_DEG_fc_",fc.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".txt",sep="")
deg.file = gsub(":",".",deg.file)
write.table(deg.table, file=deg.file, row.names=F, quote=F, sep="\t")

final.deg.file = paste(user.folder,"/DEG/",comp.name,"_DEG_stats.txt",sep="")
write.table(deg.table, file=final.deg.file, row.names=F, quote=F, sep="\t")

IPA.file = deg.table[(deg.table$genes.per.probe==1) & !is.na(deg.table$genes.per.probe),]
print(dim(IPA.file))
write.table(IPA.file, file=paste(user.folder,"/IPA/",comp.name,"_for_IPA.txt",sep=""),sep="\t", row.names=F, quote=F)

panther.background.genes = genes[(deg.table$genes.per.probe==1) & !is.na(deg.table$genes.per.probe)]
print(length(panther.background.genes))
print(length(unique(panther.background.genes)))
panther.background.file = paste(user.folder,"/GO/Input_Files/",comp.name,"_BACKGROUND_for_PANTHER.txt",sep="")
write.table(data.frame(genes=panther.background.genes),panther.background.file, quote=F, row.names=F)

panther.up.genes = genes[(deg.table$status == upID) & (deg.table$genes.per.probe==1) & !is.na(deg.table$genes.per.probe)]
print(length(panther.up.genes))
print(length(unique(panther.up.genes)))
panther.up.file = paste(user.folder,"/GO/Input_Files/",comp.name,"_UP_for_PANTHER.txt",sep="")
write.table(data.frame(genes=panther.up.genes),panther.up.file, quote=F, row.names=F)

panther.down.genes = genes[(deg.table$status == downID) & (deg.table$genes.per.probe==1) & !is.na(deg.table$genes.per.probe)]
print(length(panther.down.genes))
print(length(unique(panther.down.genes)))
panther.down.file = paste(user.folder,"/GO/Input_Files/",comp.name,"_DOWN_for_PANTHER.txt",sep="")
write.table(data.frame(genes=panther.down.genes),panther.down.file, quote=F, row.names=F)

temp.rma = rma.mat
temp.rma = temp.rma[status != "No Change", ]
deg.genes = genes[status != "No Change"]

if(length(deg.genes) > 1){
	if(length(plot.groups) > 1){
		source("heatmap.3.R")
		if((trt.group != "continuous")&(trt.group2 != "continuous")){
			grp1 = as.character(sample.description.table[,plot.groups[1]])
			grp2 = as.character(sample.description.table[,plot.groups[2]])

			group.levels = c(levels(as.factor(grp1)),levels(as.factor(grp2)))

			color.palette <- fixed.color.palatte[1:length(group.levels)]
			labelColors1 = rep("black",times=length(sample.label))
			for (i in 1:length(group.levels)){
				labelColors1[grp1 == as.character(group.levels[i])] = color.palette[i]
			}#end for (i in 1:length(group.levels))
			labelColors2 = rep("black",times=length(sample.label))
			for (i in 1:length(group.levels)){
				labelColors2[grp2 == as.character(group.levels[i])] = color.palette[i]
			}#end for (i in 1:length(group.levels))
		}else if((trt.group == "continuous")&(trt.group2 == "continuous")){
			stop("Add code for two continuous variables")
		}else if(trt.group == "continuous"){
			grp1 = as.numeric(sample.description.table[,plot.groups[1]])
			grp2 = as.character(sample.description.table[,plot.groups[2]])
		
			labelColors1 = rep("black",times=length(sample.label))
			library("RColorBrewer")
			continuous.color.breaks = 10
				
			plot.var = as.numeric(grp1)
			plot.var.min = min(plot.var, na.rm=T)
			plot.var.max = max(plot.var, na.rm=T)
				
			plot.var.range = plot.var.max - plot.var.min
			plot.var.interval = plot.var.range / continuous.color.breaks
				
			color.range = colorRampPalette(c("green","black","orange"))(n = continuous.color.breaks)
			plot.var.breaks = plot.var.min + plot.var.interval*(0:continuous.color.breaks)
			for (j in 1:continuous.color.breaks){
				#print(paste(plot.var.breaks[j],"to",plot.var.breaks[j+1]))
				labelColors1[(plot.var >= plot.var.breaks[j]) &(plot.var <= plot.var.breaks[j+1])] = color.range[j]
			}#end for (j in 1:continuous.color.breaks)
			
			group.levels = c(levels(as.factor(grp2)))
			color.palette <- fixed.color.palatte[3:(2+length(group.levels))]
			labelColors2 = rep("black",times=length(sample.label))
			for (i in 1:length(group.levels)){
				labelColors2[grp2 == as.character(group.levels[i])] = color.palette[i]
			}#end for (i in 1:length(group.levels))
		}else{
			grp1 = as.character(sample.description.table[,plot.groups[1]])
			grp2 = as.numeric(sample.description.table[,plot.groups[2]])

			group.levels = c(levels(as.factor(grp1)))
			color.palette <- fixed.color.palatte[1:(length(group.levels))]
			labelColors1 = rep("black",times=length(sample.label))
			for (i in 1:length(group.levels)){
				labelColors1[grp1 == as.character(group.levels[i])] = color.palette[i]
			}#end for (i in 1:length(group.levels))
			
			labelColors2 = rep("black",times=length(sample.label))
			library("RColorBrewer")
			continuous.color.breaks = 10
				
			plot.var = as.numeric(grp2)
			plot.var.min = min(plot.var, na.rm=T)
			plot.var.max = max(plot.var, na.rm=T)
				
			plot.var.range = plot.var.max - plot.var.min
			plot.var.interval = plot.var.range / continuous.color.breaks
				
			color.range = colorRampPalette(c("purple","black","cyan"))(n = continuous.color.breaks)
			plot.var.breaks = plot.var.min + plot.var.interval*(0:continuous.color.breaks)
			for (j in 1:continuous.color.breaks){
				#print(paste(plot.var.breaks[j],"to",plot.var.breaks[j+1]))
				labelColors2[(plot.var >= plot.var.breaks[j]) &(plot.var <= plot.var.breaks[j+1])] = color.range[j]
			}#end for (j in 1:continuous.color.breaks)
		}
		
		std.expr = apply(temp.rma, 1, standardize.arr)
		if(length(deg.genes) < 25){
			colnames(std.expr) = deg.genes
		} else {
			colnames(std.expr) = rep("", length(deg.genes))
		}
		rownames(std.expr) = sample.label

		column_annotation <- as.matrix(deg.genes)
		colnames(column_annotation) <- c("")

		row_annotation <- data.frame(label1 = labelColors1, label2 = labelColors2)
		row_annotation = as.matrix(t(row_annotation))
		rownames(row_annotation) <- c(plot.groups)

		heatmap.file <- paste(comp.name,"_",pvalue.method,"_DEG_fc_",fc.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".png",sep="")
		heatmap.file = gsub(":",".",heatmap.file)
		png(file = heatmap.file)
		heatmap.3(std.expr, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
					RowSideColors=row_annotation, trace="none", margins = c(10,15),RowSideColorsSize=4, dendrogram="both")
		if((trt.group != "continuous")&(trt.group2 != "continuous")){
					legend("topright", legend=group.levels,
							col=color.palette,
							pch=15, cex=0.7)
		}else if((trt.group == "continuous")&(trt.group2 == "continuous")){
			stop("Add code for two continuous variables")
		}else if(trt.group == "continuous"){
			legend("right",legend=c(round(plot.var.max,digits=1),rep("",length(color.range)-2),round(plot.var.min,digits=1)),
								col=rev(color.range),  pch=15, y.intersp = 0.4, cex=0.8, pt.cex=1.5)
			legend("topright", legend=group.levels, col=color.palette, pch=15, cex=0.7)
		}else{
			legend("right",legend=c(round(plot.var.max,digits=1),rep("",length(color.range)-2),round(plot.var.min,digits=1)),
								col=rev(color.range),  pch=15, y.intersp = 0.4, cex=0.8, pt.cex=1.5)
			legend("topright", legend=group.levels, col=color.palette, pch=15, cex=0.7)
		}
		dev.off()
			
		if(interaction.flag != "no"){
			temp.fc.table = as.matrix(fc.table)
			if (((trt.group == "continuous") & (trt.group2 == "continuous")) | ((trt.group != "continuous") & (trt.group2 != "continuous"))){
				temp.fc.table = temp.fc.table[,-ncol(temp.fc.table)]
			}
			temp.fc.table = temp.fc.table[status != "No Change", ]
			if(length(deg.genes) < 25){
				rownames(temp.fc.table) = deg.genes
			} else {
				rownames(temp.fc.table) = rep("",times=length(deg.genes))
			}
			colnames(temp.fc.table) = gsub(".:.",":",gsub("fold.change.","",colnames(temp.fc.table)))
		
			temp.fc.table[temp.fc.table < -10] = -10
			temp.fc.table[temp.fc.table > 10] = 10
		
			heatmap.file <- paste("fold_change_",comp.name,"_",pvalue.method,"_DEG_fc_",fc.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".png",sep="")
			heatmap.file = gsub(":",".",heatmap.file)
			png(file = heatmap.file)
			heatmap.2(temp.fc.table, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
						trace="none", margins = c(20,5), cexCol=1.5)
			dev.off()
		}#end if(interaction.flag != "no")
		
	} else {
		labelColors = rep("black",times=length(sample.label))
		if(trt.group == "continuous"){
			library("RColorBrewer")
			continuous.color.breaks = 10
			
			plot.var = as.numeric(sample.description.table[,plot.groups])
			plot.var.min = min(plot.var, na.rm=T)
			plot.var.max = max(plot.var, na.rm=T)
			
			plot.var.range = plot.var.max - plot.var.min
			plot.var.interval = plot.var.range / continuous.color.breaks
			
			color.range = colorRampPalette(c("green","black","orange"))(n = continuous.color.breaks)
			plot.var.breaks = plot.var.min + plot.var.interval*(0:continuous.color.breaks)
			for (j in 1:continuous.color.breaks){
				#print(paste(plot.var.breaks[j],"to",plot.var.breaks[j+1]))
				labelColors[(plot.var >= plot.var.breaks[j]) &(plot.var <= plot.var.breaks[j+1])] = color.range[j]
			}#end for (j in 1:continuous.color.breaks)
		}else{
			group.levels = levels(as.factor(sample.description.table[,plot.groups]))
			color.palette = fixed.color.palatte[1:length(group.levels)]
			for (i in 1:length(group.levels)){
				labelColors[grp == as.character(group.levels[i])] = color.palette[i]
			}#end for (i in 1:length(group.levels))
		}

		std.expr = apply(temp.rma, 1, standardize.arr)
		if(length(deg.genes) < 25){
			colnames(std.expr) = deg.genes
		} else {
			colnames(std.expr) = rep("", length(deg.genes))
		}
		rownames(std.expr) = sample.label
		
		heatmap.file <- paste(comp.name,"_DEG_fc_",fc.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".png",sep="")
		heatmap.file = gsub(":",".",heatmap.file)
		png(file = heatmap.file)
		heatmap.2(std.expr, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
					 RowSideColors=labelColors, trace="none", margins = c(5,15))

		if(trt.group == "continuous"){
			legend("right",legend=c(round(plot.var.max,digits=1),rep("",length(color.range)-2),round(plot.var.min,digits=1)),
								col=rev(color.range),  pch=15, y.intersp = 0.4, cex=0.8, pt.cex=1.5)
		}else{
			legend("topright", group.levels, col=color.palette, pch=15)
		}
		dev.off()
	}#end else
}#end if(length(deg.genes) > 1)

if((length(deg.genes) > 1) & (run.gostat == "yes")){
	if(gostat.type == "GO.db"){
		print("Collecting GO annotations")
		library(GO.db)
		library(org.Hs.eg.db)
		orgdb = org.Hs.eg.db
		orgdb.GOIDs = keys(orgdb, keytype="GOALL")
		hg.goid.table = select(orgdb, keys=orgdb.GOIDs, columns=c("SYMBOL"), keytype="GOALL")
		
		GO.ID = keys(GO.db)
		GO.term = c()
		GO.type = c()
		GO.genes = c()
		GO.total.gene.count = c()
		
		for (i in 1:length(GO.ID)){
			GO.term[i]=as.character(Term(GO.ID[i]))
			GO.type[i]=as.character(Ontology(GO.ID[i]))
			
			if (GO.ID[i] %in% orgdb.GOIDs){
				temp.symbols = unique(as.character(hg.goid.table$SYMBOL[hg.goid.table$GO == GO.ID[i]]))
				GO.genes[i]=paste(temp.symbols, collapse=",")
				GO.total.gene.count[i]=length(temp.symbols)
			}else{
				GO.genes[i]=NA
				GO.total.gene.count[i]=NA
			}
		}#end for (i in 1:length(GO.ID))
		
		GO.info = data.frame(GO.ID=GO.ID, GO.term=GO.term, GO.type=GO.type, total.gene.count = GO.total.gene.count)
		print(dim(GO.info))
		GO.info = GO.info[!is.na(GO.info$total.gene.count) & (GO.info$total.gene.count > 20),]
		print(dim(GO.info))
		
		background.genes = unique(as.character(panther.background.genes))
		
		print("Performing GO enrichment for up-regulated genes")
		up.genes = unique(as.character(panther.up.genes))
		if (length(up.genes) > 1){
			up.pval = c()
			up.num.genes = c()
			up.matched.genes = c()
			
			for (i in 1:nrow(GO.info)){
				full.gene.set = unlist(strsplit(GO.genes[GO.ID == GO.info$GO.ID[i]],split=","))
				matched.up = up.genes[match(full.gene.set, up.genes, nomatch=0)]
				up.num.genes[i]=length(matched.up)
				if (up.num.genes[i] == 0){
					up.matched.genes[i]=NA
					up.pval[i]=NA
				}else{
					up.matched.genes[i]=paste(matched.up,collapse=",")

					mat = matrix(c(length(matched.up), length(up.genes)-length(matched.up),
						       length(full.gene.set),length(background.genes)-length(full.gene.set)),ncol=2)
					result = fisher.test(mat, alternative="greater")
					up.pval[i]=result$p.value
				}#end else
			}#end for (i in 1:nrow(GO.info))
			
			up.fdr=p.adjust(up.pval,"fdr")
			
			up.go.table = data.frame(GO.info, num.gene.list = up.num.genes,
									over.enrichment.pvalue=up.pval, over.enrichment.fdr=up.fdr,
									genes=up.matched.genes)
			up.go.table = up.go.table[order(up.go.table$over.enrichment.pvalue),]
			up.go.file = paste(user.folder,"/GO/",comp.name,"_FE_GO_enrichment_UP.txt",sep="")
			write.table(up.go.table, up.go.file, row.names=F, sep="\t")
		}#end if (length(up.genes) > 1)
		
		print("Performing GO enrichment for down-regulated genes")
		down.genes = unique(as.character(panther.down.genes))
		if (length(down.genes) > 1){
			down.pval = c()
			down.num.genes = c()
			down.matched.genes = c()
			
			for (i in 1:nrow(GO.info)){
				full.gene.set = unlist(strsplit(GO.genes[GO.ID == GO.info$GO.ID[i]],split=","))
				matched.down = down.genes[match(full.gene.set, down.genes, nomatch=0)]
				down.num.genes[i]=length(matched.down)
				if (down.num.genes[i] == 0){
					down.matched.genes[i]=NA
					down.pval[i]=NA
				}else{
					down.matched.genes[i]=paste(matched.down,collapse=",")

					mat = matrix(c(length(matched.down), length(down.genes)-length(matched.down),
						       length(full.gene.set),length(background.genes)-length(full.gene.set)),ncol=2)
					result = fisher.test(mat, alternative="greater")
					down.pval[i]=result$p.value
				}#end else
			}#end for (i in 1:nrow(GO.info))
			
			down.fdr=p.adjust(down.pval,"fdr")
			
			down.go.table = data.frame(GO.info, num.gene.list = down.num.genes,
									over.enrichment.pvalue=down.pval, over.enrichment.fdr=down.fdr,
									genes=down.matched.genes)
			down.go.table = down.go.table[order(down.go.table$over.enrichment.pvalue),]
			down.go.file = paste(user.folder,"/GO/",comp.name,"_FE_GO_enrichment_DOWN.txt",sep="")
			write.table(down.go.table, down.go.file, row.names=F, sep="\t")
		}#end if (length(down.genes) > 1)
	}else if(gostat.type == "GOstats"){
		library(GOstats)
		print("Running GOstats")
		
		library(org.Hs.eg.db)
		orgdb = org.Hs.eg.db
		orgdb.symbols = keys(orgdb, keytype="SYMBOL")
		geneID.table = select(orgdb, keys=orgdb.symbols, columns=c("ENTREZID"), keytype="SYMBOL")
		
		entrezUniverse =unique(geneID.table$ENTREZID)
		
		print("Performing GO enrichment for up-regulated genes")
		up.genes = panther.up.genes
		if (length(up.genes) > 1){
			selectedEntrezIds = geneID.table$ENTREZID[match(up.genes, geneID.table$SYMBOL,nomatch=0)]
			
			print("BP test")
			params = new("GOHyperGParams",
							geneIds=selectedEntrezIds,
							universeGeneIds=entrezUniverse,
							annotation="org.Hs.eg",
							ontology="BP",
							pvalueCutoff=1,
							conditional=TRUE,
							testDirection="over")
			bpOver = hyperGTest(params)
			bp.mat = summary(bpOver)
			bp.mat = data.frame(bp.mat,term.type = rep("BP",nrow(bp.mat)))
			#htmlReport(bpOver, file="test.html")
			names(bp.mat)=c("GOID","Pvalue","OddsRatio","ExpCount","Count","Size","Term","term.type")
			
			print("MF test")
			params = new("GOHyperGParams",
							geneIds=selectedEntrezIds,
							universeGeneIds=entrezUniverse,
							annotation="org.Hs.eg",
							ontology="MF",
							pvalueCutoff=1,
							conditional=TRUE,
							testDirection="over")
			mfOver = hyperGTest(params)
			mf.mat = summary(mfOver)
			mf.mat = data.frame(mf.mat,term.type = rep("MF",nrow(mf.mat)))
			names(mf.mat)=c("GOID","Pvalue","OddsRatio","ExpCount","Count","Size","Term","term.type")
								
			print("CC test")
			params = new("GOHyperGParams",
							geneIds=selectedEntrezIds,
							universeGeneIds=entrezUniverse,
							annotation="org.Hs.eg",
							ontology="CC",
							pvalueCutoff=1,
							conditional=TRUE,
							testDirection="over")
			ccOver = hyperGTest(params)
			cc.mat = summary(ccOver)
			cc.mat = data.frame(cc.mat,term.type = rep("CC",nrow(cc.mat)))
			names(cc.mat)=c("GOID","Pvalue","OddsRatio","ExpCount","Count","Size","Term","term.type")

			
			up.go.table = rbind(bp.mat, mf.mat, cc.mat)
			up.go.table = data.frame(up.go.table[,1:2], FDR=p.adjust(up.go.table$Pvalue,"fdr"), up.go.table[,3:8])
			up.go.table = up.go.table[order(up.go.table$Pvalue),]
			up.go.file = paste(user.folder,"/GO/",comp.name,"_GOstat_enrichment_UP.txt",sep="")
			write.table(up.go.table, up.go.file, row.names=F, sep="\t")
		}#end if (length(up.genes) > 1)
		
		print("Performing GO enrichment for down-regulated genes")
		down.genes = panther.down.genes
		if (length(down.genes) > 1){
			selectedEntrezIds = geneID.table$ENTREZID[match(down.genes, geneID.table$SYMBOL,nomatch=0)]
			
			print("BP test")
			params = new("GOHyperGParams",
							geneIds=selectedEntrezIds,
							universeGeneIds=entrezUniverse,
							annotation="org.Hs.eg",
							ontology="BP",
							pvalueCutoff=1,
							conditional=TRUE,
							testDirection="over")
			bpOver = hyperGTest(params)
			bp.mat = summary(bpOver)
			bp.mat = data.frame(bp.mat,term.type = rep("BP",nrow(bp.mat)))
			#htmlReport(bpOver, file="test.html")
			names(bp.mat)=c("GOID","Pvalue","OddsRatio","ExpCount","Count","Size","Term","term.type")
			
			print("MF test")
			params = new("GOHyperGParams",
							geneIds=selectedEntrezIds,
							universeGeneIds=entrezUniverse,
							annotation="org.Hs.eg",
							ontology="MF",
							pvalueCutoff=1,
							conditional=TRUE,
							testDirection="over")
			mfOver = hyperGTest(params)
			mf.mat = summary(mfOver)
			mf.mat = data.frame(mf.mat,term.type = rep("MF",nrow(mf.mat)))
			names(mf.mat)=c("GOID","Pvalue","OddsRatio","ExpCount","Count","Size","Term","term.type")
								
			print("CC test")
			params = new("GOHyperGParams",
							geneIds=selectedEntrezIds,
							universeGeneIds=entrezUniverse,
							annotation="org.Hs.eg",
							ontology="CC",
							pvalueCutoff=1,
							conditional=TRUE,
							testDirection="over")
			ccOver = hyperGTest(params)
			cc.mat = summary(ccOver)
			cc.mat = data.frame(cc.mat,term.type = rep("CC",nrow(cc.mat)))
			names(cc.mat)=c("GOID","Pvalue","OddsRatio","ExpCount","Count","Size","Term","term.type")

			
			down.go.table = rbind(bp.mat, mf.mat, cc.mat)
			down.go.table = data.frame(down.go.table[,1:2], FDR=p.adjust(down.go.table$Pvalue,"fdr"), down.go.table[,3:8])
			down.go.table = down.go.table[order(down.go.table$Pvalue),]
			down.go.file = paste(user.folder,"/GO/",comp.name,"_GOstat_enrichment_DOWN.txt",sep="")
			write.table(down.go.table, down.go.file, row.names=F, sep="\t")
		}#end if (length(down.genes) > 1)
		
	}else{
		stop("To run GO enrichment in R, 'R_GO_type' must be 'GOstats' or 'GO.db'")
	}
}#end if((length(deg.genes) > 1) & (run.gostat == "yes"))
