import sys
import re
import os

aptRoot = ""
pgfFile = ""
clfFile = ""
bgpFile = ""
gene_symbol_mpsFile = ""
transcript_cluster_mpsFile = ""
celList = ""
RMAoutput = ""
summarization_type = ""
DABGoutput = ""

parameterFile = "parameters.txt"
inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]
	
	if param == "CEL_input":
		celList = value

	if param == "APT_Root":
		aptRoot = value

	if param == "summary_method":
		summarization_type = value

	if param == "RMA_Cluster_Folder":
		RMAoutput = value

	if param == "DABG_Probeset_Folder":
		DABGoutput = value

	if param == "MPS_Gene_File":
		gene_symbol_mpsFile = value

	if param == "MPS_Transcript_Cluster_File":
		transcript_cluster_mpsFile = value

	if param == "PGF_File":
		pgfFile = value

	if param == "CLF_File":
		clfFile = value

	if param == "BGP_File":
		bgpFile = value
		
mpsFile = ""
if summarization_type == "transcript_cluster":
	mpsFile = transcript_cluster_mpsFile
elif summarization_type == "gene_symbol":
	mpsFile = gene_symbol_mpsFile
else:
	print "'summary_method' must be 'transcript_cluster' or 'gene_symbol'"
	sys.exit()
	
command = aptRoot + "/apt-probeset-summarize -a rma -p " + pgfFile + " -c " + clfFile + " -b " + bgpFile + " -m "+ mpsFile +" -o "+RMAoutput+" --cel-files " + celList
os.system(command)	

command = aptRoot + "/apt-probeset-summarize -a dabg -p " + pgfFile + " -c " + clfFile + " -b " + bgpFile + " -o "+DABGoutput+" --cel-files " + celList
os.system(command)	


