import sys
import re
import os

inputFolders = ["../../CEL_Files"]

parameterFile = "parameters.txt"
celList = ""

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

outHandle = open(celList, 'w')
text = "cel_files\n"
outHandle.write(text)

for folder in inputFolders:
	fileResults = os.listdir(folder)
	for testfile in fileResults:
		celResult = re.search(".CEL",testfile)
		if celResult:
			fullFile = folder + "/" + testfile
			text = fullFile + "\n"
			outHandle.write(text)