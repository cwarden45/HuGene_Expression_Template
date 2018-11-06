import sys
import re
import os

aptRoot = "/Volumes/user_data/Seq/cwarden/Array_Annotation_Files/apt-1.19.0-x86_64-apple-yosemite/bin/"
cdfFile = "/Volumes/user_data/Seq/cwarden/Array_Annotation_Files/CD_HG-U133_Plus_2/Full/HG-U133_Plus_2/LibFiles/HG-U133_Plus_2.cdf"
mps = "probeset_gene.mps"

command = aptRoot + "apt-probeset-summarize -a rma-sketch -d " + cdfFile + " -o RMA --cel-files cel_list.txt"
#os.system(command)

command = aptRoot + "apt-probeset-summarize -a rma-sketch -d " + cdfFile + " -m " + mps+ " -o custom_RMA --cel-files cel_list.txt"
os.system(command)

