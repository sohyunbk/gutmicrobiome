##################
### 17/11/12 #####
### by Sohyun ####
##################
import glob, os, sys

os.system("mkdir 3.Importing")

CurrentDir = "/disk12/11.SH_MouseGut"
outfile = open("./3.Importing/Manifest_33.txt","w")
outfile.write("# single-end PHRED 33 fastq manifest file for forward reads\n")
outfile.write("sample-id,absolute-filepath,direction\n")
for sFiles in glob.glob("./2.Pear/*.assembled.fastq"):
	FileName = os.path.split(sFiles)[1].replace(".assembled.fastq","")
	outfile.write(FileName+","+CurrentDir+sFiles[1:]+",forward\n")

outfile.close()
	
