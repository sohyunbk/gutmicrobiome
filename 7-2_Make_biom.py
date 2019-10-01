import os, sys, glob
Files = glob.glob("./7.FeatureTable/*.qza")

for sFile in Files:
	FileName = os.path.split(sFile)[1].split(".qza")[0]
	cmd  = "qiime tools export --input-path %s --output-path ./7.FeatureTable/%s"%(sFile,FileName)
	print(cmd)
	os.system(cmd)
