import os, sys, glob
Files = glob.glob("./7.FeatureTable/*.qza")

for sFile in Files:
	FileName = os.path.split(sFile)[1].split(".qza")[0]
	cmd  = "biom convert -i ./7.FeatureTable/%s/feature-table.biom -o ./7.FeatureTable/%s.tsv --to-tsv"%(FileName,FileName)
	print(cmd)
	os.system(cmd)
