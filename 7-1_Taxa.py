import os, sys
os.system("mkdir 7.FeatureTable") 
for nLevel in range(1,8):
	print(nLevel)
	sLevel = str(nLevel)
	cmd = "qiime taxa collapse \
	--i-table ./4.QualityControl/table.qza \
	--i-taxonomy ./6.taxonomy/taxonomy.qza \
	--p-level %s \
	--o-collapsed-table ./7.FeatureTable/FeatureTable_L%s.qza"%(sLevel, sLevel)
	print(cmd)
	os.system(cmd)
