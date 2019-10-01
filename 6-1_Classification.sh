qiime feature-classifier classify-sklearn \
  --i-classifier ./6.taxonomy/gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads ./4.QualityControl/rep_seqs.qza \
  --o-classification ./6.taxonomy/taxonomy.qza

qiime metadata tabulate \
  --m-input-file ./6.taxonomy/taxonomy.qza \
  --o-visualization ./6.taxonomy/taxonomy.qzv

