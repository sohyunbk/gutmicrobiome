qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ./5.Diversity/rooted-tree.qza \
  --i-table ./4.QualityControl/table.qza \
  --p-sampling-depth 1109 \
  --m-metadata-file ./0.MetaInfo/sample_metadata.tsv \
  --output-dir ./5.Diversity/core-metrics-results
