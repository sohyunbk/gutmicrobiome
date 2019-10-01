qiime diversity alpha-rarefaction \
--i-table ./4.QualityControl/table.qza \
--i-phylogeny ./5.Diversity/rooted-tree.qza \
--p-max-depth 10000 \
--m-metadata-file ./0.MetaInfo/sample_metadata.tsv \
--o-visualization ./5.Diversity/alpha-rarefaction.qzv
