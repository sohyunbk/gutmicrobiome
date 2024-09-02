qiime feature-table summarize --i-table ./4.QualityControl/table.qza --o-visualization ./4.QualityControl/table.qzv --m-sample-metadata-file ./0.MetaInfo/sample_metadata.tsv
qiime feature-table tabulate-seqs --i-data ./4.QualityControl/rep_seqs.qza --o-visualization ./4.QualityControl/rep_seqs.qzv

