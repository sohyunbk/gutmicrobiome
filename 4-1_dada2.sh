mkdir 4.QualityControl
qiime dada2 denoise-single --p-n-threads 28 --i-demultiplexed-seqs ./3.Importing/demultiplexed.qza --p-trunc-len 0  --p-trim-left 0 --o-representative-sequences ./4.QualityControl/rep_seqs.qza --o-table ./4.QualityControl/table.qza --o-denoising-stats ./3.Importing/stats-dada2.qza
