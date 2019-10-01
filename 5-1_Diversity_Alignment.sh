mkdir 5.Diversity
qiime alignment mafft \
--i-sequences ./4.QualityControl/rep_seqs.qza \
--o-alignment ./5.Diversity/aligned-rep-seqs.qza
