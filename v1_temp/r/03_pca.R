library(tidyverse)

setwd('~/Dropbox/postdoc/projects/WGSPD/')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get files from bucket
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gsutil = '/Users/robert/google-cloud-sdk/bin/gsutil'
system('gsutil cp gs://wgspd-wgs-v1_temp/qc_measures/chr20/pca_values.tsv ~/Dropbox/postdoc/projects/WGSPD/out/')
system('gsutil cp gs://wgspd-wgs-v1_temp/qc_measures/chr20/pca_scores.tsv ~/Dropbox/postdoc/projects/WGSPD/out/')
system('gsutil cp gs://wgspd-wgs-v1_temp/pheno/final_pheno_2018_01_04.tsv ~/Dropbox/postdoc/projects/WGSPD/out/')
system('gsutil cp gs://wgspd-wgs-v1_temp/qc_measures/chr20/sample_qc_info.txt ~/Dropbox/postdoc/projects/WGSPD/out/')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ibd_exc = read.table('out/ibd_greater_177.txt', stringsAsFactors = F)[,1]
relpairs = read.table('out/ibd_greater_177_pairs.txt', stringsAsFactors = F, h=T)

values = read.table('out/pca_values.tsv', stringsAsFactors = F)
scores = read.table('out/pca_scores.tsv', stringsAsFactors = F, h=T)
pops = read.table('out/sample_qc_info.txt', stringsAsFactors = F, h=T, comment.char = '', sep='\t') %>% select('ID', 'FINAL_POP')
pheno = read.table('out/final_pheno_2018_01_04.tsv', stringsAsFactors = F, h=T, comment.char = '', sep='\t')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot pcs against populations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m = left_join(scores, select(pheno, s='ID', 'FINAL_POP'), by='s')

ggplot(m, aes(pcaScores.PC1, pcaScores.PC2, col=FINAL_POP)) + geom_point()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# assign PC scores to relatives
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

imputed = left_join(select(relpairs, 'i', s='j'), scores) %>% ddply('i', function(x) apply(x[,-1:-2], 2, function(y) mean(y, na.rm=T)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add PCs to phenotype file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pheno2 = left_join(pheno, bind_rows(dplyr::rename(scores, ID='s'), dplyr::rename(imputed, ID='i')))

write.table('out/final_pheno_2018_01_04_PCs.tsv', col.names=T, row.names=F, quote=F, sep='\t')
system(paste0(gsutil, ' cp ~/Dropbox/postdoc/projects/WGSPD/out/final_pheno_2018_01_04_PCs.tsv gs://wgspd-wgs-v1_temp/pheno/'))





