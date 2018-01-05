library(tidyverse)

setwd('~/Dropbox/postdoc/projects/WGSPD/')

gsutil = '/Users/robert/google-cloud-sdk/bin/gsutil'
system(paste0(gsutil, ' cp gs://wgspd-wgs-v1_temp/qc_measures/allchr/sample_sex_fstat.txt ~/Dropbox/postdoc/projects/WGSPD/out/'))
system(paste0(gsutil, ' cp gs://wgspd-wgs-v1_temp/qc_measures/chr20/ibd_results.txt ~/Downloads/'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ibd = read.table('~/Downloads/ibd_results.txt', stringsAsFactors = F, h=T)
sex = read.table('out/sample_sex_fstat.txt', stringsAsFactors = F, h=T)
pheno = read.table('out/final_pheno_2018_01_04.tsv', stringsAsFactors = F, h=T, comment.char = '', sep='\t')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot IBD pairs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ggplot(ibd[ibd$Z1 > .1,], aes(Z0, Z1)) + geom_point()


ggplot(sample_n(ibd, 1e4), aes(ibs0/(ibs0+ibs1+ibs2), PI_HAT)) + geom_point()

hist(ibd$PI_HAT, 100, col='black'); abline(v=.177)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define samples to exclude
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define individuals who have 1st or 2nd degree relatives and should be excluded
relpairs = ibd %>% dplyr::filter(PI_HAT > .177)
exc = unique(relpairs[,1])

write.table(relpairs, 'out/ibd_greater_177_pairs.txt', quote=F, row.names=F, col.names=T, sep='\t')
write.table(exc, 'out/ibd_greater_177.txt', quote=F, row.names=F, col.names=F)

system(paste0(gsutil, ' cp /Users/robert/Dropbox/postdoc/projects/WGSPD/out/ibd_greater_177.txt gs://wgspd-wgs-v1_temp/qc_measures/chr20/'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sex check
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m = left_join(sex, select(pheno, 'ID', 'SEX'))
table(m$isFemale, m$SEX)



