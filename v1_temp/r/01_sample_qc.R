library(tidyverse)

setwd('~/Dropbox/postdoc/projects/WGSPD/')

#tjphen = read.table('data/wgspd-genomes-phenotypes.tsv', h=T, sep='\t', stringsAsFactors = F)
#ids = read.table('data/wgspd_ids_13500.txt', stringsAsFactors = F)[,1]


sampleqc = read.table('~/Downloads/sample_qc_info.txt', h=T, stringsAsFactors = F, sep='\t')

sampleqc2 = sampleqc %>% arrange(RESEARCH_PROJECT, FINAL_POP)

sampleqc %>% select(ID) %>% mutate(ID2=gsub(' ', '_', ID)) %>% write.table(file='data/sample_ids_all_spaces_to_underscores.txt', quote=F, col.names=T, row.names=F, sep='\t')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot chimeras, call rate, coverage
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



p1 = ggplot(sampleqc2, aes(x=1:nrow(sampleqc), y=CHIMERAS, colour=RESEARCH_PROJECT, shape=dropped_chimera)) + geom_point() + theme(legend.position = 'none') + geom_hline(yintercept=.05, linetype=2) + ylab('Chimera rate') + xlab('')

# , shape=keep
p2 = ggplot(sampleqc2, aes(x=1:nrow(sampleqc), y=callRate, colour=RESEARCH_PROJECT)) + geom_point() + theme(legend.position = 'none') + geom_hline(yintercept=.95, linetype=2) + ylab('Call rate') + xlab('')

p3 = ggplot(sampleqc2, aes(x=1:nrow(sampleqc), y=dpMean, colour=RESEARCH_PROJECT, shape=dropped_DP)) + geom_point() + theme(legend.position = 'none') + geom_hline(yintercept=15, linetype=2) + ylab('Coverage') + xlab('')

p0 = ggplot(sampleqc2, aes(x=1:nrow(sampleqc), y=dpMean, colour=RESEARCH_PROJECT, shape=dropped_DP)) + geom_point() + theme() + geom_hline(yintercept=15, linetype=2) + ylab('Coverage') + xlab('')

p4 = g_legend(p0)

grid.arrange(p1, p2, p3, p4)
ggsave('plots/01_sample_qc_chimera_callrate_coverage.pdf')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot Ti/Tv, Het/Hom, Ins/Del
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sampleqc3 = sampleqc2 %>%
  group_by(TITLE) %>%
  mutate(index=row_number(),
         median_rTiTv=median(rTiTv, na.rm=T),
         mad_rTiTv=median(abs(rTiTv-median_rTiTv), na.rm=T),
         median_rHetHomVar=median(rHetHomVar, na.rm=T),
         mad_rHetHomVar=median(abs(rHetHomVar-median_rHetHomVar), na.rm=T),
         median_rInsertionDeletion=median(rInsertionDeletion, na.rm=T),
         mad_rInsertionDeletion=median(abs(rInsertionDeletion-median_rInsertionDeletion), na.rm=T))

ggplot(sampleqc3, aes(x=index, y=rTiTv, colour=FINAL_POP)) + geom_point() + theme(legend.position = 'none') + geom_hline(aes(yintercept=median_rTiTv + 4*mad_rTiTv), linetype=2) + geom_hline(aes(yintercept=median_rTiTv - 4*mad_rTiTv), linetype=2) + ylab('Ti/Tv ratio') + xlab('') + facet_wrap(~ TITLE, scale='free_x')

ggsave('plots/01_sample_qc_titv.pdf')


ggplot(sampleqc3, aes(x=index, y=rHetHomVar, colour=FINAL_POP)) + geom_point() + theme(legend.position = 'none') + geom_hline(aes(yintercept=median_rHetHomVar + 4*mad_rHetHomVar), linetype=2) + geom_hline(aes(yintercept=median_rHetHomVar - 4*mad_rHetHomVar), linetype=2) + ylab('Het/Hom ratio') + xlab('') + facet_wrap(~ TITLE, scale='free_x')

ggsave('plots/01_sample_qc_hethomvar.pdf')


ggplot(sampleqc3, aes(x=index, y=rInsertionDeletion, colour=FINAL_POP)) + geom_point() + theme(legend.position = 'none') + geom_hline(aes(yintercept=median_rInsertionDeletion + 4*mad_rInsertionDeletion), linetype=2) + geom_hline(aes(yintercept=median_rInsertionDeletion - 4*mad_rInsertionDeletion), linetype=2) + ylab('Ins/Del ratio') + xlab('') + facet_wrap(~ TITLE, scale='free_x')

ggsave('plots/01_sample_qc_insdel.pdf')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot Ti/Tv, Het/Hom, Ins/Del
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sampleqc4 = sampleqc3 %>% mutate(exc_callrate = callRate < .95,
                                 exc_chimeras = CHIMERAS > .05,
                                 exc_coverge = dpMean < 15,
                                 exc_titv = rTiTv > median_rInsertionDeletion + 4*mad_rInsertionDeletion |
                                            rTiTv < median_rInsertionDeletion - 4*mad_rInsertionDeletion,
                                 exc_hethom = rHetHomVar > median_rHetHomVar + 4*mad_rHetHomVar |
                                              rHetHomVar < median_rHetHomVar - 4*mad_rHetHomVar,
                                 exc_insdel = rInsertionDeletion > median_rInsertionDeletion + 4*mad_rInsertionDeletion |
                                              rInsertionDeletion > median_rInsertionDeletion + 4*mad_rInsertionDeletion)

sampleqc4 = sampleqc4 %>% mutate(keep_robqc = keep & !(exc_callrate |
                                                       exc_chimeras |
                                                       exc_coverge |
                                                       exc_titv |
                                                       exc_hethom |
                                                       exc_insdel))



write.table(select(sampleqc4, ID) %>% filter(keep_robqc), '~/Dropbox/postdoc/projects/WGSPD/qc/01_sample_qc_keep.txt', quote=F, col.names=F, row.names=F)
write.table(select(sampleqc4, ID) %>% filter(!keep_robqc), '~/Dropbox/postdoc/projects/WGSPD/qc/01_sample_qc_exclude.txt', quote=F, col.names=F, row.names=F)


