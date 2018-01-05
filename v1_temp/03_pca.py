from hail import *
import sys
import timeit

start = timeit.default_timer()

chrom = str(sys.argv[1])

hc = HailContext(log='/hail.log', default_reference='GRCh37', min_block_size=4096)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# input
vds_ldpruned_common_file = 'gs://wgspd-wgs-v1_temp/qced/vds/' + chrom + '/ldpruned_common.vds'
#sample_qc_file = 'gs://wgspd-wgs-v1_temp/qc_measures/' + chrom + '/samples_to_keep_after_qc.txt'
sample_sex_fstat_file = 'gs://wgspd-wgs-v1_temp/qc_measures/' + chrom + '/sample_sex_fstat.txt'
# check if this file is different from other sex sample file
mhc_chr8inv_file = 'gs://hail-common/MHC_invchr8_longLDreg.txt'
rel_exclusion_file = 'gs://wgspd-wgs-v1_temp/qc_measures/' + chrom + '/ibd_greater_177.txt'



# output
pca_value_file = 'gs://wgspd-wgs-v1_temp/qc_measures/' + chrom + '/pca_values.tsv'
pca_score_file = 'gs://wgspd-wgs-v1_temp/qc_measures/' + chrom + '/pca_scores.tsv'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data and filter samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## interval list
mhc_chr8inv = KeyTable.import_interval_list(mhc_chr8inv_file)
##
rel_exclusion = hc.import_table(rel_exclusion_file, no_header=True).key_by('f0')

vds = hc.read(vds_ldpruned_common_file).filter_samples_expr('sa.meta.dropped_HetStat||sa.meta.dropped_sexStat||sa.meta.dropped_chimera||sa.meta.dropped_contam||sa.meta.dropped_DP||sa.meta.dropped_IS||sa.meta.dropped_hardfilters', keep=False)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pca
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

eigenvalues, scores, _  = (vds
#.rename_samples('gs://wgspd-wgs-v1/qc_measures/sample_ids_all_spaces_to_underscores.txt')
.filter_samples_table(rel_exclusion, keep=False)
#.annotate_samples_expr('sa = drop(sa, ID)')
#.annotate_samples_table(sample_sex_fstat_file,'ID',impute = True,code = 'sa=merge(sa, table)')
.filter_variants_table(mhc_chr8inv, keep=False)
.filter_variants_expr('v.contig == "X" || v.contig == "Y"', keep=False)
.filter_intervals(Interval.parse('6:25M-35M'), keep=False)
.filter_intervals(Interval.parse('8:7M-15M'), keep=False)
#.filter_samples_expr('sa.first_2nd_related==1', keep=False)
.pca('if (isDefined(g.GT)) 1 else 0', k=10)
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write output
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

with hadoop_write(pca_value_file) as f:
	for val in eigenvalues:
		print>>f, val
scores.flatten().export(pca_score_file)

# print runtime

stop = timeit.default_timer()
print "runtime: " + str(stop - start) + " seconds"

