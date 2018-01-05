import hail
import sys
import timeit

start = timeit.default_timer()

chrom = str(sys.argv[1])

hc = hail.HailContext(log='/hail.log', default_reference='GRCh37', min_block_size=4096)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# input
rawvds_file = 'gs://wgspd-wgs-v1_temp/raw/vds/' + chrom + '/wgspd_wgs_v1.vds'
splitvds_file = 'gs://wgspd-wgs-v1_temp/raw/vds/' + chrom + '/wgspd_wgs_v1_split.vds'

# output
sample_qc_info_file = 'gs://wgspd-wgs-v1_temp/qc_measures/' + chrom + '/sample_qc_info.txt'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

wgspd_wgs_v1_temp_raw_split = hc.read(splitvds_file)


(wgspd_wgs_v1_temp_raw_split
#.rename_samples('gs://wgspd-wgs-v1/qc_measures/sample_ids_all_spaces_to_underscores.txt')
#.filter_variants_intervals('gs://hail-common/LCR.interval_list', keep = False)
.sample_qc()
.samples_table().select(['ID=s', \
	'CHIMERAS=sa.meta.PCT_CHIMERAS', \
	'PROJECT_OR_COHORT=sa.meta.project_or_cohort', \
	'RESEARCH_PROJECT=sa.meta.research_project', \
	'TITLE=sa.meta.Title', \
	'PCR_FREE=sa.meta.pcr_free', \
	'POPULATION=sa.meta.population', \
	'FINAL_POP=sa.meta.final_pop', \
	'QC_POP=sa.meta.qc_pop', \
	'SEX=sa.meta.sex', \
	'callRate=sa.qc.callRate', \
	'nCalled=sa.qc.nCalled', \
	'nSingleton=sa.qc.nSingleton', \
	'dpMean=sa.qc.dpMean', \
	'gqMean=sa.qc.gqMean', \
	'rTiTv=sa.qc.rTiTv', \
	'rHetHomVar=sa.qc.rHetHomVar', \
	'rInsertionDeletion=sa.qc.rInsertionDeletion', \
	'dropped_HetStat=sa.meta.dropped_HetStat', \
	'dropped_sexStat=sa.meta.dropped_sexStat', \
	'dropped_chimera=sa.meta.dropped_chimera', \
	'dropped_contam=sa.meta.dropped_contam', \
	'dropped_DP=sa.meta.dropped_DP', \
	'dropped_IS=sa.meta.dropped_IS', \
	'dropped_fam_exac=sa.meta.dropped_fam_exac', \
	'dropped_fam=sa.meta.dropped_fam', \
	'dropped_hardfilters=sa.meta.dropped_hardfilters', \
	'keep=sa.meta.keep' \
	]).export(sample_qc_info_file))




# print runtime

stop = timeit.default_timer()

print "runtime: " + str(stop - start) + " seconds"











