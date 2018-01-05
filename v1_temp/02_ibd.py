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
splitvds_file = 'gs://wgspd-wgs-v1_temp/raw/vds/' + chrom + '/wgspd_wgs_v1_split.vds'
idmap_file = 'gs://wgspd-wgs-v1_temp/qc_measures/sample_ids_all_spaces_to_underscores.txt'
sample_qc_info_file = 'gs://wgspd-wgs-v1_temp/qc_measures/' + chrom + '/sample_qc_info.txt'
rawvds_file = 'gs://wgspd-wgs-v1_temp/raw/vds/' + chrom + '/wgspd_wgs_v1.vds'
rawvds_file = 'gs://wgspd-wgs-v1_temp/raw/vds/wgspd_wgs_v1.vds/'

# output
sample_sex_fstat_file = 'gs://wgspd-wgs-v1_temp/qc_measures/' + chrom + '/sample_sex_fstat.txt'
ibd_results_file = 'gs://wgspd-wgs-v1_temp/qc_measures/' + chrom + '/ibd_results.txt'
vds_ldpruned_common_file = 'gs://wgspd-wgs-v1_temp/qced/vds/' + chrom + '/ldpruned_common.vds'

mapping_table = hc.import_table(idmap_file)
mapping_dict = {row.ID: row.ID2 for row in mapping_table.collect()}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data and filter samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Save a high-quality variants dataset in plink and perform sex-check ##
wgspd_wgs_v1_temp_raw_split = hc.read(splitvds_file).filter_samples_expr('sa.meta.dropped_HetStat||sa.meta.dropped_sexStat||sa.meta.dropped_chimera||sa.meta.dropped_contam||sa.meta.dropped_DP||sa.meta.dropped_IS||sa.meta.dropped_hardfilters', keep=False)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# variant QC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print "variant QC..."
vds = (wgspd_wgs_v1_temp_raw_split
.rename_samples(mapping_dict)
#.filter_variants_intervals('gs://hail-common/LCR.interval_list', keep = False)
.variant_qc()
.filter_genotypes('(g.AD[0] + g.AD[1]) / g.DP < 0.9 || (g.GT.isHomRef() && (g.AD[0] / g.DP < 0.9 || g.GQ < 20)) || (g.GT.isHet() && (g.AD[1] / g.DP < 0.20 || g.PL[0] < 20)) || (g.GT.isHomVar() && (g.AD[1] / g.DP < 0.9 || g.PL[0] < 20)) || g.DP > 200', keep=False)
.variant_qc()
.filter_variants_expr('va.info.QD > 4 && va.qc.callRate > 0.99 && va.qc.dpMean > 20 && va.qc.AF > 0.05 && va.filters.isEmpty() && va.qc.AF < 0.95', keep=True)
)

print "variant count after QC:"
print vds.count_variants()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sex imputation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print "sex imputation..."
(vds.impute_sex(maf_threshold=0.05)
.samples_table()
.select(['ID=s','sexFstat=sa.imputesex.Fstat','isFemale=sa.imputesex.isFemale'])
.export(sample_sex_fstat_file))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ld pruning
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print "LD pruning..."
vds = vds.ld_prune(40)

print "writing VDS..."
vds.write(vds_ldpruned_common_file, overwrite=True)
# calculate ibd

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IBD analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#vds = hc.read(vds_ldpruned_common_file)
vds.ibd(min=0.1).flatten().rename({'ibd.Z0': 'Z0', 'ibd.Z1': 'Z1', 'ibd.Z2': 'Z2', 'ibd.PI_HAT': 'PI_HAT'}).export(ibd_results_file)


# print runtime
stop = timeit.default_timer()
print "runtime: " + str(stop - start) + " seconds"


