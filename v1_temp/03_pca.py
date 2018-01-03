import hail
import sys
import timeit

start = timeit.default_timer()

chrom = str(sys.argv[1])

hc = hail.HailContext(log='/hail.log', default_reference='GRCh38', min_block_size=4096)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# input
splitvds_file = 'gs://wgspd-wgs-v2/raw/vds/' + chrom + '/wgspd_wgs_v2_split.vds'
sample_qc_file = 'gs://wgspd-wgs-v2/qc_measures/' + chrom + '/samples_to_keep_after_qc.txt'
variant_ldprune_file = 'gs://wgspd-wgs-v2/raw/' + chrom + '/wgspd_wgs_v2_idependent_SNPs.prune.in'
sample_sex_fstat_file = 'gs://wgspd-wgs-v2/qc_measures/' + chrom + '/sex_related_annot.txt'
# check if this file is different from other sex sample file
mhc_chr8inv_file = 'gs://hail-common/MHC_invchr8_longLDreg.txt'


# output
pca_output_file = 'gs://wgspd-wgs-v2/qc_measures/' + chrom + '/pca.tsv'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Save a high-quality variants dataset in plink and perform sex-check ##
vds = hc.read(splitvds_file)


(wgspd_wgs_v2_raw_split
#.rename_samples('gs://wgspd-wgs-v1/qc_measures/sample_ids_all_spaces_to_underscores.txt')
.filter_samples_list(sample_qc_file, keep=True)
.annotate_samples_expr('sa = drop(sa, ID)')
.annotate_samples_table(sample_sex_fstat_file,'ID',impute = True,code = 'sa=merge(sa, table)')
.filter_variants_list(variant_ldprune_file, keep=True)
.filter_variants_intervals(mhc_chr8inv_file, keep=False)
.filter_variants_expr('v.contig == "X" ', keep=False)
.filter_samples_expr('sa.first_2nd_related==1', keep=False)
.pca('sa.pca')
.export_samples(pca_output_file,'id=s.id, sa.pca.*')



# print runtime

stop = timeit.default_timer()
print "runtime: " + str(stop - start) + " seconds"

