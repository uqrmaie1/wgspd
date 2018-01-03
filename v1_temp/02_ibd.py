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

# output
sample_sex_fstat_file = 'gs://wgspd-wgs-v1_temp/qc_measures/' + chrom + '/sample_sex_fstat.txt'
plink_common_file = 'gs://wgspd-wgs-v1_temp/raw/' + chrom + '/wgspd_wgs_v1_temp_HQ_common'
ibd_results_file = 'gs://wgspd-wgs-v1_temp/qc_measures/' + chrom + '/ibd_results.txt'


mapping_table = hc.import_table(idmap_file)
mapping_dict = {row.ID: row.ID2 for row in mapping_table.collect()}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Save a high-quality variants dataset in plink and perform sex-check ##
wgspd_wgs_v2_raw_split = hc.read(splitvds_file)

vds = (wgspd_wgs_v2_raw_split
.rename_samples(mapping_dict)
#.filter_variants_intervals('gs://hail-common/LCR.interval_list', keep = False)
.variant_qc()
.filter_genotypes('(g.AD[0] + g.AD[1]) / g.DP < 0.9 || (g.GT.isHomRef() && (g.AD[0] / g.DP < 0.9 || g.GQ < 20)) || (g.GT.isHet() && (g.AD[1] / g.DP < 0.20 || g.PL[0] < 20)) || (g.GT.isHomVar() && (g.AD[1] / g.DP < 0.9 || g.PL[0] < 20)) || g.DP > 200', keep=False)
.variant_qc()
.filter_variants_expr('va.info.QD > 4 && va.qc.callRate > 0.99 && va.qc.dpMean > 20 && va.qc.AF > 0.05 && va.filters.isEmpty() && va.qc.AF < 0.95', keep=True)
#.export_plink(plink_common_file)
.impute_sex(maf_threshold=0.05))

(vds
.samples_table()
.select(['ID=s','sexFstat=sa.imputesex.Fstat','isFemale=sa.imputesex.isFemale'])
.export(sample_sex_fstat_file))
#.export_samples(sample_sex_fstat_file,'ID=s.id,sexFstat=sa.imputesex.Fstat,isFemale=sa.imputesex.isFemale')


# calculate ibd
vds.ibd(min=0.1).flatten().rename({'ibd.Z0': 'Z0', 'ibd.Z1': 'Z1', 'ibd.Z2': 'Z2', 'ibd.PI_HAT': 'PI_HAT'}).export(ibd_results_file)


# print runtime

stop = timeit.default_timer()
print "runtime: " + str(stop - start) + " seconds"





## Create a set of independent SNPs in plink ##
## This is done on a VM, that's why the paths are different ##
# plink -bfile /mnt/data/aganna/bucket/raw/wgspd_wgs_v1_HQ_common \
# --snps-only \
# --indep-pairwise 100 50 0.2 \
# --out /mnt/data/aganna/bucket/raw/wgspd_wgs_v1_idependent_SNPs

# plink -bfile /mnt/data/aganna/bucket/raw/wgspd_wgs_v1_HQ_common \
# --extract /mnt/data/aganna/bucket/raw/wgspd_wgs_v1_idependent_SNPs.prune.in \
# --hwe 0.000001 \
# --make-bed --out /mnt/data/aganna/bucket/raw/wgspd_wgs_v1_HQ_common_pruned

# ## Calculate kinship matrix with King ##
# king -b /mnt/data/aganna/bucket/raw/wgspd_wgs_v1_HQ_common_pruned.bed --kinship --related --ibs --prefix /mnt/data/aganna/bucket/qc_measures/relatedness

# ## Save file for the IBD plot below ##
# awk '{if($16 > 0.04) print $0}' /mnt/data/aganna/bucket/qc_measures/relatedness.ibs > /mnt/data/aganna/bucket/qc_measures/related_for_plot

# ## Save file with > 2nd degree relatives for tagging and exclusions ##
# awk '{if($16 > 0.177) print $0}' /mnt/data/aganna/bucket/qc_measures/relatedness.ibs > /mnt/data/aganna/bucket/qc_measures/1_2_degree_rel