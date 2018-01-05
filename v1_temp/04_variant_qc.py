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
splitvds_file = 'gs://wgspd-wgs-v1_temp/raw/vds/' + chrom + '/wgspd_wgs_v1_split.vds'
lcr_file = 'gs://hail-common/LCR.interval_list'
sample_qc_file = 'gs://wgspd-wgs-v1_temp/qc_measures/' + chrom + '/samples_to_keep_after_qc.txt'
sample_sex_fstat_file = 'gs://wgspd-wgs-v1_temp/qc_measures/' + chrom + '/sample_sex_fstat.txt'

# output
qced_file = 'gs://wgspd-wgs-v1_temp/qced/vds/' + chrom + '/wgspd_wgs_v1_temp_split.vds'
variant_qc_stats_file = 'gs://wgspd-wgs-v1_temp/qc_measures/' + chrom + '/variant_qc_numbers.txt'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define constants
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mincallrate = '0.98'
hwep = '0.000001'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

vds = hc.read(splitvds_file)

## LCR
lcr = KeyTable.import_interval_list(lcr_file)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# variant qc
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# annotate samples with project IDs

# annotate samples with result from sex imputation

# filter samples
vds = vds.filter_samples_expr('sa.meta.dropped_HetStat||sa.meta.dropped_sexStat||sa.meta.dropped_chimera||sa.meta.dropped_contam||sa.meta.dropped_DP||sa.meta.dropped_IS||sa.meta.dropped_hardfilters', keep=False)
vds = vds.filter_variants_table(lcr, keep = False)


step0 = step1 = step2 = step31 = step32 = step4 = 0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# step0: AC > 0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print "starting step 0..."
step00 = vds.count_variants()
vds = vds.variant_qc().filter_variants_expr('va.qc.AC > 0')
#step0 = vds.count_variants()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# step1: PASS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print "starting step 1..."
vds = vds.filter_variants_expr('va.filters.isEmpty()')
#step1 = vds.count_variants()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# step2: QD > 4
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print "starting step 2..."
vds = vds.filter_variants_expr('va.info.QD > 4')
#step2 = vds.count_variants()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# step3.1: set missing:
#    DP > 400
#    Heterozygous AND [ (AD alt/DP) < 20% OR PL ref < 20 OR (AD ref + AD alt) / DP < 90% ]
#    Homozygous ref AND [ GQ < 20 OR (AD ref / DP) < 0.9 ]
#    Homozygous alt AND [ PL ref < 20 OR (AD alt / DP) < 0.9 ]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print "starting step 3.1..."
vds = (vds.variant_qc().filter_genotypes(
     """
     g.DP > 400 ||
     (g.GT.isHomRef() && (g.AD[0] / g.DP < 0.9 || g.GQ < 20)) ||
     (g.GT.isHomVar() && (g.AD[1] / g.DP < 0.9 || g.PL[0] < 20)) ||
     (g.GT.isHet() && ( (g.AD[0] + g.AD[1]) / g.DP < 0.9 || g.AD[1] / g.DP < 0.20 || g.PL[0] < 20 || (!sa.fam.isFemale && ( v.inXNonPar || v.inYNonPar )) )) || (v.inYNonPar && sa.fam.isFemale)
     """, keep=False))
#step31 = vds.count_variants()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# step3.2: population specific call rate > .98
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print "starting step 3.2..."
l = set(vds.samples_table().select('sa.meta.population').collect())
pops = [i.population for i in l]

for p in pops:
     vds = vds.annotate_variants_expr('va.callRate.' + p + ' = gs.filter(g => sa.meta.population == "' + p + '").fraction(g => isDefined(g.GT))')

for p in pops:
     vds = vds.filter_variants_expr('va.callRate.' + p + ' > ' + mincallrate)

#step32 = vds.count_variants()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# step4: population specific HWE P-value > 1e-06
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print "starting step 4..."
for p in pops:
     vds = vds.annotate_variants_expr('va.hwe.' + p + ' = '
     'let nHomRef = gs.filter(g => sa.meta.population == "' + p + '" && g.GT.isHomRef()).count().toInt32() and '
     'nHet = gs.filter(g => sa.meta.population == "' + p + '" && g.GT.isHet()).count().toInt32() and '
     'nHomVar = gs.filter(g => sa.meta.population == "' + p + '" && g.GT.isHomVar()).count().toInt32() in '
     'hwe(nHomRef, nHet, nHomVar)')


for p in pops:
     vds = vds.filter_variants_expr('va.hwe.' + p + '.pHWE > ' + hwep)

#step4 = vds.count_variants()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write output file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print "writing output files..."
with hadoop_write(variant_qc_stats_file) as f:
     for k, v in {'step0':step0, 'step1':step1, 'step2':step2, 'step31':step31, 'step32':step32, 'step4':step4}.iteritems():
          print>>f, k + '\t' + str(v)


#vds.annotate_global('global.stats', [step0, step1, step2, step3, step4, numgsdpo400, fracgsabhomref, fracgsabhomvar, fracgsabhet], expr.types.TArray(expr.types.TFloat64()))

vds.write(qced_file)


# print runtime

stop = timeit.default_timer()
print "runtime: " + str(stop - start) + " seconds"



