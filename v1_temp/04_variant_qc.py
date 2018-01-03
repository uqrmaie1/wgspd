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
lcr_file = 'gs://hail-common/LCR.interval_list'
sample_qc_file = 'gs://wgspd-wgs-v2/qc_measures/' + chrom + '/samples_to_keep_after_qc.txt'
sample_sex_fstat_file = 'gs://wgspd-wgs-v2/qc_measures/' + chrom + '/sample_sex_fstat.txt'

# output
qced_file = 'gs://wgspd-wgs-v2/qced/vds/' + chrom + '/wgspd_wgs_v2_split.vds'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define constants
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mincallrate = '0.98'
hwep = '0.000001'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

vds = hc.read(splitvds_file)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# variant qc
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# annotate samples with project IDs

# annotate samples with result from sex imputation

# filter samples
vds = vds.filter_samples_list(sample_qc_file)
vds = vds.filter_variants_intervals(lcr_file, keep = False)


vds = (vds.variant_qc().annotate_variants_expr(
     """
     va.step0 = va.qc.AC > 0,
     va.step1 = va.qc.AC > 0 && va.filters.isEmpty(),
     va.step2 = va.qc.AC > 0 && va.filters.isEmpty() && va.info.QD > 4
     """))

step0 = vds.query_variants("variants.filter(v => va.step0).count()")
step1 = vds.query_variants("variants.filter(v => va.step1).count()")
step2 = vds.query_variants("variants.filter(v => va.step2).count()")

vds = vds.filter_variants_expr('va.qc.AC > 0 && va.filters.isEmpty() && va.info.QD > 4', keep=True)

vds = (vds.annotate_variants_expr(
     """
     va.numgsdpo400=gs.filter(g => g.DP > 400).count(),
     va.numgsabhomref=gs.filter(g => (g.GT.isHomRef() && (g.AD[0] / g.DP < 0.9 || g.GQ < 20))).count(),
     va.numgsabhomvar=gs.filter(g => (g.GT.isHomVar() && (g.AD[1] / g.DP < 0.9 || g.PL[0] < 20))).count(),
     va.numgsabhet=gs.filter(g => (g.GT.isHet() && ( (g.AD[0] + g.AD[1]) / g.DP < 0.9 || g.AD[1] / g.DP < 0.20 || g.PL[0] < 20 || (!sa.imputesex.isFemale && ( v.inXNonPar || v.inYNonPar )) )) || (v.inYNonPar && sa.imputesex.isFemale)).count(),
     va.allgs=gs.count(),
     va.allgshomref=gs.filter(g => g.GT.isHomRef()).count(),
     va.allgshomvar=gs.filter(g => g.GT.isHomVar()).count(),
     va.allgshet=gs.filter(g => g.GT.isHet()).count()
     """))

numgsdpo400 = vds.query_variants("variants.map(v => va.numgsdpo400).sum()/variants.map(v => va.allgs).sum()")
fracgsabhomref = vds.query_variants("variants.map(v => va.numgsabhomref).sum()/variants.map(v => va.allgshomref).sum()")
fracgsabhomvar = vds.query_variants("variants.map(v => va.numgsabhomvar).sum()/variants.map(v => va.allgshomvar).sum()")
fracgsabhet = vds.query_variants("variants.map(v => va.numgsabhet).sum()/variants.map(v => va.allgshet).sum()")

vds = (vds.filter_genotypes(
     """
     g.DP > 400 ||
     (g.GT.isHomRef() && (g.AD[0] / g.DP < 0.9 || g.GQ < 20)) ||
     (g.GT.isHomVar() && (g.AD[1] / g.DP < 0.9 || g.PL[0] < 20)) ||
     (g.GT.isHet() && ( (g.AD[0] + g.AD[1]) / g.DP < 0.9 || g.AD[1] / g.DP < 0.20 || g.PL[0] < 20 || (!sa.imputesex.isFemale && ( v.inXNonPar || v.inYNonPar )) )) || (v.inYNonPar && sa.imputesex.isFemale)
     """, keep = False))


# project stratified HWE and isCalled

# step3: variants with callrate > 0.98 in all projects
# step4: variants which pass step3 && pHWE > 0.000001 in all projects  

#numvar_step3 = variants.filter(v => va.step3).count(),
#numvar_step4 = variants.filter(v => va.step4).count(),
#callRate098 = variants.filter(v => va.callRate < 0.98).count(),
#MAC = variants.filter(v => va.qc.AC == 0).count()


vds = vds.filter_variants_expr('va.qc.pHWE > '+hwep+' && va.qc.AC > 0 && va.qc.callRate > '+mincallrate+' && va.callRate_GPC_afr > '+mincallrate+' && va.callRate_Fin > '+mincallrate+'  && va.callRate_GPC_nfe > '+mincallrate+'  && va.callRate_Est > '+mincallrate+'  && va.callRate_other > '+mincallrate+' && va.pHWE_GPC_afr > '+hwep+' && va.pHWE_Fin > '+hwep+' && va.pHWE_GPC_nfe > '+hwep+' && va.pHWE_Est > '+hwep+' && va.pHWE_other > '+hwep+' && v.contig != "X" ', keep=True)

vds.annotate_global('global.stats', [step0, step1, step2, step3, step4, numgsdpo400, fracgsabhomref, fracgsabhomvar, fracgsabhet], hail.expr.types.TArray(hail.expr.types.TFloat64()))

vds.write(qced_file)


# print runtime

stop = timeit.default_timer()
print "runtime: " + str(stop - start) + " seconds"



