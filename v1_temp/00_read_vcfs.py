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
# if VCF
rawvcf_file = 'gs://wgspd-wgs-v1/raw/vcf/wgspd_wgs_v1.vcf.bgz/'
# if VDS
rawvds_input_file = 'gs://wgspd-wgs-v1_temp/raw/vds/wgspd_wgs_v1.vds/'

# output
rawvds_file = 'gs://wgspd-wgs-v1_temp/raw/vds/' + chrom + '/wgspd_wgs_v1.vds'
splitvds_file = 'gs://wgspd-wgs-v1_temp/raw/vds/' + chrom + '/wgspd_wgs_v1_split.vds'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#filelist = [rawvcf_file + str(i) + '_preQC.vcf.gz' for i in range(20360+1)]
# if VCF
#vds = hc.import_vcf(rawvcf_file + '*bgz', force_bgz=True)
# if VDS
vds = hc.read(rawvds_input_file)
print "Import done!"

if chrom == 'chr20':
  vds = vds.filter_variants_expr('v.contig == "20"', keep=True)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# split multi
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
splitvds = vds.split_multi_generic(
   'va.aIndex = aIndex, va.wasSplit = wasSplit',
   '''g = let
         newgt = downcode(g.GT, aIndex) and
         newad = if (isDefined(g.AD))
             let sum = g.AD.sum() and adi = g.AD[aIndex] in [sum - adi, adi]
           else
             NA: Array[Int] and
         newpl = if (isDefined(g.PL))
             range(3).map(i => range(g.PL.length).filter(j => downcode(Call(j), aIndex) == Call(i)).map(j => g.PL[j]).min())
           else
             NA: Array[Int] and
         newgq = gqFromPL(newpl)
    in { GT: newgt, AD: newad, DP: g.DP, GQ: newgq, PL: newpl }''')





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write vds
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

splitvds.write(splitvds_file, overwrite=True)
vds.write(rawvds_file, overwrite=True)



# print runtime

stop = timeit.default_timer()

print "runtime: " + str(stop - start) + " seconds"

















