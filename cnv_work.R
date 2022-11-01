library(rCNV)
library(tidyverse)
#https://piyalkarum.github.io/rCNV/articles/rCNV.html tutorial
#removing multiallelic: bcftools view -m2 -M2 -v snps ~/Documents/neanderthal/chag/chag_vcfs/chr22.noRB.vcf > ~/Documents/chr22_chag_biallelic.vcf


vcf = readVCF(vcf.file.path = "chr22_chag_biallelic.vcf", verbose = F)
  #says that it can be .gz, but it LIES 
  #imports the VCF file - for this, the file = chag_22

vcf$vcf[1:10, 1:10]
  #prints the first 10 rows and cols of imported vcf

h.table = hetTgen(vcf = vcf, info.type = "AD", verbose = F)
  #issue is that there is circular logic
  #w/ normal VCF, can't remove multiallelic loci w/out h.table
  #can get AD, but not missing data with this - bcftools removes
  #h/ever, can't do with norm vcf, as there are biallelic SNPz
#so this is where i can get to with norm VCFs

#####missingness removed in vcf in terminal
setwd("~/Documents/neanderthal/chag/chag_vcfs")

vcf = readVCF(vcf.file.path = "chr22_nomissing.vcf", verbose = F)

ad.tab = hetTgen(vcf = vcf, info.type = "AD")
  #gives warning of: vcf file contains multi-allelic variants: only bi-allelic SNPs allowed
  #Use maf() to remove non-bi-allilic snps but not an error
ad.tab = maf(h.table = ad.tab, AD=T, verbose = T)
  #using this, get warning: NAs introduced by coercion

gt<-hetTgen(vcf = vcf,info.type = "GT")
ad.tab<-ad.correct(ad.tab,gt.table = gt)

#normalize depth table with cpm.normal()
ad.nor<-cpm.normal(ad.tab,method="MedR")
ad.nor[1:6,1:6]

A.info<-allele.info(X=ad.tab,x.norm = ad.nor, plot.allele.cov = TRUE)
# X is the corrected non-normalized allele depth table and x.norm is the normalized allele depth table
#Error: Error in txtProgressBar(min = 0, max = pb_Total, width = 50, style = 3) : 
#must have 'max' > 'min'

head(A.info)


#####biallelic only vcf
vcf = readVCF(vcf.file.path = "chr22_chag_biallelic.vcf", verbose = F)

h.table = hetTgen(vcf = vcf, info.type = "AD", verbose = F)

gt = hetTgen(vcf = vcf, info.type = "GT", verbose = F)

ad.tab = ad.correct(h.table, gt.table = gt)
  #geno table is needed for correcting mismatches, as above

ad.nor = cpm.normal(ad.tab, method = "MedR")
ad.nor[1:6,1:6]
  #normalise depth table
  #cannot print with above command - incorrect dimensions error

#to detect cnvs, generating allele information table is needed which in turn needs:
#allele ratios, proportions of hom/hets per SNP, depth ratios, z-score, chi square, f-excess of het

##
mss = get.miss(data = vcf, verbose = F, plot = F)
