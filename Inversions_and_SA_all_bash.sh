############### ############### ############### 
############### Inversions and SA #############
############### ############### ############### 


###############  Using git ############### 

#Clone github repo
git clone https://github.com/filipluca/Inversions_and_SA.git

#Modify readme
echo "Winter project" >> README.md

#Add R script to github repo
git add -A; git commit -m "R commands, 22 July"; git push


############### ############### ###############
############### Inversions and candidate genes overlap #############
############### ############### ###############

#See R code


############### ############### ###############
############### Inversions and coupling LD in LHM #############
############### ############### ###############

## ANTAGONISM INDEX

## GWAS on antagonism index
~/Downloads/ldak5.mac --linear ant --pheno pheno.txt --bfile f3c.lhm.snp --grm kinsm_no_outlier --mpheno 6; 
rm ant.pvalues; rm ant.coeff; rm ant.progress; rm ant.score; rm ant.summaries;

## Extract top 2,372 SNPs
sort -k 7g ant.assoc | sed -n 2,2373p > ant.top
cut -d " " -f2 ant.top > ant.top.predictors;
rm ant.assoc; 

## Calculate LD in plink for subset of predictors that are candidates, regardless of distance
## Subset all predictors in plink to keep only top candidates
~/Downloads/plink_mac_20190304/plink --bfile f3c.lhm.snp --extract ant.top.predictors --make-bed --out ant.top;
## Calculate r for each pair of SNP
~/Downloads/plink_mac_20190304/plink --r --bfile ant.top --ld-window 2373 --ld-window-kb 100000 --out ant.top
rm ant.top.log; rm ant.top.nosex;rm ant.top.bed; rm ant.top.fam; rm ant.top.bim;

## Calculate LD in plink for subset of predictors that are candidates AND missense, regardless of distance
## Subset all predictors in plink to keep only top candidates AND missense
sort ant.top.predictors > ant.top.predictors.sorted; sort missense_SNPs.txt > missense_SNPs.txt.sorted; comm -12 ant.top.predictors.sorted missense_SNPs.txt.sorted > ant.top.predictors.missense
~/Downloads/plink_mac_20190304/plink --bfile f3c.lhm.snp --extract ant.top.predictors.missense --make-bed --out ant.top.missense;
## Calculate r for each pair of SNP
~/Downloads/plink_mac_20190304/plink --r --bfile ant.top.missense --ld-window 2373 --ld-window-kb 100000 --out ant.top.missense;
rm ant.top.missense.log; rm ant.top.missense.nosex;rm ant.top.missense.bed; rm ant.top.missense.fam; rm ant.top.missense.bim;


## CONCORDANT INDEX

## GWAS on concordant index
~/Downloads/ldak5.mac --linear conc --pheno pheno.txt --bfile f3c.lhm.snp --grm kinsm_no_outlier --mpheno 5; 
rm conc.pvalues; rm conc.coeff; rm conc.progress; rm conc.score; rm conc.summaries;

## Extract top 2,372 SNPs
sort -k 7g conc.assoc | sed -n 2,2373p > conc.top;
cut -d " " -f2 conc.top > conc.top.predictors;
rm conc.assoc; 

## Calculate LD in plink for subset of predictors that are candidates, regardless of distance
## Subset all predictors in plink to keep only top candidates
~/Downloads/plink_mac_20190304/plink --bfile f3c.lhm.snp --extract conc.top.predictors --make-bed --out conc.top;
## Calculate r for each pair of SNP
~/Downloads/plink_mac_20190304/plink --r --bfile conc.top --ld-window 2373 --ld-window-kb 100000 --out conc.top;
rm conc.top.log; rm conc.top.nosex;rm conc.top.bed; rm conc.top.fam; rm conc.top.bim;


## Calculate LD in plink for subset of predictors that are candidates AND missense, regardless of distance
## Subset all predictors in plink to keep only top candidates AND missense
sort conc.top.predictors > conc.top.predictors.sorted; sort missense_SNPs.txt > missense_SNPs.txt.sorted; comm -12 conc.top.predictors.sorted missense_SNPs.txt.sorted > conc.top.predictors.missense
~/Downloads/plink_mac_20190304/plink --bfile f3c.lhm.snp --extract conc.top.predictors.missense --make-bed --out conc.top.missense;
## Calculate r for each pair of SNP
~/Downloads/plink_mac_20190304/plink --r --bfile conc.top.missense --ld-window 2373 --ld-window-kb 100000 --out conc.top.missense;
rm conc.top.missense.log; rm conc.top.missense.nosex; rm conc.top.missense.bed; rm conc.top.missense.fam; rm conc.top.missense.bim;





############### ############### ###############
############### Inversions and coupling LD in the DGRP/DPGP3 #############
############### ############### ###############

## DGRP VCF, reference 6
#Downloaded from https://zenodo.org/record/837947#.XVDoRZMzbUI

#Subset dgrp2_dm6_dbSNP.vcf.bim file to only include positions that have effect size information from GWAS
#See R code for getting 'dgrp_overlapping_predictors.txt' file
~/Downloads/plink_mac_20190304/plink --bfile dgrp2_dm6_dbSNP.vcf --extract dgrp_overlapping_predictors.txt --make-bed --out dgrp2_dm6_dbSNP.overlapping_gwas --allow-extra-chr

## The SNP IDs are nevertheless annoyingly formatted. To modify them, do the following:

#First convert bed to vcf
~/Downloads/plink_mac_20190304/plink --bfile dgrp2_dm6_dbSNP.overlapping_gwas --recode vcf --out dgrp2_dm6_dbSNP.overlapping_gwas --allow-extra-chr

#Then modify SNP IDs of vcf file
grep "^#" dgrp2_dm6_dbSNP.overlapping_gwas.vcf > dgrp2_dm6_dbSNP.overlapping_gwas.header;
grep -v "^#" dgrp2_dm6_dbSNP.overlapping_gwas.vcf | cut -f4- > dgrp2_dm6_dbSNP.overlapping_gwas.body;
cut -f1,2,3 dgrp2_dm6_dbSNP.overlapping_gwas.vcf > dgrp2_dm6_dbSNP.overlapping_gwas.cols123;
#See R code for getting 'dgrp2_dm6_dbSNP.overlapping_gwas.cols123new' file
paste dgrp2_dm6_dbSNP.overlapping_gwas.cols123new dgrp2_dm6_dbSNP.overlapping_gwas.body > dgrp2_dm6_dbSNP.overlapping_gwas.newbody;
cat dgrp2_dm6_dbSNP.overlapping_gwas.header dgrp2_dm6_dbSNP.overlapping_gwas.newbody > dgrp2_dm6_dbSNP.overlapping_gwas.new.vcf;
rm dgrp2_dm6_dbSNP.overlapping_gwas.header; rm dgrp2_dm6_dbSNP.overlapping_gwas.body; rm dgrp2_dm6_dbSNP.overlapping_gwas.newbody; rm dgrp2_dm6_dbSNP.overlapping_gwas.cols123*; rm dgrp2_dm6_dbSNP.overlapping_gwas.vcf;
rm dgrp2_dm6_dbSNP.overlapping_gwas.bed; rm dgrp2_dm6_dbSNP.overlapping_gwas.bim; rm dgrp2_dm6_dbSNP.overlapping_gwas.fam; rm dgrp2_dm6_dbSNP.overlapping_gwas.log; rm dgrp2_dm6_dbSNP.overlapping_gwas.nosex;

#Then convert modified vcf back to bed format
~/Downloads/plink_mac_20190304/plink --vcf dgrp2_dm6_dbSNP.overlapping_gwas.new.vcf --make-bed --out dgrp2_dm6_dbSNP.overlapping_gwas.new
rm dgrp2_dm6_dbSNP.overlapping_gwas.new.log; rm dgrp2_dm6_dbSNP.overlapping_gwas.new.nosex;



## ZI VCF, reference 6
#Downloaded as 'r6_modified_no_dp_filter.vcf', i.e. snp-sites vcf from Nexus, converted to r6 coordinates, filtered for biallelic sites


## Convert vcf for each chromosome arm to bim format

#Subset dgrp2_dm6_dbSNP.vcf.bim file to only include positions that have effect size information from GWAS
#See R code for getting 'dgrp_overlapping_predictors.txt' file
~/Downloads/plink_mac_20190304/plink --bfile dgrp2_dm6_dbSNP.vcf --extract dgrp_overlapping_predictors.txt --make-bed --out dgrp2_dm6_dbSNP.overlapping_gwas --allow-extra-chr

## The SNP IDs are nevertheless annoyingly formatted. To modify them, do the following:

#First convert bed to vcf
~/Downloads/plink_mac_20190304/plink --bfile dgrp2_dm6_dbSNP.overlapping_gwas --recode vcf --out dgrp2_dm6_dbSNP.overlapping_gwas --allow-extra-chr

#Then modify SNP IDs of vcf file
grep "^#" dgrp2_dm6_dbSNP.overlapping_gwas.vcf > dgrp2_dm6_dbSNP.overlapping_gwas.header;
grep -v "^#" dgrp2_dm6_dbSNP.overlapping_gwas.vcf | cut -f4- > dgrp2_dm6_dbSNP.overlapping_gwas.body;
cut -f1,2,3 dgrp2_dm6_dbSNP.overlapping_gwas.vcf > dgrp2_dm6_dbSNP.overlapping_gwas.cols123;
#See R code for getting 'dgrp2_dm6_dbSNP.overlapping_gwas.cols123new' file
paste dgrp2_dm6_dbSNP.overlapping_gwas.cols123new dgrp2_dm6_dbSNP.overlapping_gwas.body > dgrp2_dm6_dbSNP.overlapping_gwas.newbody;
cat dgrp2_dm6_dbSNP.overlapping_gwas.header dgrp2_dm6_dbSNP.overlapping_gwas.newbody > dgrp2_dm6_dbSNP.overlapping_gwas.new.vcf;
rm dgrp2_dm6_dbSNP.overlapping_gwas.header; rm dgrp2_dm6_dbSNP.overlapping_gwas.body; rm dgrp2_dm6_dbSNP.overlapping_gwas.newbody; rm dgrp2_dm6_dbSNP.overlapping_gwas.cols123*; rm dgrp2_dm6_dbSNP.overlapping_gwas.vcf;
rm dgrp2_dm6_dbSNP.overlapping_gwas.bed; rm dgrp2_dm6_dbSNP.overlapping_gwas.bim; rm dgrp2_dm6_dbSNP.overlapping_gwas.fam; rm dgrp2_dm6_dbSNP.overlapping_gwas.log; rm dgrp2_dm6_dbSNP.overlapping_gwas.nosex;

#Then convert modified vcf back to bed format
~/Downloads/plink_mac_20190304/plink --vcf dgrp2_dm6_dbSNP.overlapping_gwas.new.vcf --make-bed --out dgrp2_dm6_dbSNP.overlapping_gwas.new
rm dgrp2_dm6_dbSNP.overlapping_gwas.new.log; rm dgrp2_dm6_dbSNP.overlapping_gwas.new.nosex;









## ANTAGONISM INDEX
## Calculate LD in plink for subset of predictors that are candidates, regardless of distance
## Subset all predictors in plink to keep only top candidates
~/Downloads/plink_mac_20190304/plink --bfile dgrp2_dm6_dbSNP.overlapping_gwas.new --extract ant.top.predictors --make-bed --out ant.top.dgrp;
~/Downloads/plink_mac_20190304/plink --r --bfile ant.top.dgrp --ld-window 2373 --ld-window-kb 100000 --out ant.top.dgrp
rm ant.top.dgrp.log; rm ant.top.dgrp.nosex;rm ant.top.dgrp.bed; rm ant.top.dgrp.fam; rm ant.top.dgrp.bim;

#Missense SNPs
~/Downloads/plink_mac_20190304/plink --bfile dgrp2_dm6_dbSNP.overlapping_gwas.new --extract ant.top.predictors.missense --make-bed --out ant.top.missense.dgrp;
~/Downloads/plink_mac_20190304/plink --r --bfile ant.top.missense.dgrp --ld-window 2373 --ld-window-kb 100000 --out ant.top.missense.dgrp;
rm ant.top.missense.dgrp.log; rm ant.top.missense.dgrp.nosex;rm ant.top.missense.dgrp.bed; rm ant.top.missense.dgrp.fam; rm ant.top.missense.dgrp.bim;




## CONCORDANT INDEX
~/Downloads/plink_mac_20190304/plink --bfile dgrp2_dm6_dbSNP.overlapping_gwas.new --extract conc.top.predictors --make-bed --out conc.top.dgrp;
~/Downloads/plink_mac_20190304/plink --r --bfile conc.top.dgrp --ld-window 2373 --ld-window-kb 100000 --out conc.top.dgrp
rm conc.top.dgrp.log; rm conc.top.dgrp.nosex;rm conc.top.dgrp.bed; rm conc.top.dgrp.fam; rm conc.top.dgrp.bim;

#Missense SNPs
~/Downloads/plink_mac_20190304/plink --bfile dgrp2_dm6_dbSNP.overlapping_gwas.new --extract conc.top.predictors.missense --make-bed --out conc.top.missense.dgrp;
~/Downloads/plink_mac_20190304/plink --r --bfile conc.top.missense.dgrp --ld-window 2373 --ld-window-kb 100000 --out conc.top.missense.dgrp;
rm conc.top.missense.dgrp.log; rm conc.top.missense.dgrp.nosex;rm conc.top.missense.dgrp.bed; rm conc.top.missense.dgrp.fam; rm conc.top.missense.dgrp.bim;









## 1000 PERMUTED ANTAGONISM INDICES

for i in 1000_shuffled_phenos/*; 
do ~/Downloads/ldak5.mac --linear $i --pheno $i --bfile f3c.lhm.snp --grm kinsm_no_outlier --mpheno 6; 
rm $i.pvalues; rm $i.coeff; rm $i.progress; rm $i.score; rm $i.summaries;
sort -k 7g $i.assoc | sed -n 2,2373p > $i.top;
cut -d " " -f2 $i.top > $i.top.predictors;
rm $i.assoc; 
#LHM
~/Downloads/plink_mac_20190304/plink --bfile f3c.lhm.snp --extract $i.top.predictors --make-bed --out $i.top;
~/Downloads/plink_mac_20190304/plink --r --bfile $i.top --ld-window 2373 --ld-window-kb 100000 --out $i.top;
rm $i.top.log; rm $i.top.nosex; rm $i.top.bed; rm $i.top.fam; rm $i.top.bim;
#DGRP
~/Downloads/plink_mac_20190304/plink --bfile dgrp2_dm6_dbSNP.overlapping_gwas.new --extract $i.top.predictors --make-bed --out $i.top.dgrp;
~/Downloads/plink_mac_20190304/plink --r --bfile $i.top.dgrp --ld-window 2373 --ld-window-kb 100000 --out $i.top.dgrp;
rm $i.top.dgrp.log; rm $i.top.dgrp.nosex; rm $i.top.dgrp.bed; rm $i.top.dgrp.fam; rm $i.top.dgrp.bim;
done






############### Clump SNPs ###############

#Import top 226 antagonistic SNPs (=clumped_highly_associated_r0.4_10kb.clumped)

#Clump top 2,373 concordant SNPs, as clumped_highly_associated_r0.4_10kb.clumped
#Maximum P-value for candidate concordant SNPs is 0.0025863
~/Downloads/plink --noweb --bfile f3c.lhm.snp --clump conc.assoc2 --clump-p1 0.0025863 --clump-p2 0.1 --clump-kb 10 --clump-r2 0.4 --clump-field Wald_P --out clumped_highly_associated_concordant_r0.4_10kb




