## Set working directory
setwd("Documents/inversions_and_SA/")

#################################################################
## Subset dgrp2 dataset to only include sites that are in GWAS
#################################################################

dgrp2.bim <- read.table("dgrp2_dm6_dbSNP.vcf.bim")
dgrp2.bim.snps <- subset(dgrp2.bim,V5 %in% c("A","T","C","G") & V6 %in% c("A","T","C","G"))
rm(dgrp2.bim)

#Make column which names SNPs as they are named in LHm dataset (Chromplink_Basepair)
levels(dgrp2.bim.snps$V1)[1] <- "X"
levels(dgrp2.bim.snps$V1) <- c(5,NA,1,2,3,4,NA,NA,NA,NA)
dgrp2.bim.snps$Predictor <- paste(dgrp2.bim.snps$V1,dgrp2.bim.snps$V4,sep="_")

#Subset dgrp SNPs to only include overlapping GWAS sites
assoc <- read.table("Assoc_Data.txt",h=T)
dgrp2.bim.snps.overlapping <- subset(dgrp2.bim.snps,Predictor %in% assoc$Predictor) 
rm(dgrp2.bim.snps)
dgrp2.bim.snps.overlapping <- merge(dgrp2.bim.snps.overlapping,assoc[c("Predictor","A1","A2")],by="Predictor",all.x=T)
dgrp2.bim.snps.overlapping$V5 <- as.character(dgrp2.bim.snps.overlapping$V5)
dgrp2.bim.snps.overlapping$V6 <- as.character(dgrp2.bim.snps.overlapping$V6)
dgrp2.bim.snps.overlapping$A1 <- as.character(dgrp2.bim.snps.overlapping$A1)
dgrp2.bim.snps.overlapping$A2 <- as.character(dgrp2.bim.snps.overlapping$A2)
dgrp2.bim.snps.overlapping$Same_alleles <- with(dgrp2.bim.snps.overlapping,ifelse((V5==A1 & V6==A2) | (V5==A2 & V6==A1),1,0))
dgrp2.bim.snps.overlapping.same.alleles <- subset(dgrp2.bim.snps.overlapping,Same_alleles==1)
rm(dgrp2.bim.snps.overlapping)
#Write out list of predictors to potentially calculate LD over
write.table(dgrp2.bim.snps.overlapping.same.alleles$V2,"dgrp_overlapping_predictors.txt",quote=F,col.names=F,row.names=F)


#################################################################
## Rename dgrp SNP IDs
#################################################################

dgrp2.snp.ids <- read.table("dgrp2_dm6_dbSNP.overlapping_gwas.cols123")
levels(dgrp2.snp.ids$V1) <- c(5,1,2,3,4)
dgrp2.snp.ids$V3 <- paste(dgrp2.snp.ids$V1,dgrp2.snp.ids$V2,sep="_")
write.table(dgrp2.snp.ids,"dgrp2_dm6_dbSNP.overlapping_gwas.cols123new",sep="\t",quote=F,row.names=F,col.names = F)


#################################################################
## Write missense SNPs
#################################################################

assoc <- read.table("Assoc_Data.txt",h=T)
write.table(assoc$Predictor[assoc$Consequence=="missense_variant"],"missense_SNPs.txt",quote=F,row.names=F,col.names=F)

#################################################################
## Import inversions
#################################################################

## Import inversion dataframe
inv <- read.table("inversions_plus_genomic_coordinates.txt",head=T)

#Make column with definitive inversion breakpoints, i.e. genomic if possible, but otherwise cyto->seq converted
inv$Start_definitive <- with(inv,ifelse(!is.na(inv$Start_r6),Start_r6,Start_converted_midpoint))
inv$End_definitive <- with(inv,ifelse(!is.na(inv$End_r6),End_r6,End_converted_midpoint))

#Remove extraneous columns
inv <- inv[c("Name","Common","Rare","Chromosome","Start_definitive","End_definitive")]

#Inversion length and total chromosome length
inv$Length <- inv$End_definitive-inv$Start_definitive
#Chromosome sizes downloaded from http://genome.ucsc.edu/cgi-bin/hgTracks?db=dm6&chromInfoPage=
inv$Max <- ifelse(inv$Chromosome=="2L",23513712,ifelse(inv$Chromosome=="2R",25286936,ifelse(inv$Chromosome=="3L",28110227,ifelse(inv$Chromosome=="3R",32079331,ifelse(inv$Chromosome=="X",23542271,NA)))))
#Remove pericentric inversions
inv <- subset(inv,Chromosome %in% c("2L","2R","3L","3R","X"))

#Subset total inversion dataframe to include only common paracentric inversions 
inv.common <- subset(inv,Common==1)
inv.common <- inv.common[order(inv.common$Chromosome),]
#Subset total inversion dataframe to include only rare paracentric inversions
inv.rare <- subset(inv,Rare==1)
inv.Chakraborty <- subset(inv,is.na(Common))

#################################################################
## Import candidate SNP info
#################################################################

## ANTAGONISM INDEX #####

## Import assoc info
ant.top.assoc <- read.table("ant.top")
#Order by chromosome and position
ant.top.assoc <- ant.top.assoc[order(ant.top.assoc$V1,ant.top.assoc$V3),]
#Rename important columns
names(ant.top.assoc)[1:8] <- c("Chromosome","Predictor","Basepair","A1","A2","Wald_Stat","Wald_P","Effect")
#Make chromosome column non-numeric
ant.top.assoc$Chromosome_numeric <- ant.top.assoc$Chromosome
ant.top.assoc$Chromosome <- as.factor(ant.top.assoc$Chromosome)
levels(ant.top.assoc$Chromosome) <- c("2L","2R","3L","3R","X")

#Check if minor allele in LHM is minor allele in DGRP
ant.top.assoc <- merge(ant.top.assoc,dgrp2.bim.snps.overlapping.same.alleles[c("Predictor","V5","V6")],by="Predictor",all.x=T)
names(ant.top.assoc)[17:18] <- c("A1_DGRP","A2_DGRP")
ant.top.assoc$Same_minor_alleles <- with(ant.top.assoc,ifelse(A1==A1_DGRP & A2==A2_DGRP,1,ifelse(A1==A2_DGRP & A2==A1_DGRP,0,NA)))

#Remove extraneous columns
ant.top.assoc <- ant.top.assoc[,c("Chromosome","Basepair","Predictor","Effect","A1","A2","A1_DGRP","A2_DGRP","Same_minor_alleles")]


## ANTAGONISM INDEX, MISSENSE ONLY #####

## Import assoc info, with only missense candidate SNPS
ant.top.missense.assoc <- subset(assoc,Consequence=="missense_variant" & Cand_parametric==1)

#Check if minor allele in LHM is minor allele in DGRP
ant.top.missense.assoc <- merge(ant.top.missense.assoc,dgrp2.bim.snps.overlapping.same.alleles[c("Predictor","V5","V6")],by="Predictor",all.x=T)
names(ant.top.missense.assoc)[91:92] <- c("A1_DGRP","A2_DGRP")
ant.top.missense.assoc$Same_minor_alleles <- with(ant.top.missense.assoc,ifelse(A1==A1_DGRP & A2==A2_DGRP,1,ifelse(A1==A2_DGRP & A2==A1_DGRP,0,NA)))

#Remove extraneous columns
ant.top.missense.assoc <- ant.top.missense.assoc[,c("Chromosome","Basepair","Predictor","Effect","A1","A2","A1_DGRP","A2_DGRP","Same_minor_alleles")]



## CONCORDANT INDEX #####


## Import assoc info, with only missense candidate SNPS
conc.top.missense.assoc <- subset(assoc,Consequence=="missense_variant" & Cand_concordant==1)

#Check if minor allele in LHM is minor allele in DGRP
conc.top.missense.assoc <- merge(conc.top.missense.assoc,dgrp2.bim.snps.overlapping.same.alleles[c("Predictor","V5","V6")],by="Predictor",all.x=T)
names(conc.top.missense.assoc)[91:92] <- c("A1_DGRP","A2_DGRP")
conc.top.missense.assoc$Same_minor_alleles <- with(conc.top.missense.assoc,ifelse(A1==A1_DGRP & A2==A2_DGRP,1,ifelse(A1==A2_DGRP & A2==A1_DGRP,0,NA)))

#Remove extraneous columns
conc.top.missense.assoc <- conc.top.missense.assoc[,c("Chromosome","Basepair","Predictor","Effect","A1","A2","A1_DGRP","A2_DGRP","Same_minor_alleles")]


#################################################################
## If candidate SNP falls within inversion, record it
#################################################################

## ANTAGONISM INDEX #####

## If common inversion falls within a SNP (j), the overlap column is the name of the inversion
for (i in 1:21){
  ant.top.assoc[[paste("is_common_inv",i,sep = "_")]] <- 0
for (j in 1:nrow(ant.top.assoc)){
  #Common inversion
    if(inv.common$Chromosome[i]==ant.top.assoc$Chromosome[j]){
      ant.top.assoc[[paste("is_common_inv",i,sep = "_")]][j] <- ifelse((ant.top.assoc$Basepair[j]>=inv.common$Start_definitive[i] & ant.top.assoc$Basepair[j]<=inv.common$End_definitive[i]),1,ant.top.assoc[[paste("is_common_inv",i,sep = "_")]][j])
    }
}
  print(i)
}
## If rare inversion falls within a SNP (j), the overlap column is the name of the inversion
for (i in 1:318){
  ant.top.assoc[[paste("is_rare_inv",i,sep = "_")]] <- 0
  for (j in 1:nrow(ant.top.assoc)){
    #Common inversion
    if(inv.rare$Chromosome[i]==ant.top.assoc$Chromosome[j]){
      ant.top.assoc[[paste("is_rare_inv",i,sep = "_")]][j] <- ifelse((ant.top.assoc$Basepair[j]>=inv.rare$Start_definitive[i] & ant.top.assoc$Basepair[j]<=inv.rare$End_definitive[i]),1,ant.top.assoc[[paste("is_rare_inv",i,sep = "_")]][j])
    }
  }
  print(i)
}

## If Chakraborty inversion falls within a SNP (j), the overlap column is the name of the inversion
for (i in 1:27){
  ant.top.assoc[[paste("is_Chakraborty_inv",i,sep = "_")]] <- 0
  for (j in 1:nrow(ant.top.assoc)){
    #Common inversion
    if(inv.Chakraborty$Chromosome[i]==ant.top.assoc$Chromosome[j]){
      ant.top.assoc[[paste("is_Chakraborty_inv",i,sep = "_")]][j] <- ifelse((ant.top.assoc$Basepair[j]>=inv.Chakraborty$Start_definitive[i] & ant.top.assoc$Basepair[j]<=inv.Chakraborty$End_definitive[i]),1,ant.top.assoc[[paste("is_Chakraborty_inv",i,sep = "_")]][j])
    }
  }
  print(i)
}

## ANTAGONISM INDEX, ONLY MISSENSE SNPS #####

## If common inversion falls within a SNP (j), the overlap column is the name of the inversion
for (i in 1:21){
  ant.top.missense.assoc[[paste("is_common_inv",i,sep = "_")]] <- 0
  for (j in 1:nrow(ant.top.missense.assoc)){
    #Common inversion
    if(inv.common$Chromosome[i]==ant.top.missense.assoc$Chromosome[j]){
      ant.top.missense.assoc[[paste("is_common_inv",i,sep = "_")]][j] <- ifelse((ant.top.missense.assoc$Basepair[j]>=inv.common$Start_definitive[i] & ant.top.missense.assoc$Basepair[j]<=inv.common$End_definitive[i]),1,ant.top.missense.assoc[[paste("is_common_inv",i,sep = "_")]][j])
    }
  }
  print(i)
}
## If rare inversion falls within a SNP (j), the overlap column is the name of the inversion
for (i in 1:318){
  ant.top.missense.assoc[[paste("is_rare_inv",i,sep = "_")]] <- 0
  for (j in 1:nrow(ant.top.missense.assoc)){
    #Common inversion
    if(inv.rare$Chromosome[i]==ant.top.missense.assoc$Chromosome[j]){
      ant.top.missense.assoc[[paste("is_rare_inv",i,sep = "_")]][j] <- ifelse((ant.top.missense.assoc$Basepair[j]>=inv.rare$Start_definitive[i] & ant.top.missense.assoc$Basepair[j]<=inv.rare$End_definitive[i]),1,ant.top.missense.assoc[[paste("is_rare_inv",i,sep = "_")]][j])
    }
  }
  print(i)
}

## If Chakraborty inversion falls within a SNP (j), the overlap column is the name of the inversion
for (i in 1:27){
  ant.top.missense.assoc[[paste("is_Chakraborty_inv",i,sep = "_")]] <- 0
  for (j in 1:nrow(ant.top.missense.assoc)){
    #Common inversion
    if(inv.Chakraborty$Chromosome[i]==ant.top.missense.assoc$Chromosome[j]){
      ant.top.missense.assoc[[paste("is_Chakraborty_inv",i,sep = "_")]][j] <- ifelse((ant.top.missense.assoc$Basepair[j]>=inv.Chakraborty$Start_definitive[i] & ant.top.missense.assoc$Basepair[j]<=inv.Chakraborty$End_definitive[i]),1,ant.top.missense.assoc[[paste("is_Chakraborty_inv",i,sep = "_")]][j])
    }
  }
  print(i)
}

## CONCORDANT INDEX #####

## If common inversion falls within a SNP (j), the overlap column is the name of the inversion
for (i in 1:21){
  conc.top.assoc[[paste("is_common_inv",i,sep = "_")]] <- 0
  for (j in 1:nrow(conc.top.assoc)){
    #Common inversion
    if(inv.common$Chromosome[i]==conc.top.assoc$Chromosome[j]){
      conc.top.assoc[[paste("is_common_inv",i,sep = "_")]][j] <- ifelse((conc.top.assoc$Basepair[j]>=inv.common$Start_definitive[i] & conc.top.assoc$Basepair[j]<=inv.common$End_definitive[i]),1,conc.top.assoc[[paste("is_common_inv",i,sep = "_")]][j])
    }
  }
  print(i)
}

## If rare inversion falls within a SNP (j), the overlap column is the name of the inversion
for (i in 1:318){
  conc.top.assoc[[paste("is_rare_inv",i,sep = "_")]] <- 0
  for (j in 1:nrow(conc.top.assoc)){
    #Common inversion
    if(inv.rare$Chromosome[i]==conc.top.assoc$Chromosome[j]){
      conc.top.assoc[[paste("is_rare_inv",i,sep = "_")]][j] <- ifelse((conc.top.assoc$Basepair[j]>=inv.rare$Start_definitive[i] & conc.top.assoc$Basepair[j]<=inv.rare$End_definitive[i]),1,conc.top.assoc[[paste("is_rare_inv",i,sep = "_")]][j])
    }
  }
  print(i)
}

## If Chakraborty inversion falls within a SNP (j), the overlap column is the name of the inversion
for (i in 1:27){
  conc.top.assoc[[paste("is_Chakraborty_inv",i,sep = "_")]] <- 0
  for (j in 1:nrow(conc.top.assoc)){
    #Common inversion
    if(inv.Chakraborty$Chromosome[i]==conc.top.assoc$Chromosome[j]){
      conc.top.assoc[[paste("is_Chakraborty_inv",i,sep = "_")]][j] <- ifelse((conc.top.assoc$Basepair[j]>=inv.Chakraborty$Start_definitive[i] & conc.top.assoc$Basepair[j]<=inv.Chakraborty$End_definitive[i]),1,conc.top.assoc[[paste("is_Chakraborty_inv",i,sep = "_")]][j])
    }
  }
  print(i)
}

## CONCORDANT INDEX, ONLY MISSENSE SNPS #####

## If common inversion falls within a SNP (j), the overlap column is the name of the inversion
for (i in 1:21){
  conc.top.missense.assoc[[paste("is_common_inv",i,sep = "_")]] <- 0
  for (j in 1:nrow(conc.top.missense.assoc)){
    #Common inversion
    if(inv.common$Chromosome[i]==conc.top.missense.assoc$Chromosome[j]){
      conc.top.missense.assoc[[paste("is_common_inv",i,sep = "_")]][j] <- ifelse((conc.top.missense.assoc$Basepair[j]>=inv.common$Start_definitive[i] & conc.top.missense.assoc$Basepair[j]<=inv.common$End_definitive[i]),1,conc.top.missense.assoc[[paste("is_common_inv",i,sep = "_")]][j])
    }
  }
  print(i)
}

## If rare inversion falls within a SNP (j), the overlap column is the name of the inversion
for (i in 1:318){
  conc.top.missense.assoc[[paste("is_rare_inv",i,sep = "_")]] <- 0
  for (j in 1:nrow(conc.top.missense.assoc)){
    #Common inversion
    if(inv.rare$Chromosome[i]==conc.top.missense.assoc$Chromosome[j]){
      conc.top.missense.assoc[[paste("is_rare_inv",i,sep = "_")]][j] <- ifelse((conc.top.missense.assoc$Basepair[j]>=inv.rare$Start_definitive[i] & conc.top.missense.assoc$Basepair[j]<=inv.rare$End_definitive[i]),1,conc.top.missense.assoc[[paste("is_rare_inv",i,sep = "_")]][j])
    }
  }
  print(i)
}

## If Chakraborty inversion falls within a SNP (j), the overlap column is the name of the inversion
for (i in 1:27){
  conc.top.missense.assoc[[paste("is_Chakraborty_inv",i,sep = "_")]] <- 0
  for (j in 1:nrow(conc.top.missense.assoc)){
    #Common inversion
    if(inv.Chakraborty$Chromosome[i]==conc.top.missense.assoc$Chromosome[j]){
      conc.top.missense.assoc[[paste("is_Chakraborty_inv",i,sep = "_")]][j] <- ifelse((conc.top.missense.assoc$Basepair[j]>=inv.Chakraborty$Start_definitive[i] & conc.top.missense.assoc$Basepair[j]<=inv.Chakraborty$End_definitive[i]),1,conc.top.missense.assoc[[paste("is_Chakraborty_inv",i,sep = "_")]][j])
    }
  }
  print(i)
}

#################################################################
## Coupling LD for each inversion, LHM
#################################################################


## ANTAGONISM INDEX #####

## Import LD info
ant.top.ld <- read.table("ant.top.ld",head=T)

#Obtain allelic effect of locus A ("SNP_A") from assoc info dataframe
tmp1 <- merge(ant.top.ld,ant.top.assoc[c("Predictor","Effect")],by.x="SNP_A",by.y="Predictor",all.x=T)
#Obtain allelic effect of locus B ("SNP_B") from assoc info dataframe
ant.top.ld <- merge(tmp1,ant.top.assoc,by.x="SNP_B",by.y="Predictor",all.x=T)
names(ant.top.ld)[c(8,11)] <- c("Effect_A","Effect_B")
ant.top.ld <- ant.top.ld[order(ant.top.ld$SNP_A,ant.top.ld$SNP_B),]

## Create coupling LD column
#If allelic effects are in the same direction (both positive, or both negative), R is measuring coupling LD, so keep r2 as is; if the signs are opposite, R is measuring repulsion LD, so flip the sign of R
ant.top.ld$Coupling_R <- ifelse((ant.top.ld$Effect_A>0 & ant.top.ld$Effect_B>0) | (ant.top.ld$Effect_A<0 & ant.top.ld$Effect_B<0),ant.top.ld$R,-ant.top.ld$R)

#Double check that plink calculates LD between major/minor alleles, not between reference/alternative!  

## Create distance column
ant.top.ld$Distance <- ant.top.ld$BP_B-ant.top.ld$BP_A

## Plot LD
plot(ant.top.ld$Distance[ant.top.ld$Distance<1000],ant.top.ld$Coupling_R[ant.top.ld$Distance<1000],ylim=c(0,1))


## CONCORDANT INDEX #####

## Import LD info
conc.top.ld <- read.table("conc.top.ld",head=T)

#Obtain allelic effect of locus A ("SNP_A") from assoc info dataframe
tmp1 <- merge(conc.top.ld,conc.top.assoc[c("Predictor","Effect")],by.x="SNP_A",by.y="Predictor",all.x=T)
#Obtain allelic effect of locus B ("SNP_B") from assoc info dataframe
conc.top.ld <- merge(tmp1,conc.top.assoc,by.x="SNP_B",by.y="Predictor",all.x=T)
names(conc.top.ld)[c(8,11)] <- c("Effect_A","Effect_B")
conc.top.ld <- conc.top.ld[order(conc.top.ld$SNP_A,conc.top.ld$SNP_B),]

## Create coupling LD column
#If allelic effects are in the same direction (both positive, or both negative), R is measuring coupling LD, so keep r2 as is; if the signs are opposite, R is measuring repulsion LD, so flip the sign of R
conc.top.ld$Coupling_R <- ifelse((conc.top.ld$Effect_A>0 & conc.top.ld$Effect_B>0) | (conc.top.ld$Effect_A<0 & conc.top.ld$Effect_B<0),conc.top.ld$R,-conc.top.ld$R)

#Double check that plink calculates LD between major/minor alleles, not between reference/alternative!  

## Create distance column
conc.top.ld$Distance <- conc.top.ld$BP_B-conc.top.ld$BP_A

## Plot LD
plot(conc.top.ld$Distance[conc.top.ld$Distance<1000],conc.top.ld$Coupling_R[conc.top.ld$Distance<1000],ylim=c(0,1))


#################################################################
## Coupling LD for each inversion, DGRP
#################################################################

## ANTAGONISM INDEX

## Import LD info
ant.top.dgrp.ld <- read.table("ant.top.dgrp.ld",head=T)

#Obtain allelic effect of locus A ("SNP_A") from assoc info dataframe
tmp1 <- merge(ant.top.dgrp.ld,ant.top.assoc[c("Predictor","Effect","Same_minor_alleles","A1","A2","A1_DGRP","A2_DGRP")],by.x="SNP_A",by.y="Predictor",all.x=T)
#Obtain allelic effect of locus B ("SNP_B") from assoc info dataframe
ant.top.dgrp.ld <- merge(tmp1,ant.top.assoc,by.x="SNP_B",by.y="Predictor",all.x=T)
names(ant.top.dgrp.ld)[c(8:13,16:21)] <- c("Effect_A","Same_minor_alleles_A","A1_A","A2_A","A1_DGRP_A","A2_DGRP_A","Effect_B","A1_B","A2_B","A1_DGRP_B","A2_DGRP_B","Same_minor_alleles_B")
ant.top.dgrp.ld <- ant.top.dgrp.ld[order(ant.top.dgrp.ld$SNP_A,ant.top.dgrp.ld$SNP_B),]

## Create coupling LD column
#If allelic effects are in the same direction (both positive, or both negative), R is measuring coupling LD, so keep r2 as is; if the signs are opposite, R is measuring repulsion LD, so flip the sign of R
ant.top.dgrp.ld$Naive_coupling_R <- ifelse((ant.top.dgrp.ld$Effect_A>0 & ant.top.dgrp.ld$Effect_B>0) | (ant.top.dgrp.ld$Effect_A<0 & ant.top.dgrp.ld$Effect_B<0),ant.top.dgrp.ld$R,-ant.top.dgrp.ld$R)
#If minor alleles are the same in LHm and DGRP, then 'naive' coupling R value is correct; else the sign of R should be flipped
ant.top.dgrp.ld$Coupling_R <- with(ant.top.dgrp.ld,ifelse((Same_minor_alleles_A==1 & Same_minor_alleles_B==1) | (Same_minor_alleles_A==0 & Same_minor_alleles_B==0),R,ifelse((Same_minor_alleles_A==0 & Same_minor_alleles_B==1) | (Same_minor_alleles_A==1 & Same_minor_alleles_B==0),-R,NA)))
#Double check that plink calculates LD between major/minor alleles, not between reference/alternative!!  

## Create distance column
ant.top.dgrp.ld$Distance <- ant.top.dgrp.ld$BP_B-ant.top.dgrp.ld$BP_A

## Calculate r2
ant.top.dgrp.ld$R2 <- (ant.top.dgrp.ld$Coupling_R)^2
##Calculate residual r
ant.top.dgrp.ld$Exp_Coupling_R <- exp(ant.top.dgrp.ld$Coupling_R)
ant.top.dgrp.ld$R_resid[!is.na(ant.top.dgrp.ld$Coupling_R)] <- glm(data=ant.top.dgrp.ld,Exp_Coupling_R~Distance)$residuals

## Are SNPs situated in any common ivnersion
ant.top.dgrp.ld$is_common_inv_any <- ifelse(rowSums(ant.top.dgrp.ld[22:42])>0,1,0)


## ANTAGONISM INDEX, MISSENSE SNPS

## Import LD info
ant.top.missense.dgrp.ld <- read.table("ant.top.missense.dgrp.ld",head=T)

#Obtain allelic effect of locus A ("SNP_A") from assoc info dataframe
tmp1 <- merge(ant.top.missense.dgrp.ld,ant.top.missense.assoc[c("Predictor","Effect","Same_minor_alleles","A1","A2","A1_DGRP","A2_DGRP")],by.x="SNP_A",by.y="Predictor",all.x=T)
#Obtain allelic effect of locus B ("SNP_B") from assoc info dataframe
ant.top.missense.dgrp.ld <- merge(tmp1,ant.top.missense.assoc,by.x="SNP_B",by.y="Predictor",all.x=T)
names(ant.top.missense.dgrp.ld)[c(8:13,16:21)] <- c("Effect_A","Same_minor_alleles_A","A1_A","A2_A","A1_DGRP_A","A2_DGRP_A","Effect_B","A1_B","A2_B","A1_DGRP_B","A2_DGRP_B","Same_minor_alleles_B")
ant.top.missense.dgrp.ld <- ant.top.missense.dgrp.ld[order(ant.top.missense.dgrp.ld$SNP_A,ant.top.missense.dgrp.ld$SNP_B),]

## Create coupling LD column
#If allelic effects are in the same direction (both positive, or both negative), R is measuring coupling LD, so keep r2 as is; if the signs are opposite, R is measuring repulsion LD, so flip the sign of R
ant.top.missense.dgrp.ld$Naive_coupling_R <- ifelse((ant.top.missense.dgrp.ld$Effect_A>0 & ant.top.missense.dgrp.ld$Effect_B>0) | (ant.top.missense.dgrp.ld$Effect_A<0 & ant.top.missense.dgrp.ld$Effect_B<0),ant.top.missense.dgrp.ld$R,-ant.top.missense.dgrp.ld$R)
#If minor alleles are the same in LHm and DGRP, then 'naive' coupling R value is correct; else the sign of R should be flipped
ant.top.missense.dgrp.ld$Coupling_R <- with(ant.top.missense.dgrp.ld,ifelse((Same_minor_alleles_A==1 & Same_minor_alleles_B==1) | (Same_minor_alleles_A==0 & Same_minor_alleles_B==0),R,ifelse((Same_minor_alleles_A==0 & Same_minor_alleles_B==1) | (Same_minor_alleles_A==1 & Same_minor_alleles_B==0),-R,NA)))
#Double check that plink calculates LD between major/minor alleles, not between reference/alternative!!  

## Create distance column
ant.top.missense.dgrp.ld$Distance <- ant.top.missense.dgrp.ld$BP_B-ant.top.missense.dgrp.ld$BP_A

## Calculate r2
ant.top.missense.dgrp.ld$R2 <- (ant.top.missense.dgrp.ld$Coupling_R)^2
##Calculate residual r
ant.top.missense.dgrp.ld$Exp_Coupling_R <- exp(ant.top.missense.dgrp.ld$Coupling_R)
ant.top.missense.dgrp.ld$R_resid[!is.na(ant.top.missense.dgrp.ld$Coupling_R)] <- glm(data=ant.top.missense.dgrp.ld,Exp_Coupling_R~Distance)$residuals



## CONCORDANT INDEX

## Import LD info
conc.top.dgrp.ld <- read.table("conc.top.dgrp.ld",head=T)

#Obtain allelic effect of locus A ("SNP_A") from assoc info dataframe
tmp1 <- merge(conc.top.dgrp.ld,conc.top.assoc[c("Predictor","Effect","Same_minor_alleles","A1","A2","A1_DGRP","A2_DGRP")],by.x="SNP_A",by.y="Predictor",all.x=T)
#Obtain allelic effect of locus B ("SNP_B") from assoc info dataframe
conc.top.dgrp.ld <- merge(tmp1,conc.top.assoc,by.x="SNP_B",by.y="Predictor",all.x=T)
names(conc.top.dgrp.ld)[c(8:13,16:21)] <- c("Effect_A","Same_minor_alleles_A","A1_A","A2_A","A1_DGRP_A","A2_DGRP_A","Effect_B","A1_B","A2_B","A1_DGRP_B","A2_DGRP_B","Same_minor_alleles_B")
conc.top.dgrp.ld <- conc.top.dgrp.ld[order(conc.top.dgrp.ld$SNP_A,conc.top.dgrp.ld$SNP_B),]

## Create coupling LD column
#If allelic effects are in the same direction (both positive, or both negative), R is measuring coupling LD, so keep r2 as is; if the signs are opposite, R is measuring repulsion LD, so flip the sign of R
conc.top.dgrp.ld$Naive_coupling_R <- ifelse((conc.top.dgrp.ld$Effect_A>0 & conc.top.dgrp.ld$Effect_B>0) | (conc.top.dgrp.ld$Effect_A<0 & conc.top.dgrp.ld$Effect_B<0),conc.top.dgrp.ld$R,-conc.top.dgrp.ld$R)
#If minor alleles are the same in LHm and DGRP, then 'naive' coupling R value is correct; else the sign of R should be flipped
conc.top.dgrp.ld$Coupling_R <- with(conc.top.dgrp.ld,ifelse((Same_minor_alleles_A==1 & Same_minor_alleles_B==1) | (Same_minor_alleles_A==0 & Same_minor_alleles_B==0),R,ifelse((Same_minor_alleles_A==0 & Same_minor_alleles_B==1) | (Same_minor_alleles_A==1 & Same_minor_alleles_B==0),-R,NA)))
#Double check that plink calculates LD between major/minor alleles, not between reference/alternative!!  

## Create distance column
conc.top.dgrp.ld$Distance <- conc.top.dgrp.ld$BP_B-conc.top.dgrp.ld$BP_A

## Calculate r2
conc.top.dgrp.ld$R2 <- (conc.top.dgrp.ld$Coupling_R)^2

##Calculate residual r
conc.top.dgrp.ld$Exp_Coupling_R <- exp(conc.top.dgrp.ld$Coupling_R)
conc.top.dgrp.ld$R_resid[!is.na(conc.top.dgrp.ld$Coupling_R)] <- glm(data=conc.top.dgrp.ld,Exp_Coupling_R~Distance)$residuals

## Are SNPs situated in any common ivnersion
conc.top.dgrp.ld$is_common_inv_any <- ifelse(rowSums(conc.top.dgrp.ld[22:42])>0,1,0)


## CONCORDANT INDEX, MISSENSE SNPS

## Import LD info
conc.top.missense.dgrp.ld <- read.table("conc.top.missense.dgrp.ld",head=T)

#Obtain allelic effect of locus A ("SNP_A") from assoc info dataframe
tmp1 <- merge(conc.top.missense.dgrp.ld,conc.top.missense.assoc[c("Predictor","Effect","Same_minor_alleles","A1","A2","A1_DGRP","A2_DGRP")],by.x="SNP_A",by.y="Predictor",all.x=T)
#Obtain allelic effect of locus B ("SNP_B") from assoc info dataframe
conc.top.missense.dgrp.ld <- merge(tmp1,conc.top.missense.assoc,by.x="SNP_B",by.y="Predictor",all.x=T)
names(conc.top.missense.dgrp.ld)[c(8:13,16:21)] <- c("Effect_A","Same_minor_alleles_A","A1_A","A2_A","A1_DGRP_A","A2_DGRP_A","Effect_B","A1_B","A2_B","A1_DGRP_B","A2_DGRP_B","Same_minor_alleles_B")
conc.top.missense.dgrp.ld <- conc.top.missense.dgrp.ld[order(conc.top.missense.dgrp.ld$SNP_A,conc.top.missense.dgrp.ld$SNP_B),]

## Create coupling LD column
#If allelic effects are in the same direction (both positive, or both negative), R is measuring coupling LD, so keep r2 as is; if the signs are opposite, R is measuring repulsion LD, so flip the sign of R
conc.top.missense.dgrp.ld$Naive_coupling_R <- ifelse((conc.top.missense.dgrp.ld$Effect_A>0 & conc.top.missense.dgrp.ld$Effect_B>0) | (conc.top.missense.dgrp.ld$Effect_A<0 & conc.top.missense.dgrp.ld$Effect_B<0),conc.top.missense.dgrp.ld$R,-conc.top.missense.dgrp.ld$R)
#If minor alleles are the same in LHm and DGRP, then 'naive' coupling R value is correct; else the sign of R should be flipped
conc.top.missense.dgrp.ld$Coupling_R <- with(conc.top.missense.dgrp.ld,ifelse((Same_minor_alleles_A==1 & Same_minor_alleles_B==1) | (Same_minor_alleles_A==0 & Same_minor_alleles_B==0),R,ifelse((Same_minor_alleles_A==0 & Same_minor_alleles_B==1) | (Same_minor_alleles_A==1 & Same_minor_alleles_B==0),-R,NA)))
#Double check that plink calculates LD between major/minor alleles, not between reference/alternative!!  

## Create distance column
conc.top.missense.dgrp.ld$Distance <- conc.top.missense.dgrp.ld$BP_B-conc.top.missense.dgrp.ld$BP_A

## Calculate r2
conc.top.missense.dgrp.ld$R2 <- (conc.top.missense.dgrp.ld$Coupling_R)^2
##Calculate residual r
conc.top.missense.dgrp.ld$Exp_Coupling_R <- exp(conc.top.missense.dgrp.ld$Coupling_R)
conc.top.missense.dgrp.ld$R_resid[!is.na(conc.top.missense.dgrp.ld$Coupling_R)] <- glm(data=conc.top.missense.dgrp.ld,Exp_Coupling_R~Distance)$residuals



## PERMUTED INDEX 1

## Import LD info
perm.top.dgrp.ld <- read.table("perm.top.dgrp.ld",head=T)

#Obtain allelic effect of locus A ("SNP_A") from assoc info dataframe
tmp1 <- merge(perm.top.dgrp.ld,perm.top.assoc[c("Predictor","Effect","Same_minor_alleles","A1","A2","A1_DGRP","A2_DGRP")],by.x="SNP_A",by.y="Predictor",all.x=T)
#Obtain allelic effect of locus B ("SNP_B") from assoc info dataframe
perm.top.dgrp.ld <- merge(tmp1,perm.top.assoc,by.x="SNP_B",by.y="Predictor",all.x=T)
names(perm.top.dgrp.ld)[c(8:13,16:21)] <- c("Effect_A","Same_minor_alleles_A","A1_A","A2_A","A1_DGRP_A","A2_DGRP_A","Effect_B","A1_B","A2_B","A1_DGRP_B","A2_DGRP_B","Same_minor_alleles_B")
perm.top.dgrp.ld <- perm.top.dgrp.ld[order(perm.top.dgrp.ld$SNP_A,perm.top.dgrp.ld$SNP_B),]

## Create coupling LD column
#If allelic effects are in the same direction (both positive, or both negative), R is measuring coupling LD, so keep r2 as is; if the signs are opposite, R is measuring repulsion LD, so flip the sign of R
perm.top.dgrp.ld$Naive_coupling_R <- ifelse((perm.top.dgrp.ld$Effect_A>0 & perm.top.dgrp.ld$Effect_B>0) | (perm.top.dgrp.ld$Effect_A<0 & perm.top.dgrp.ld$Effect_B<0),perm.top.dgrp.ld$R,-perm.top.dgrp.ld$R)
#If minor alleles are the same in LHm and DGRP, then 'naive' coupling R value is correct; else the sign of R should be flipped
perm.top.dgrp.ld$Coupling_R <- with(perm.top.dgrp.ld,ifelse((Same_minor_alleles_A==1 & Same_minor_alleles_B==1) | (Same_minor_alleles_A==0 & Same_minor_alleles_B==0),R,ifelse((Same_minor_alleles_A==0 & Same_minor_alleles_B==1) | (Same_minor_alleles_A==1 & Same_minor_alleles_B==0),-R,NA)))
#Double check that plink calculates LD between major/minor alleles, not between reference/alternative!!  

## Create distance column
perm.top.dgrp.ld$Distance <- perm.top.dgrp.ld$BP_B-perm.top.dgrp.ld$BP_A

## Calculate r2
perm.top.dgrp.ld$R2 <- (perm.top.dgrp.ld$Coupling_R)^2

##Calculate residual r
perm.top.dgrp.ld$Exp_Coupling_R <- exp(perm.top.dgrp.ld$Coupling_R)
perm.top.dgrp.ld$R_resid[!is.na(perm.top.dgrp.ld$Coupling_R)] <- glm(data=perm.top.dgrp.ld,Exp_Coupling_R~Distance)$residuals

## Are SNPs situated in any common ivnersion
perm.top.dgrp.ld$is_common_inv_any <- ifelse(rowSums(perm.top.dgrp.ld[22:42])>0,1,0)



#################################################################
## Plot Coupling r for each inversion
#################################################################


#DGRP
#Plot ant and conc together
ant.top.dgrp.ld$Phenotype <- "ant"
conc.top.dgrp.ld$Phenotype <- "conc"
top.dgrp.ld <- rbind(conc.top.dgrp.ld,ant.top.dgrp.ld)

#Are SNPs situated in inversion that carries both antagonistic and concordant SNPs?
inv.contains.ant.and.conc <- names(which(apply(ant.top.dgrp.ld[22:42],2,sum)>0 & apply(conc.top.dgrp.ld[22:42],2,sum)>0))

library(ggplot2)
ggplot(subset(top.dgrp.ld,Distance<100000 & Coupling_R>0 & (is_common_inv_1==1 | is_common_inv_3==1 | is_common_inv_4==1 | is_common_inv_5==1 | is_common_inv_9==1 | is_common_inv_10==1 | is_common_inv_11==1 | is_common_inv_12==1 | is_common_inv_13==1 | is_common_inv_14==1 | is_common_inv_15==1)),aes(y=Coupling_R,x=Distance))+
  geom_smooth(method="lm",formula=y~log(x),aes(col=factor(Phenotype)),linetype = "dashed")+
  scale_color_manual(labels = c("SA", "non-SA"), values = c("purple", "orange"))+
  theme_bw()+
  ylab("Coupling LD (r)")+
  xlab("Distance (bp)")+
  ylim(c(0,1))+
  theme(axis.title = element_text(size=25),axis.text = element_text(size=25),legend.title=element_blank(),legend.text=element_text(size=25),legend.background = element_blank())


#Plot mean coupling LD for given bins of distance
top.dgrp.ld$Distance_bins_100k <- factor(cut(top.dgrp.ld$Distance,breaks = seq(1,100000,1000),labels=c(seq(1,100000,1000)[2:100])))
top.dgrp.ld$Distance_bins_10k <- factor(cut(top.dgrp.ld$Distance,breaks = seq(1,10000,100),labels=c(seq(1,10000,100)[2:100])))

#10k distance
top.dgrp.ld.bin <- aggregate(subset(top.dgrp.ld,Coupling_R>0 & (is_common_inv_1==1 | is_common_inv_3==1 | is_common_inv_4==1 | is_common_inv_5==1 | is_common_inv_9==1 | is_common_inv_10==1 | is_common_inv_11==1 | is_common_inv_12==1 | is_common_inv_13==1 | is_common_inv_14==1 | is_common_inv_15==1))$Coupling_R,by=list(subset(top.dgrp.ld,Coupling_R>0 & (is_common_inv_1==1 | is_common_inv_3==1 | is_common_inv_4==1 | is_common_inv_5==1 | is_common_inv_9==1 | is_common_inv_10==1 | is_common_inv_11==1 | is_common_inv_12==1 | is_common_inv_13==1 | is_common_inv_14==1 | is_common_inv_15==1))$Phenotype,subset(top.dgrp.ld,Coupling_R>0 & (is_common_inv_1==1 | is_common_inv_3==1 | is_common_inv_4==1 | is_common_inv_5==1 | is_common_inv_9==1 | is_common_inv_10==1 | is_common_inv_11==1 | is_common_inv_12==1 | is_common_inv_13==1 | is_common_inv_14==1 | is_common_inv_15==1))$Distance_bins_10k),FUN = function(x) c(mean = mean(x,na.rm=T)))
names(top.dgrp.ld.bin) <- c("Phenotype","Distance","Coupling_R")
top.dgrp.ld.bin$Distance <- as.numeric(as.character(top.dgrp.ld.bin$Distance))
  
ggplot(subset(top.dgrp.ld.bin),aes(y=Coupling_R,x=Distance))+
  geom_point(aes(col=factor(Phenotype)))+
  geom_line(aes(col=factor(Phenotype),group=factor(Phenotype)))+
  scale_color_manual(labels = c("SA", "control"), values = c("darkorchid4", "pink1"))+
  theme_bw()+
  ylab(expression(paste("Coupling LD (", r[w], ", DGRP)",sep="")))+
  xlab("Distance (bp)")+
  ylim(c(0,1))+
  scale_x_continuous(breaks=c(0,5000,10000))+
  theme(axis.title = element_text(size=30),plot.margin=unit(c(1.5,1.5,1.5,1.5),"cm"),legend.position = c(0.7, 0.8),axis.text = element_text(size=25),legend.title=element_blank(),legend.text=element_text(size=30),legend.background = element_blank())+
  guides(color = guide_legend(override.aes = list(size=5)))

#100k distance
top.dgrp.ld.bin.100k <- aggregate(subset(top.dgrp.ld,Coupling_R>0 & (is_common_inv_1==1 | is_common_inv_3==1 | is_common_inv_4==1 | is_common_inv_5==1 | is_common_inv_9==1 | is_common_inv_10==1 | is_common_inv_11==1 | is_common_inv_12==1 | is_common_inv_13==1 | is_common_inv_14==1 | is_common_inv_15==1))$Coupling_R,by=list(subset(top.dgrp.ld,Coupling_R>0 & (is_common_inv_1==1 | is_common_inv_3==1 | is_common_inv_4==1 | is_common_inv_5==1 | is_common_inv_9==1 | is_common_inv_10==1 | is_common_inv_11==1 | is_common_inv_12==1 | is_common_inv_13==1 | is_common_inv_14==1 | is_common_inv_15==1))$Phenotype,subset(top.dgrp.ld,Coupling_R>0 & (is_common_inv_1==1 | is_common_inv_3==1 | is_common_inv_4==1 | is_common_inv_5==1 | is_common_inv_9==1 | is_common_inv_10==1 | is_common_inv_11==1 | is_common_inv_12==1 | is_common_inv_13==1 | is_common_inv_14==1 | is_common_inv_15==1))$Distance_bins_100k),FUN = function(x) c(mean = mean(x,na.rm=T)))
names(top.dgrp.ld.bin.100k) <- c("Phenotype","Distance","Coupling_R")
top.dgrp.ld.bin.100k$Distance <- as.numeric(as.character(top.dgrp.ld.bin.100k$Distance))

ggplot(subset(top.dgrp.ld.bin.100k),aes(y=Coupling_R,x=Distance))+
  geom_point(aes(col=factor(Phenotype)))+
  geom_line(aes(col=factor(Phenotype),group=factor(Phenotype)))+
  scale_color_manual(labels = c("SA", "non-SA"), values = c("darkorchid4", "pink1"))+
  theme_bw()+
  ylab(expression(paste("Coupling LD (", r[w], ", DGRP)",sep="")))+
  xlab("Distance (bp)")+
  ylim(c(0,1))+
  theme(axis.title = element_text(size=30),axis.text = element_text(size=25),legend.title=element_blank(),legend.text=element_text(size=40),legend.background = element_blank())+
  guides(color = guide_legend(override.aes = list(size=5)))


#10k distance, for a promising inversion
top.dgrp.ld.bin <- aggregate(subset(top.dgrp.ld,Coupling_R>0 & is_common_inv_14==1)$Coupling_R,by=list(subset(top.dgrp.ld,Coupling_R>0 & is_common_inv_14==1)$Phenotype,subset(top.dgrp.ld,Coupling_R>0 & is_common_inv_14==1)$Distance_bins_10k),FUN = function(x) c(mean = mean(x,na.rm=T)))
names(top.dgrp.ld.bin) <- c("Phenotype","Distance","Coupling_R")
top.dgrp.ld.bin$Distance <- as.numeric(as.character(top.dgrp.ld.bin$Distance))

ggplot(subset(top.dgrp.ld.bin),aes(y=Coupling_R,x=Distance))+
  geom_point(aes(col=factor(Phenotype)))+
  geom_line(aes(col=factor(Phenotype),group=factor(Phenotype)))+
  scale_color_manual(labels = c("SA", "control"), values = c("darkorchid4", "pink1"))+
  theme_bw()+
  ylab(expression(paste("Coupling LD (", r[w], "), DGRP population.",sep="")))+
  xlab("Distance (bp)")+
  ylim(c(0,1))+
  scale_x_continuous(breaks=c(0,5000,10000))+
  theme(axis.title = element_text(size=30),plot.margin=unit(c(1.5,1.5,1.5,1.5),"cm"),legend.position = c(0.7, 0.8),axis.text = element_text(size=25),legend.title=element_blank(),legend.text=element_text(size=30),legend.background = element_blank())+
  guides(color = guide_legend(override.aes = list(size=5)))



#Promising inversions
ggplot(subset(top.dgrp.ld,Distance<10000 & Coupling_R>0),aes(y=Coupling_R,x=Distance))+
  geom_point(size=0.1)+
  geom_smooth(method="lm",formula=y~log(x),aes(col=factor(Phenotype)))+
  theme_bw()

#Promising inversions
ggplot(subset(top.dgrp.ld,Distance<10000 & is_common_inv_9==1 & Coupling_R>0),aes(y=Coupling_R,x=Distance))+
  geom_point(size=0.1,aes(col=factor(Phenotype)))+
  #geom_point(size=0.1)+
  geom_smooth(method="lm",formula=y~log(x),aes(col=factor(Phenotype)))+
  theme_bw()

ggplot(subset(top.dgrp.ld,Distance<10000 & is_common_inv_14==1 & Coupling_R>0),aes(y=Coupling_R,x=Distance))+
  #geom_point(size=0.1,aes(col=factor(Phenotype)))+
  geom_point(size=0.1)+
  geom_smooth(method="lm",formula=y~log(x),aes(col=factor(Phenotype)))+
  theme_bw()

ggplot(subset(top.dgrp.ld,Distance<100000 & is_common_inv_12==1 & Coupling_R>0),aes(y=Coupling_R,x=Distance))+
  #geom_point(size=0.1,aes(col=factor(Phenotype)))+
  geom_point(size=0.1)+
  geom_smooth(method="lm",formula=y~log(x),aes(col=factor(Phenotype)))+
  theme_bw()


#DGRP, missense
#Plot ant and conc together
ant.top.missense.dgrp.ld$Phenotype <- "ant"
conc.top.missense.dgrp.ld$Phenotype <- "conc"
top.missense.dgrp.ld <- rbind(conc.top.missense.dgrp.ld,ant.top.missense.dgrp.ld)

library(ggplot2)

#Promising inversions
ggplot(subset(top.missense.dgrp.ld,Distance<10000 & Coupling_R>0),aes(y=Coupling_R,x=Distance))+
  geom_point(size=0.1,aes(col=factor(Phenotype)))+
  geom_smooth(method="lm",formula=y~log(x),aes(col=factor(Phenotype)))+
  theme_bw()

#Promising inversions
ggplot(subset(top.missense.dgrp.ld,Distance<10000 & is_common_inv_9==1 & Coupling_R>0),aes(y=Coupling_R,x=Distance))+
  geom_point(size=0.1,aes(col=factor(Phenotype)))+
  #geom_point(size=0.1)+
  geom_smooth(method="lm",formula=y~log(x),aes(col=factor(Phenotype)))+
  theme_bw()

ggplot(subset(top.missense.dgrp.ld,Distance<1000000 & is_common_inv_14==1 & Coupling_R>0),aes(y=Coupling_R,x=Distance))+
  geom_point(size=0.1,aes(col=factor(Phenotype)))+
  geom_smooth(method="lm",formula=y~log(x),aes(col=factor(Phenotype)))+
  theme_bw()



#################################################################
## Statistical test: coupling r ~ ant vs conc
#################################################################


