## Combine gene coordinate data with antagonistic and concordant SNPs from GWAS 

setwd("~/Documents/inversions_and_SA/")

############################################################
## Import table of antagonistic and concordant SNPs from https://zenodo.org/record/2623225
############################################################

assoc <- read.table("Assoc_Data.txt",h=T)

## Create antagonistic and concordant tables
#Define top 2,372 concordant SNPs as equivalent to 2,372 antagonistic SNPs
assoc$Cand_concordant <- ifelse(assoc$Wald_P_concordant<=sort(assoc$Wald_P_concordant)[2372],1,0)
#Keep only relevant columns
assoc.ant <- subset(assoc,Cand_parametric==1)[c("Chromosome","Basepair")]
assoc.conc <- subset(assoc,Cand_concordant==1)[c("Chromosome","Basepair")]
#Keep only candidate genes with missense SNPs
assoc.ant.missense <- subset(assoc,Cand_parametric==1 & Consequence=="missense_variant")[c("Chromosome","Basepair")]
assoc.conc.missense <- subset(assoc,Cand_concordant==1 & Consequence=="missense_variant")[c("Chromosome","Basepair")]

write.table(assoc.ant,"antagonistic_SNPs.txt",quote=F,sep="\t",row.names=F)
write.table(assoc.conc,"concordant_SNPs.txt",quote=F,sep="\t",row.names=F)

############################################################
## Import genes
############################################################

genes <- read.table("ensgenes.txt")
names(genes) <- c("Gene_name","Chromosome_numeric","Start","End")
genes$Chromosome <- factor(genes$Chromosome_numeric)
levels(genes$Chromosome) <- c("2L","2R","3L","3R","X")
#Extend gene coordinates
genes$Start_minus_5k <- genes$Start-5000
genes$End_plus_5k <- genes$End+5000

############################################################
## Import inversions
############################################################

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

############################################################
## Does a gene overlap an antagonistic/concordant SNP?
############################################################

#If the antagonistic SNP (i) falls within a given gene (j)'s coordinates, the antagonistic column is 1
#If the concordant SNP (k) falls within a given gene (j)'s coordinates, the concordant column is 1
#To speed things up, ignore this procedure for any genes where current value of the antagonistic column is already 1, or where the antagonistic SNP and genes are on different chromosomes (and do the same for concordant genes)
genes$is_ant <- 0
genes$is_conc <- 0
genes$is_ant_missense <- 0
genes$is_conc_missense <- 0

for (j in 1:nrow(genes)){
  for (i in 1:nrow(assoc.ant)){
     if(genes$is_ant[j]==0 & assoc.ant$Chromosome[i]==genes$Chromosome[j]){
        genes$is_ant[j] <- ifelse(assoc.ant$Basepair[i]>=genes$Start_minus_5k[j] & assoc.ant$Basepair[i]<=genes$End_plus_5k[j],1,genes$is_ant[j])
     }
  }
  for (k in 1:nrow(assoc.conc)){
      if(genes$is_conc[j]==0 & assoc.conc$Chromosome[k]==genes$Chromosome[j]){
        genes$is_conc[j] <- ifelse(assoc.conc$Basepair[k]>=genes$Start_minus_5k[j] & assoc.conc$Basepair[k]<=genes$End_plus_5k[j],1,genes$is_conc[j])
      }
  }
  #for (i in 1:nrow(assoc.ant.missense)){
   # if(genes$is_ant_missense[j]==0 & assoc.ant.missense$Chromosome[i]==genes$Chromosome[j]){
   #   genes$is_ant_missense[j] <- ifelse(assoc.ant.missense$Basepair[i]>=genes$Start_minus_5k[j] & assoc.ant.missense$Basepair[i]<=genes$End_plus_5k[j],1,genes$is_ant_missense[j])
   # }
 # }
  #for (i in 1:nrow(assoc.conc.missense)){
   # if(genes$is_conc_missense[j]==0 & assoc.conc.missense$Chromosome[i]==genes$Chromosome[j]){
    #  genes$is_conc_missense[j] <- ifelse(assoc.conc.missense$Basepair[i]>=genes$Start_minus_5k[j] & assoc.conc.missense$Basepair[i]<=genes$End_plus_5k[j],1,genes$is_conc_missense[j])
   # }
  #}
  print(j)
}

write.table(genes,"antagonistic_and_concordant_genes.txt",quote=F,sep="\t",row.names=F)

genes <- read.table("antagonistic_and_concordant_genes.txt",h=T)

#How many genes are antagonistic/concordant
sum(genes$is.ant)
sum(genes$is.conc)

#Subsets of genes which only include antagonistic and concordant genes
genes.ant <- subset(genes,is.ant==1)
genes.conc <- subset(genes,is.conc==1)

############################################################
## How does fraction of genome covered fraction of the genome is covered by each inversion? 
############################################################

## Fraction of genome covered by each inversion
#Total genome content = 23513712+25286936+28110227+32079331+23542271=132532477
inv$Fraction_of_genome <- inv$Length/132532477

## Expected number of antagonistic genes within a given inversion
inv$Expected_antagonistic_genes <- inv$Fraction_of_genome*sum(genes$is.ant==1)
inv$Expected_concordant_genes <- inv$Fraction_of_genome*sum(genes$is.conc==1)
inv$Expected_antagonistic_genes_sim1 <- rpois(nrow(inv),inv$Expected_antagonistic_genes)
inv$Expected_concordant_genes_sim1 <- rpois(nrow(inv),inv$Expected_concordant_genes)

## Observed number of antagonistic genes within a given inversion
inv$Observed_antagonistic_genes <- 0
inv$Observed_concordant_genes <- 0

inv <- inv[order(inv$Chromosome,inv$Start_definitive),]

#Record number of antagonistic genes that overlap each inversion
for (k in 1:nrow(inv)){
  #Antagonistic genes
  for (j in 1:nrow(genes.ant)){
    if(inv$Chromosome[k]==genes.ant$Chromosome[j]){
      inv$Observed_antagonistic_genes[k] <- ifelse((inv$Start_definitive[k]>=genes.ant$Start_minus_5k[j] & inv$Start_definitive[k]<=genes.ant$End_plus_5k[j])|(inv$End_definitive[k]>=genes.ant$Start_minus_5k[j] & inv$End_definitive[k]<=genes.ant$End_plus_5k[j])|(inv$Start_definitive[k]<=genes.ant$Start_minus_5k[j] & inv$End_definitive[k]>=genes.ant$End_plus_5k[j]),1+inv$Observed_antagonistic_genes[k],inv$Observed_antagonistic_genes[k])
    }
  }
  #Concordant genes
  for (j in 1:nrow(genes.conc)){
    if(inv$Chromosome[k]==genes.conc$Chromosome[j]){
      inv$Observed_concordant_genes[k] <- ifelse((inv$Start_definitive[k]>=genes.conc$Start_minus_5k[j] & inv$Start_definitive[k]<=genes.conc$End_plus_5k[j])|(inv$End_definitive[k]>=genes.conc$Start_minus_5k[j] & inv$End_definitive[k]<=genes.conc$End_plus_5k[j])|(inv$Start_definitive[k]<=genes.conc$Start_minus_5k[j] & inv$End_definitive[k]>=genes.conc$End_plus_5k[j]),1+inv$Observed_concordant_genes[k],inv$Observed_concordant_genes[k])
    }
  }
  print(k)
}

## Ratio of observed to expected for each inversion
inv$Ratio_antagonistic_observed_to_expected <- log2(inv$Observed_antagonistic_genes/inv$Expected_antagonistic_genes)
inv$Ratio_concordant_observed_to_expected <- log2(inv$Observed_concordant_genes/inv$Expected_concordant_genes)

inv$Ratio_antagonistic_observed_to_expected_sim1 <- log2(inv$Observed_antagonistic_genes/inv$Expected_antagonistic_genes_sim1)
inv$Ratio_concordant_observed_to_expected_sim1 <- log2(inv$Observed_concordant_genes/inv$Expected_concordant_genes_sim1)


##Histogram of observed and expected antagonistic/concordant genes in each inversion; considering all inversions
library(ggplot2)
library(reshape2)

inv_plot <- melt(inv[c(1,10:15)])
inv_plot2 <- subset(inv_plot,variable!="Expected_concordant_genes" & variable!="Expected_antagonistic_genes"  & variable!="Expected_concordant_genes_sim1")
levels(inv_plot2$variable)
inv_plot2$variable <- factor(inv_plot2$variable, levels = c("Expected_antagonistic_genes_sim1","Observed_concordant_genes","Observed_antagonistic_genes"))

p1 <- ggplot(inv_plot2,aes(value,fill=variable))+
  geom_histogram(alpha=0.75,position="identity")+
  scale_fill_manual(values=c("grey75","pink1","darkorchid4"),labels = c("control1","control2","SA"))+
  theme_bw()+
  theme(axis.title = element_text(size=30),plot.margin=unit(c(1.5,1.5,1.5,1.5),"cm"),legend.position = c(0.7, 0.8),axis.text = element_text(size=25),legend.title=element_blank(),legend.text=element_text(size=30),legend.background = element_blank())+
  xlab("Genes per inversion")+
  ylab("Count")
p1

p2 <- ggplot(subset(inv_plot,variable=="Ratio_antagonistic_observed_to_expected_sim1" | variable=="Ratio_concordant_observed_to_expected_sim1"),aes(value,fill=variable))+
  geom_histogram(alpha=0.5,position="identity")+
  scale_fill_manual(values=c("red","blue"))+
  theme_bw()+
  xlab("log2 ratio (obs/exp)")+
  xlim(c(-6,6))+
  geom_vline(xintercept=0)

library(gridExtra)
grid.arrange(p1,p2)


##Histogram of observed and expected antagonistic/concordant genes in each inversion; considering common inversions
inv_plot <- melt(inv[inv$Common==1,c(1,10:15)])
inv_plot2 <- subset(inv_plot,variable!="Expected_concordant_genes" & variable!="Expected_antagonistic_genes"  & variable!="Expected_concordant_genes_sim1")
levels(inv_plot2$variable)
inv_plot2$variable <- factor(inv_plot2$variable, levels = c("Expected_antagonistic_genes_sim1","Observed_concordant_genes","Observed_antagonistic_genes"))

p1 <- ggplot(inv_plot2,aes(value,fill=variable))+
  geom_histogram(alpha=0.75,position="identity")+
  scale_fill_manual(values=c("grey75","pink1","darkorchid4"),labels = c("control1","control2","SA"))+
  theme_bw()+
  theme(axis.title = element_text(size=30),plot.margin=unit(c(1.5,1.5,1.5,1.5),"cm"),legend.position = c(0.7, 0.8),axis.text = element_text(size=25),legend.title=element_blank(),legend.text=element_text(size=30),legend.background = element_blank())+
  xlab("Genes per common inversion")+
  ylab("Count")
p1


p2 <- ggplot(subset(inv_plot,variable=="Ratio_antagonistic_observed_to_expected" | variable=="Ratio_concordant_observed_to_expected"),aes(value,fill=variable))+
  geom_histogram(alpha=0.5,position="identity")+
  scale_fill_manual(values=c("red","blue"))+
  theme_bw()+
  xlab("log2 ratio (obs/exp)")+
  xlim(c(-6,6))+
  geom_vline(xintercept=0)

library(gridExtra)
grid.arrange(p1,p2)

##Histogram of observed and expected antagonistic/concordant genes in each inversion; considering common inversions
inv_plot <- melt(inv[inv$Rare==1,c(1,10:15)])
inv_plot2 <- subset(inv_plot,variable!="Expected_concordant_genes" & variable!="Expected_antagonistic_genes"  & variable!="Expected_concordant_genes_sim1")
levels(inv_plot2$variable)
inv_plot2$variable <- factor(inv_plot2$variable, levels = c("Expected_antagonistic_genes_sim1","Observed_concordant_genes","Observed_antagonistic_genes"))

p1 <- ggplot(inv_plot2,aes(value,fill=variable))+
  geom_histogram(alpha=0.75,position="identity")+
  scale_fill_manual(values=c("grey75","pink1","darkorchid4"),labels = c("control1","control2","SA"))+
  theme_bw()+
  theme(axis.title = element_text(size=30),plot.margin=unit(c(1.5,1.5,1.5,1.5),"cm"),legend.position = c(0.7, 0.8),axis.text = element_text(size=25),legend.title=element_blank(),legend.text=element_text(size=30),legend.background = element_blank())+
  xlab("Genes per rare inversion")+
  ylab("Count")
p1


p2 <- ggplot(subset(inv_plot,variable=="Ratio_antagonistic_observed_to_expected" | variable=="Ratio_concordant_observed_to_expected"),aes(value,fill=variable))+
  geom_histogram(alpha=0.5,position="identity")+
  scale_fill_manual(values=c("red","blue"))+
  theme_bw()+
  xlab("log2 ratio (obs/exp)")+
  xlim(c(-6,6))+
  geom_vline(xintercept=0)

library(gridExtra)
grid.arrange(p1,p2)


############################################################
## Do genes overlap with inversions? Simple analysis
############################################################

## If the inversion (k=common; l=rare) falls within a gene (j), the overlap column is 1
genes$is_common_inv <- 0
genes$is_rare_inv <- 0

for (j in 1:nrow(genes)){
  #Common inversions
  for (k in 1:nrow(inv.common)){
    if(genes$is_common_inv[j]==0 & inv.common$Chromosome[k]==genes$Chromosome[j]){
      genes$is_common_inv[j] <- ifelse((inv.common$Start_definitive[k]>=genes$Start_minus_5k[j] & inv.common$Start_definitive[k]<=genes$End_plus_5k[j])|(inv.common$End_definitive[k]>=genes$Start_minus_5k[j] & inv.common$End_definitive[k]<=genes$End_plus_5k[j])|(inv.common$Start_definitive[k]<=genes$Start_minus_5k[j] & inv.common$End_definitive[k]>=genes$End_plus_5k[j]),1,genes$is_common_inv[j])
    }
  }
  #Rare inversions
  for (l in 1:nrow(inv.rare)){
    if(genes$is_rare_inv[j]==0 & inv.rare$Chromosome[l]==genes$Chromosome[j]){
      genes$is_rare_inv[j] <- ifelse((inv.rare$Start_definitive[l]>=genes$Start_minus_5k[j] & inv.rare$Start_definitive[l]<=genes$End_plus_5k[j])|(inv.rare$End_definitive[l]>=genes$Start_minus_5k[j] & inv.rare$End_definitive[l]<=genes$End_plus_5k[j])|(inv.rare$Start_definitive[l]<=genes$Start_minus_5k[j] & inv.rare$End_definitive[l]>=genes$End_plus_5k[j]),1,genes$is_rare_inv[j])
    }
  }
  print(j)
}

#All inversions
genes$is_inv <- ifelse(genes$is_rare_inv==1 | genes$is_common_inv==1,1,0)

#How many antagonistic genes overlap with common inversions?
sum(genes$is.ant & genes$is_common_inv)#224
sum(genes$is.conc & genes$is_common_inv)#558
#How many antagonistic genes overlap with rare inversions? (this really needs to be done using a sample of rare inversions that matches the distribution and/or length of common inversions, see below)
sum(genes$is.ant & genes$is_rare_inv)#510
sum(genes$is.conc & genes$is_rare_inv)#671

############################################################
## Do genes overlap with inversions? Analysis matched for distribution of common inversion across chromosome arms
############################################################

## Sample rare inversions to match the distribution of common inversions

table(inv.common$Common,inv.common$Chromosome)
#Number of common inversions on each chromosome: 2L(5),2R(3),3L(5),3R(5),X(3)

out1 <- matrix(nrow=1000,ncol=2)
set.seed(123)
for (i in 1:1000){
  #Make empty vector of rare inversions
  genes.ant$is_rare_inv_matched <- 0
  genes.conc$is_rare_inv_matched <- 0
  #Sample matched set of rare inversions
  inv.rare.matched.rownames <- c(sample(rownames(subset(inv.rare,Chromosome=="2L")),size =  5),sample(rownames(subset(inv.rare,Chromosome=="2R")),size = 3),sample(rownames(subset(inv.rare,Chromosome=="3L")),size = 5),sample(rownames(subset(inv.rare,Chromosome=="3R")),size = 5),sample(rownames(subset(inv.rare,Chromosome=="X")),size = 3))
  inv.rare.matched <- inv.rare[inv.rare.matched.rownames,]
  #Check if antagonistic genes overlap with set of rare inversions
  for (l in 1:nrow(inv.rare.matched)){
    for (j in 1:nrow(genes.ant)){
      if(inv.rare.matched$Chromosome[l]==genes.ant$Chromosome[j]){
        genes.ant$is_rare_inv_matched[j] <- ifelse((inv.rare.matched$Start_definitive[l]>=genes.ant$Start_minus_5k[j] & inv.rare.matched$Start_definitive[l]<=genes.ant$End_plus_5k[j])|(inv.rare.matched$End_definitive[l]>=genes.ant$Start_minus_5k[j] & inv.rare.matched$End_definitive[l]<=genes.ant$End_plus_5k[j])|(inv.rare.matched$Start_definitive[l]<=genes.ant$Start_minus_5k[j] & inv.rare.matched$End_definitive[l]>=genes.ant$End_plus_5k[j]),1,genes.ant$is_rare_inv_matched[j])
      }
    }
    #Check if concordant genes overlap with the set of rare inversions
    for (k in 1:nrow(genes.conc)){
      if(inv.rare.matched$Chromosome[l]==genes.conc$Chromosome[k]){
        genes.conc$is_rare_inv_matched[k] <- ifelse((inv.rare.matched$Start_definitive[l]>=genes.conc$Start_minus_5k[k] & inv.rare.matched$Start_definitive[l]<=genes.conc$End_plus_5k[k])|(inv.rare.matched$End_definitive[l]>=genes.conc$Start_minus_5k[k] & inv.rare.matched$End_definitive[l]<=genes.conc$End_plus_5k[k])|(inv.rare.matched$Start_definitive[l]<=genes.conc$Start_minus_5k[k] & inv.rare.matched$End_definitive[l]>=genes.conc$End_plus_5k[k]),1,genes.conc$is_rare_inv_matched[k])
      }
    }
    out1[i,1] <- sum(genes.ant$is_rare_inv_matched)
    out1[i,2] <- sum(genes.conc$is_rare_inv_matched)
  }
  print(i)
}


############################################################
## Do genes overlap with inversions? Analysis matched for distribution and length of common inversion across chromosome arms
############################################################

## Sample inversions breakpoints to match the distribution of common inversions
#Number of common inversions on each chromosome: 2L(5),2R(3),3L(5),3R(5),X(3)

out2 <- matrix(nrow=1000,ncol=2)
set.seed(123)
for (i in 1:1000){
  #Make empty vector of rare inversions
  genes.ant$is_rare_inv_matched2 <- 0
  genes.conc$is_rare_inv_matched2 <- 0
  #Sample set of rare inversions that are matched by i) chromosomal distribution, and ii) common inversion length
  inv.rare.matched.rownames <- c(sample(rownames(subset(inv.rare,Chromosome=="2L")),size =  5),sample(rownames(subset(inv.rare,Chromosome=="2R")),size = 3),sample(rownames(subset(inv.rare,Chromosome=="3L")),size = 5),sample(rownames(subset(inv.rare,Chromosome=="3R")),size = 5),sample(rownames(subset(inv.rare,Chromosome=="X")),size = 3))
  inv.rare.matched <- inv.rare[inv.rare.matched.rownames,]
  #Choose a made-up end breakpoint
  #First make one possible end breakpoint at the proximal end of the chromosome
  inv.rare.matched$End_madeup_plus <- inv.rare.matched$Start_definitive+inv.common$Length
  #Then make one possible end breakpoint at the distal end of the chromosome
  inv.rare.matched$End_madeup_minus <- inv.rare.matched$Start_definitive-inv.common$Length
  #If one of them does not fall within chromosome coordinates, do not choose it; if both do, choose either one at random
  inv.rare.matched$End_madeup_chosen <- with(inv.rare.matched,ifelse(End_madeup_plus>inv.common$Max & End_madeup_minus>0,End_madeup_minus,ifelse(End_madeup_plus<inv.common$Max & End_madeup_minus<1,End_madeup_plus,ifelse(End_madeup_plus<inv.common$Max & End_madeup_minus>0,apply(inv.rare.matched[8:9], 1, sample, size = 1),End_madeup_plus))))
  #Check if antagonistic genes overlap with set of rare inversions
  for (l in 1:nrow(inv.rare.matched)){
  for (j in 1:nrow(genes.ant)){
      if(inv.rare.matched$Chromosome[l]==genes.ant$Chromosome[j]){
        genes.ant$is_rare_inv_matched2[j] <- ifelse((inv.rare.matched$Start_definitive[l]>=genes.ant$Start_minus_5k[j] & inv.rare.matched$Start_definitive[l]<=genes.ant$End_plus_5k[j])|(inv.rare.matched$End_madeup_chosen[l]>=genes.ant$Start_minus_5k[j] & inv.rare.matched$End_madeup_chosen[l]<=genes.ant$End_plus_5k[j])|(inv.rare.matched$Start_definitive[l]<=genes.ant$Start_minus_5k[j] & inv.rare.matched$End_madeup_chosen[l]>=genes.ant$End_plus_5k[j]),1,genes.ant$is_rare_inv_matched2[j])
      }
    }
  #Check if concordant genes overlap with the set of rare inversions
  for (k in 1:nrow(genes.conc)){
      if(inv.rare.matched$Chromosome[l]==genes.conc$Chromosome[k]){
        genes.conc$is_rare_inv_matched2[k] <- ifelse((inv.rare.matched$Start_definitive[l]>=genes.conc$Start_minus_5k[k] & inv.rare.matched$Start_definitive[l]<=genes.conc$End_plus_5k[k])|(inv.rare.matched$End_madeup_chosen[l]>=genes.conc$Start_minus_5k[k] & inv.rare.matched$End_madeup_chosen[l]<=genes.conc$End_plus_5k[k])|(inv.rare.matched$Start_definitive[l]<=genes.conc$Start_minus_5k[k] & inv.rare.matched$End_madeup_chosen[l]>=genes.conc$End_plus_5k[k]),1,genes.conc$is_rare_inv_matched2[k])
      }
  }
  out2[i,1] <- sum(genes.ant$is_rare_inv_matched2)
  out2[i,2] <- sum(genes.conc$is_rare_inv_matched2)
  }
  print(i)
}




