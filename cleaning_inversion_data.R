## Converting breakpoints to genomic coordinates
#Set working directory
setwd("~/Documents/inversions_and_SA/")

#Import conversion table
ref <- read.csv("DmelMapTable.160615c_v2.csv")

#Import inversion data table
dd <- read.csv("inversion_raw_table.csv")
#Remove some useless columns
dd <- dd[1:19]

#Add genomic coordinate information to inversion data
dd2 <- merge(dd,ref[,c("Cytogenetic.map.position","Sequence.coordinates..release.6.")],all.x=T,by.x="Start",by.y="Cytogenetic.map.position")
dd3 <- merge(dd2,ref[,c("Cytogenetic.map.position","Sequence.coordinates..release.6.")],all.x=T,by.x="End",by.y="Cytogenetic.map.position")

## Clean up genomic coordinate columns
#Breakpoint start
for (i in 1:nrow(dd3)){
dd3$Start_part1[i] <- strsplit(as.character(dd3$Sequence.coordinates..release.6..x[i]),split="..", fixed = TRUE)[[1]][1]
dd3$Start_part2[i] <- strsplit(as.character(dd3$Sequence.coordinates..release.6..x[i]),split="..", fixed = TRUE)[[1]][2]
dd3$Start_Chrom[i] <- strsplit(as.character(dd3$Start_part1[i]),split=":", fixed = TRUE)[[1]][1]
dd3$Start_part1[i] <- strsplit(as.character(dd3$Start_part1[i]),split=":", fixed = TRUE)[[1]][2]
}
dd3$Start_converted_midpoint <- round((as.numeric(dd3$Start_part1)+as.numeric(dd3$Start_part2))/2)

#Breakpoint end
for (i in 1:nrow(dd3)){
  dd3$End_part1[i] <- strsplit(as.character(dd3$Sequence.coordinates..release.6..y[i]),split="..", fixed = TRUE)[[1]][1]
  dd3$End_part2[i] <- strsplit(as.character(dd3$Sequence.coordinates..release.6..y[i]),split="..", fixed = TRUE)[[1]][2]
  dd3$End_Chrom[i] <- strsplit(as.character(dd3$End_part1[i]),split=":", fixed = TRUE)[[1]][1]
  dd3$End_part1[i] <- strsplit(as.character(dd3$End_part1[i]),split=":", fixed = TRUE)[[1]][2]
}
dd3$End_converted_midpoint <- round((as.numeric(dd3$End_part1)+as.numeric(dd3$End_part2))/2)


#Clean columns with missing start/end values
#Should be 28 of them, i.e. 27 Chakraborty inversions + single additional inversion from RCD 2012
#Manually make modifications to input files when this does not match up (e.g. adding 98F band in DmelMapTable.160615c_v2.csv file, because it was not there in the first place)
nrow(subset(dd3,is.na(End_Chrom)))
nrow(subset(dd3,is.na(Start_Chrom)))

#Remove extraneous columns
inv <- dd3[,c(1:19,24:25,28:29)]

#Write table
write.table(inv,"inversions_plus_genomic_coordinates.txt",row.names=F)


