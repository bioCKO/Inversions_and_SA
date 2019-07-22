## Converting breakpoints to genomic coordinates

#Set working directory
setwd("~/Documents/inversions_and_SA/")

#Import conversion table
ref <- read.csv("DmelMapTable.160615c.csv")

#Import inversion data table
dd <- read.csv("test.csv")

#Add genomic coordinate information to inversion data
dd2 <- merge(dd,ref[,c("Cytogenetic.map.position","Sequence.coordinates..release.6.")],all.x=T,by.x="Start",by.y="Cytogenetic.map.position")
dd3 <- merge(dd2,ref[,c("Cytogenetic.map.position","Sequence.coordinates..release.6.")],all.x=T,by.x="End",by.y="Cytogenetic.map.position")

#Re-order by cytological starting position
dd3 <- dd3[order(dd3$Start),]

#Clean up genomic coordinate columns
#Breakpoint start
dd3$Start_part1 <- ifelse(!is.na(dd3$Sequence.coordinates..release.6..x) , do.call(rbind,strsplit(as.character(dd3$Sequence.coordinates..release.6..x),split="..", fixed = TRUE))[,1] , NA)
dd3$Start_part2 <- ifelse(!is.na(dd3$Sequence.coordinates..release.6..x),do.call(rbind,strsplit(as.character(dd3$Sequence.coordinates..release.6..x),split="..", fixed = TRUE))[,2],NA)
dd3$Start_Chrom <- ifelse(!is.na(dd3$Start_part1),do.call(rbind,strsplit(as.character(dd3$Start_part1),split=":", fixed = TRUE))[,1],NA)
dd3$Start_part1 <- ifelse(!is.na(dd3$Start_part1),do.call(rbind,strsplit(as.character(dd3$Start_part1),split=":", fixed = TRUE))[,2],NA)
dd3$Start_midpoint <- round((as.numeric(dd3$Start_part1)+as.numeric(dd3$Start_part2))/2)

#Breakpoint end
dd3$End_part1 <- ifelse(!is.na(dd3$Sequence.coordinates..release.6..y),do.call(rbind,strsplit(as.character(dd3$Sequence.coordinates..release.6..y),split="..", fixed = TRUE))[,1],NA)
dd3$End_part2 <- ifelse(!is.na(dd3$Sequence.coordinates..release.6..y),do.call(rbind,strsplit(as.character(dd3$Sequence.coordinates..release.6..y),split="..", fixed = TRUE))[,2],NA)
dd3$End_Chrom <- ifelse(!is.na(dd3$End_part1),do.call(rbind,strsplit(as.character(dd3$End_part1),split=":", fixed = TRUE))[,1],NA)
dd3$End_part1 <- ifelse(!is.na(dd3$End_part1),do.call(rbind,strsplit(as.character(dd3$End_part1),split=":", fixed = TRUE))[,2],NA)
dd3$End_midpoint <- round((as.numeric(dd3$End_part1)+as.numeric(dd3$End_part2))/2)

