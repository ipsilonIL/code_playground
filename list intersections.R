##### rearrange column positions in a csv file
rm(list=ls())
setwd("/Users/piazza/Desktop/")

library("stringr", lib.loc="\\\\nas-green-1/share-green-appv-1-$/r_library")
library("stringr")

##########################################################################################GetIN
peptidelist_A <- read.csv('testA.csv') 
peptidelist_B <- read.csv('testB.csv') 
head(peptidelist_A)
head(peptidelist_B)

##########################################################################################
#create a list of unique identifier for each row of the tables
IDpeptidelistA<-paste(peptidelist_A$PROTEIN,peptidelist_A$PEPTIDE,sep=".")
IDpeptidelistB<-paste(peptidelist_B$PROTEIN,peptidelist_B$PEPTIDE,sep=".")

##########################################################################################
#Find intersections, common unique identifiers
interseca<-subset(peptidelist_A,IDpeptidelistA%in%IDpeptidelistB)

#intersecaINV<-subset(peptidelist_B,IDpeptidelistB%in%IDpeptidelistA)

##########################################################################################
#Find intersections, common proteins, counts the retained fraction of each dataset 
commonPROTEINS<-as.character(unique(interseca$PROTEIN))
retainedA<-length(commonPROTEINS)/length(unique(peptidelist_A$PROTEIN))
retainedB<-length(commonPROTEINS)/length(unique(peptidelist_B$PROTEIN))
paste("% OF RETAINED A:", retainedA, sep = " ")
paste("% OF RETAINED B:", retainedB, sep = " ")

##########################################################################################
#Find intersections, common peptides, counts the retained fraction of each dataset 
retainedA_pep<-nrow(interseca)/nrow(peptidelist_A)
retainedB_pep<-nrow(interseca)/nrow(peptidelist_B)
paste("% OF RETAINED peptides A:", retainedA_pep, sep = " ")
paste("% OF RETAINED peptides B:", retainedB_pep, sep = " ")

##########################################################################################GetOUT
#new table in the same dataset format containing only common protein-peptide pairs
write.csv(interseca, file = 'commonPeptides.csv', row.names=FALSE, quote=FALSE)
#table containing list of common proteins 
write.csv(commonPROTEINS, file = 'commonProteins.csv', row.names=FALSE, quote=FALSE)

####Tschuss