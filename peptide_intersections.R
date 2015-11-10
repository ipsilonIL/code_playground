##### rearrange column positions in a csv file
rm(list=ls())
setwd("/Users/piazza/Desktop/")

library("stringr", lib.loc="\\\\nas-green-1/share-green-appv-1-$/r_library")
library("stringr")

##########################################################################################GetIN
peptidelist_A <- read.csv('FBP_1mM.csv') 
peptidelist_B <- read.csv('FBP_25mM.csv') 
head(peptidelist_A)
head(peptidelist_B)

##########################################################################################
#create a list of unique identifier for each row of the tables. N.B: STRIPPED peptide sequences. 
IDpeptidelistA<-paste(peptidelist_A$proteinName,peptidelist_A$peptide,sep=".")
IDpeptidelistB<-paste(peptidelist_B$proteinName,peptidelist_B$peptide,sep=".")
StripIDpeptidelistA<-unique(IDpeptidelistA)
StripIDpeptidelistB<-unique(IDpeptidelistB)
##########################################################################################
#Find intersections, common unique identifiers
#interseca<-subset(peptidelist_A,IDpeptidelistA%in%IDpeptidelistB)
intersecaSTRI<-subset(peptidelist_A,StripIDpeptidelistA%in%StripIDpeptidelistB)

#intersecaINV<-subset(peptidelist_B,IDpeptidelistB%in%IDpeptidelistA)

##########################################################################################
#Subsets peptides validated by DIA
DIAvalid<-intersecaSTRI$DIA.validation==1
commonValPep<-intersecaSTRI[DIAvalid, ]

##########################################################################################
#Find intersections, counts the % of DIA validated proteins at C1 and C2

#commonPROTEINS<-as.character(unique(commonValPep$proteinName))
#commonDescriptors<-as.character(unique(commonValPep$proteinDescription))
uniques<-commonValPep[!duplicated(commonValPep$proteinName), ]
uniqueProteins<-uniques[ ,3:4]
commonPROTEINS<-as.character(uniqueProteins[,1])

retainedA<-length(commonPROTEINS)/length(unique(peptidelist_A$proteinName))
retainedB<-length(commonPROTEINS)/length(unique(peptidelist_B$proteinName))
paste("% of DIA validated proteins at C1:", retainedA, sep = " ")
paste("% of DIA validated proteins at C2:", retainedB, sep = " ")

##########################################################################################
#Find intersections, counts the % of DIA validated peptides at C1 and C2
retainedA_pep<-nrow(commonValPep)/length(StripIDpeptidelistA)
retainedB_pep<-nrow(commonValPep)/length(StripIDpeptidelistB)
paste("% of DIA validated peptides at C1:", retainedA_pep, sep = " ")
paste("% OF DIA validated peptides at C2:", retainedB_pep, sep = " ")

##########################################################################################
#Find intersections, common peptides from the DDA list, counts the retained fraction of each dataset 
retainedDDA<-nrow(intersecaSTRI)/length(StripIDpeptidelistA)
ExtraDDA_B<-1-(nrow(intersecaSTRI)/length(StripIDpeptidelistB))
paste("% peptides retained at higher concentration:", retainedDDA, sep = " ")
paste("% peptides found at C2 but not at C1:", ExtraDDA_B, sep = " ")

##########################################################################################
#Find intersections, common DIA validated peptides, counts the retained fraction of each dataset 
DIAvalidpeptidelist_A_ix<-peptidelist_A$DIA.validation==1
DIAvalidpeptidelist_B_ix<-peptidelist_B$DIA.validation==1
DIAvalidpeptidelist_A<-peptidelist_A[DIAvalidpeptidelist_A_ix, ]
DIAvalidpeptidelist_B<-peptidelist_B[DIAvalidpeptidelist_B_ix, ]

retainedDIA<-nrow(commonValPep)/length(DIAvalidpeptidelist_A)
ExtraDIA_B<-1-(nrow(commonValPep)/length(DIAvalidpeptidelist_B))
paste("% DIA validated peptides retained at higher concentration:", retainedDIA, sep = " ")
paste("% DIA validated peptides found at C2 but not at C1:", ExtraDIA_B, sep = " ")

##########################################################################################GetOUT
#new table in the same dataset format containing only common protein-peptide pairs
write.csv(commonValPep, file = 'commonPeptides.csv', row.names=FALSE, quote=FALSE)
#table containing list of common proteins 
write.csv(commonPROTEINS, file = 'commonProteins.csv', row.names=FALSE, quote=FALSE)


newtable <- data.frame(commonPROTEINS, commonDescriptors)

####Tschuss