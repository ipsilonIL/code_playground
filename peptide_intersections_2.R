##### rearrange column positions in a csv file
rm(list=ls())
setwd("/Users/piazza/crap for R/")

library("stringr", lib.loc="\\\\nas-green-1/share-green-appv-1-$/r_library")
library("stringr")

##########################################################################################GetIN
peptidelist_A <- read.csv('L-Phe_2mM.csv') 
peptidelist_B <- read.csv('L-Phe_20mM.csv') 
head(peptidelist_A)
head(peptidelist_B)

##########################################################################################
#create a list of unique identifier for each row of the tables. N.B: STRIPPED peptide sequences. 
IDpeptidelistA<-paste(peptidelist_A$proteinName,peptidelist_A$peptide,sep=".")
IDpeptidelistB<-paste(peptidelist_B$proteinName,peptidelist_B$peptide,sep=".")
peptidelist_A<-cbind(peptidelist_A, IDpeptidelistA)
peptidelist_B<-cbind(peptidelist_B, IDpeptidelistB)
  
    #Stripped peptide list returns a table with unique peptides - stripped-
StrippeptidelistA_ix<-!duplicated(peptidelist_A$IDpeptidelistA)
StrippeptidelistA<-peptidelist_A[StrippeptidelistA_ix , ]

StrippeptidelistB_ix<-!duplicated(peptidelist_B$IDpeptidelistB)
StrippeptidelistB<-peptidelist_B[StrippeptidelistB_ix , ]

##########################################################################################
##Subset of stripped peptides at C1 validated by DIA
DIAvalidStrippeptidelistA_ix<-StrippeptidelistA$DIA.validation==1
DIAvalidStrippeptidelistA<-StrippeptidelistA[DIAvalidStrippeptidelistA_ix, ]
##Subset of stripped peptides at C2 validated by DIA
DIAvalidStrippeptidelistB_ix<-StrippeptidelistB$DIA.validation==1
DIAvalidStrippeptidelistB<-StrippeptidelistB[DIAvalidStrippeptidelistB_ix, ]

##########################################################################################
#Find intersections, common unique stripped peptides bw C1 and C2
intersecaSTRI<-subset(StrippeptidelistA,StrippeptidelistA$IDpeptidelistA%in%StrippeptidelistB$IDpeptidelistB)
#intersecaINV<-subset(peptidelist_B,IDpeptidelistB%in%IDpeptidelistA)

#Subsets of common unique stripped peptides bw C1 and C2 AND validated by DIA
DIAvalid<-intersecaSTRI$DIA.validation==1
commonValPep<-intersecaSTRI[DIAvalid, ]

##########################################################################################
#Counts the % of DIA validated PROTEINS over shotgun hits found at C1 and C2
uniquesvA<-DIAvalidStrippeptidelistA[!duplicated(DIAvalidStrippeptidelistA$proteinName), ]
uniqueProteinsvA<-uniquesvA[ ,3:4]
commonPROTEINS_vA<-as.character(uniqueProteinsvA[,1])

uniquesvB<-DIAvalidStrippeptidelistB[!duplicated(DIAvalidStrippeptidelistB$proteinName), ]
uniqueProteinsvB<-uniquesvB[ ,3:4]
commonPROTEINS_vB<-as.character(uniqueProteinsvB[,1])

uniquesA<-StrippeptidelistA[!duplicated(StrippeptidelistA$proteinName), ]
uniqueProteinsA<-uniquesA[ ,3:4]
commonPROTEINS_A<-as.character(uniqueProteinsA[,1])

uniquesB<-StrippeptidelistB[!duplicated(StrippeptidelistB$proteinName), ]
uniqueProteinsB<-uniquesB[ ,3:4]
commonPROTEINS_B<-as.character(uniqueProteinsB[,1])

#commonPROTEINS<-as.character(unique(commonValPep$proteinName))
#commonDescriptors<-as.character(unique(commonValPep$proteinDescription))

retainedA<-length(commonPROTEINS_vA)/length(commonPROTEINS_A)
retainedB<-length(commonPROTEINS_vB)/length(commonPROTEINS_B)
testo<-paste("% of DIA validated shotgun protein hits at C1:", retainedA, sep = " ")
write(testo, file="PeptideNumbers", append=TRUE, sep="\t")
testo<-paste("% of DIA validated shotgun protein hits at C2:", retainedB, sep = " ")
write(testo, file="PeptideNumbers", append=TRUE, sep="\t")
#########################################################################################
#Find intersections, counts the % of common DIA validated PROTEINS bw C1 and C2
uniques<-commonValPep[!duplicated(commonValPep$proteinName), ]
uniqueProteins<-uniques[ ,3:4]
commonPROTEINS<-as.character(uniqueProteins[,1])

IntretainedA<-length(commonPROTEINS)/length(commonPROTEINS_vA)
IntretainedB<-length(commonPROTEINS)/length(commonPROTEINS_vB)
testo<-paste("% of common protein over total protein hits at C1 (DIA validated):", IntretainedA, sep = " ")
write(testo, file="PeptideNumbers", append=TRUE, sep="\t")
testo<-paste("% of common DIA validated protein over total protein hits at C2 (DIA validated):", IntretainedB, sep = " ")
write(testo, file="PeptideNumbers", append=TRUE, sep="\t")
##########################################################################################
#Find intersections, counts the % of DIA validated PEPTIDES at C1 and C2
retainedA_pep<-nrow(commonValPep)/nrow(DIAvalidStrippeptidelistA)
retainedB_pep<-nrow(commonValPep)/nrow(DIAvalidStrippeptidelistB)
testo<-paste("% OF DIA validated shotgun peptide hits at C1:", retainedA_pep, sep = " ")
write(testo, file="PeptideNumbers", append=TRUE, sep="\t")
testo<-paste("% OF DIA validated shotgun peptide hits at C2:", retainedB_pep, sep = " ")
write(testo, file="PeptideNumbers", append=TRUE, sep="\t")
##########################################################################################
#Find intersections, common peptides from the DDA list, counts the retained fraction of each dataset 
retainedDDA<-nrow(intersecaSTRI)/nrow(StrippeptidelistA)
ExtraDDA_B<-1-(nrow(intersecaSTRI)/nrow(StrippeptidelistB))
testo<-paste("% peptides retained at higher concentration (shotgun data only):", retainedDDA, sep = " ")
write(testo, file="PeptideNumbers", append=TRUE, sep="\t")
testo<-paste("% peptides found at C2 but not at C1 (shotgun data only):", ExtraDDA_B, sep = " ")
write(testo, file="PeptideNumbers", append=TRUE, sep="\t")
##########################################################################################
#Find intersections, common DIA validated peptides, counts the retained fraction of each dataset 
retainedDIA<-nrow(commonValPep)/nrow(DIAvalidStrippeptidelistA)
ExtraDIA_B<-1-(nrow(commonValPep)/nrow(DIAvalidStrippeptidelistB))
testo<-paste("% DIA validated peptides retained at higher concentration:", retainedDIA, sep = " ")
write(testo, file="PeptideNumbers", append=TRUE, sep="\t")
testo<-paste("% DIA validated peptides found at C2 but not at C1:", ExtraDIA_B, sep = " ")
write(testo, file="PeptideNumbers", append=TRUE, sep="\t")
##########################################################################################GetOUT
#Returns a table of stripped peptides from C1:
write.csv(StrippeptidelistA, file = 'StrippedC1.csv', row.names=FALSE, quote=FALSE)
#Returns a table of stripped peptides from C2:
write.csv(StrippeptidelistB, file = 'StrippedC2.csv', row.names=FALSE, quote=FALSE)
#Returns a table of -DIA validated - stripped peptides from C1:
write.csv(DIAvalidStrippeptidelistA, file = 'StrippedDIAValC1.csv', row.names=FALSE, quote=FALSE)
#Returns a table of -DIA validated -stripped peptides from C2:
write.csv(DIAvalidStrippeptidelistB, file = 'StrippedDIAValC2.csv', row.names=FALSE, quote=FALSE)

#tableS containing list of common proteins 
#new table in the same dataset format containing only common protein-peptide pairs
write.csv(commonValPep, file = 'commonStrValPeptides.csv', row.names=FALSE, quote=FALSE)
#table containing list of common proteins 
write.csv(uniqueProteins, file = 'commonValProteins.csv', row.names=FALSE, quote=FALSE)

####Tschuss