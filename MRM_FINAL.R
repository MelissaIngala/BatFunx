############## MRMs on host phylo, microbiomes, diet matrices

install.packages("ecodist")
library(ecodist)

########## BRAY-CURTIS DISSIMILIARITY MATRIX OF METAGENOME FUNCTIONS
#compute Bray Curtis distance matrix on merged dataset (drop control first)
subs.physeq<-subset_samples(physeq = physeq, Binary!="Control")

merged.bray.fun<-phyloseq::distance(mergedspp, method = "bray")
writemydf<-as.data.frame(as.matrix(merged.bray.fun))
write.csv(writemydf,file = "batsppdists_corrected.csv") #correct taxon names manually

A <- read.table(file="batsppdists_corrected.csv", sep = ",", row.names = 1, header=T)
merged.fun.dist <- data.matrix(A)

#########REORDER
merged.fun.dist2<-merged.fun.dist[order(rownames(merged.fun.dist)), order(colnames(merged.fun.dist))]


########### PATRISTIC DISTANCES - PHYLOGENY
#compute patristic distance matrix for MRMs
library(ape)
PatristicDistMatrix<-cophenetic(tree_6237)
phylo.dist<-as.data.frame(as.matrix(PatristicDistMatrix))
write.csv(phylo.dist,file = "batphylodists.csv")
#######REORDER: ESSENTIAL FOR MRM
PatristicDistMatrix<-PatristicDistMatrix[order(rownames(PatristicDistMatrix)), order(colnames(PatristicDistMatrix))]
plot(hclust(PatristicDistMatrix), cex=0.4)

####### ECOLOGICAL AND BODY MASS- ELTONTRAITS
#convert diet percentages into distance matrix
diets<-read.delim(file = "BatDiet_EltonTraits.txt", sep = "\t", header = T, row.names = 2)
bat.diets.percents<-diets[,3:10]
require(scales)
bat.diets.dist<-vegan::vegdist(bat.diets.percents, metric='jaccard')
options(repr.plot.width=9, repr.plot.height=4)
plot(hclust(bat.diets.dist), cex=0.4)
bat.diets.matrix<-as.matrix(bat.diets.dist)
##########REORDER
bat.diets.matrix2<-bat.diets.matrix[order(rownames(bat.diets.matrix)), order(colnames(bat.diets.matrix))]

######################### CONSTRUCT AND RUN MRMS #######################

hostecoevomrm<-MRM(formula = as.dist(merged.fun.dist2) ~ as.dist(bat.diets.matrix2) + as.dist(PatristicDistMatrix),
    nperm = 100000, method = 'linear')

#not significant for host phylogeny or diet
#possible that many functions are simply housekeeping, not driving diffs

#Pull out functions significant according to Binary LEfSE
LDA.signifpwy<-read.csv(file = 'Binary_Pathways_annotated.csv', header = T)
signifpwys<- LDA.signifpwy$MetaCyc_Pathway

physeq.subset <- subset_taxa(mergedspp, rownames(tax_table(mergedspp)) %in% signifpwys)
bray.curtis.signif<-phyloseq::distance(physeq.subset, method = "bray")
signifdf<-as.data.frame(as.matrix(bray.curtis.signif))
write.csv(signifdf,file = "subset-signif-funx-braydists.csv") #correct taxon names manually

B <- read.table(file="subset-signif-funx-braydists.csv", sep = ",", row.names = 1, header=TRUE)
bray.dist <- data.matrix(B)
  
subshostecoevomrm<-MRM(formula = as.dist(bray.dist) ~ as.dist(bat.diets.matrix) + as.dist(PatristicDistMatrix) + as.dist(mass),
                   nperm = 15000, method = 'linear')