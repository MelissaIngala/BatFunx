########### Ingala et al. 2021
###### You are more than what you eat: differential enrichment of 
# microbiome functions across bat dietary guilds.
# R version: 4.0.0

#Load packages
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
library(phyloseq)
#packageVersion("phyloseq")
library("ggplot2")
#packageVersion("ggplot2")
library(vegan)
library(picante)
library(ape)
#install.packages("devtools") # if not already installed
library("devtools")
#devtools::install_github("biomformat", "joey711")
library(biomformat)
#biocLite("DESeq2")
library(viridis)

#Import the BIOM file, metadata and taxonomy included
dat <- read_biom(biom_file = "~/Downloads/merged-table.w_smd.biom")
otu_table <- as.data.frame(as.matrix(biom_data(dat)))
taxonomy <- observation_metadata(dat)
metadata <- sample_metadata(dat)

#Merge into phyloseq obj
SAM<-sample_data(metadata)
TAX<-tax_table(as.matrix(taxonomy))
OTU<- otu_table(otu_table, taxa_are_rows=TRUE)
physeq<-merge_phyloseq(OTU, TAX, SAM)
colnames(tax_table(physeq))<-"MetaCyc_pathway"

#Relative abundance
#Calculate Relative Abundance
relative  = transform_sample_counts(physeq, function(OTU) OTU / sum(OTU))

#Visually inspect
#Barplot of top 25 OTUs
Top25OTUs = names(sort(taxa_sums(relative), TRUE)[1:25])
comparetop25 = prune_taxa(Top25OTUs, relative)
plot_bar(comparetop25, fill = "MetaCyc_pathway", title = "Functional Profile by Sample")

###################### ORDINATION #########################
#Ordination plots of microbiomes 
#Scale
set.seed(42)

#Hellinger transform the data
require(microbiome)
physeq.trans<-transform(x = physeq, transform = "hellinger", target = "OTU", shift = 0, scale = 1)
physeq.trans<-subset_samples(physeq = physeq.trans, Binary!="Control")
physeq.sub<-subset_samples(physeq = physeq, Binary!="Control")

funx.ord <- ordinate(
  physeq = physeq.trans, 
  method = "PCoA", 
  distance = "bray"
)

# Calculate distance matrix
study.bray <- phyloseq::distance(physeq.trans, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq.trans))

########## Plot ordination
#pdf(file = "Hellinger-funx-PCOA.pdf", width = 6, height = 5)
plot_ordination(
  physeq = physeq.trans,
  ordination = funx.ord,
  axes = c(1,2), 
  color = "FeedingNiche", 
  title = "Ordination of Metagenome Functions") +
  scale_color_viridis_d(option="plasma", direction = -1)+
  geom_point(aes(color = FeedingNiche, shape = Binary), size = 2.5) +
  theme_bw() +
  stat_ellipse() + geom_jitter()
#dev.off()

# Adonis test
#install.packages("remotes")
#remotes::install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)
sampledf$FeedingNiche<-as.factor(sampledf$FeedingNiche)
sampledf$Binary<-as.factor(sampledf$Binary)

adonis(study.bray ~ FeedingNiche * HostSpecies + HostGenus + HostFamily, data = sampledf) #F = 8.6712, r2 = 0.06251, p = 0.001 *** 
adonis(study.bray ~ Binary * HostSpecies + HostGenus + HostFamily, data = sampledf) #F = 9.2791, r2 = 0.02676, P = 0.001 ***
adonis_res<-adonis.pair(study.bray, sampledf$FeedingNiche, nper = 1000, corr.method = "fdr")
write.table(x = adonis_res, file = "permanova_dietniche.csv", sep = ",")
beta <- betadisper(study.bray, sampledf$FeedingNiche)
permutest(beta) #p = 0.061

#microbiome functions are significantly but weakly different.
#this is probably due to large proportion of shared "housekeeping" functions
#See LEfSe analyses for discriminating functions analysis.

################### HOST PHYLOGENY #####################
#This is a trimmed host phylogeny from VertLife (Upham et al. 2019)
library(phytools)
############## MRMs on host phylo, microbiomes, diet matrices
install.packages("ecodist")
library(ecodist)

########## BRAY-CURTIS DISSIMILIARITY MATRIX OF METAGENOME FUNCTIONS
#compute Bray Curtis distance matrix on merged dataset (drop control first)
subs.physeq<-subset_samples(physeq = physeq, Binary!="Control")

merged.bray.fun<-phyloseq::distance(mergedspp, method = "bray")
#writemydf<-as.data.frame(as.matrix(merged.bray.fun))
#write.csv(writemydf,file = "batsppdists_corrected.csv") #correct taxon names manually

A <- read.table(file="batsppdists_corrected.csv", sep = ",", row.names = 1, header=T)
merged.fun.dist <- data.matrix(A)

#########REORDER
merged.fun.dist2<-merged.fun.dist[order(rownames(merged.fun.dist)), order(colnames(merged.fun.dist))]

########### PATRISTIC DISTANCES - PHYLOGENY
#compute patristic distance matrix for MRMs
library(ape)
PatristicDistMatrix<-cophenetic(tree_6237)
phylo.dist<-as.data.frame(as.matrix(PatristicDistMatrix))
#write.csv(phylo.dist,file = "batphylodists.csv")

#######REORDER: ESSENTIAL FOR MRM
PatristicDistMatrix<-PatristicDistMatrix[order(rownames(PatristicDistMatrix)), order(colnames(PatristicDistMatrix))]

####### ECOLOGICAL SPP LEVEL DIETS- ELTONTRAITS
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
#write.csv(signifdf,file = "subset-signif-funx-braydists.csv") #correct taxon names manually

B <- read.table(file="subset-signif-funx-braydists.csv", sep = ",", row.names = 1, header=TRUE)
bray.dist <- data.matrix(B)

subshostecoevomrm<-MRM(formula = as.dist(bray.dist) ~ as.dist(bat.diets.matrix) + as.dist(PatristicDistMatrix),
                       nperm = 15000, method = 'linear')

################# Perform Phylogenetic Comparative Anlayses
#### modeling pathways as traits of the host

install.packages("geiger")
library(geiger)
library(phytools)
library(phyloseq)
library(viridis)
library(microbiome)

#Import phylogenetic tree (Upham et al. 2019)
trees<-read.nexus("~/Documents/BigBatProj/tree-pruner-7beec7de-2578-4bce-b87b-1e80d044f018/output.nex")
#rf.tree<-averageTree(trees,method="branch.score.difference")
tree_6237<-trees$tree_9764
tree_6237<-force.ultrametric(tree_6237)
plotTree(tree_6237, node.numbers = F)

#Create vector of pathways most influential in RF decision tree model gain
top.pwys<-c("P125-PWY", "PWY-6612", "LACTOSECAT-PWY", "OANTIGEN-PWY", "P164-PWY",
            "DTDPRHAMSYN-PWY", "FASYN-ELONG-PWY", "FOLSYN-PWY", "PWY-1269",
            "BRANCHED-CHAIN-AA-SYN-PWY")

#Prune merged functional feature table to just these pathways for the analysis
pwy_of_int = subset_taxa(mergedspp, MetaCyc_pathway %in% top.pwys)
plot_bar(physeq = pwy_of_int, fill = "MetaCyc_pathway")
#Transform data to normalize, center log ratio
pwys.clr <- microbiome::transform(pwy_of_int, "clr")

#manually correct species names to match phylogeny (for now, must update txonomy later)
#write.csv(x = as.data.frame(pwys.clr@otu_table), file = "Bat_mb_traits.csv")
Bat_traits<-read.csv(file = "Bat_mb_traits.csv", header = T, sep = ",", row.names = 1)

require(geiger)
tree<-treedata(tree_6237,Bat_traits,sort=T,warnings=T)$phy
pwys<-treedata(tree_6237,Bat_traits,sort=T,warnings=T)$dat
bat_td<-treedata(tree_6237,Bat_traits,sort=T,warnings=T)

######### From Capunitan et al 2020
setwd(dir = "Geiger_outputs/AnalysesOct1920/")
geiger_results <- list()
for (variable in c("P125.PWY", "PWY.6612", "LACTOSECAT.PWY", "OANTIGEN.PWY", "P164.PWY",
                   "DTDPRHAMSYN.PWY", "FASYN.ELONG.PWY", "FOLSYN.PWY", "PWY.1269",
                   "BRANCHED.CHAIN.AA.SYN.PWY")){
  
  states<-bat_td$data[,variable] 
  tree<-bat_td$phy
  states<-states[tree$tip.label]
  
  fit_bm<-fitContinuous(tree,states,model="BM") 
  fit_ou<-fitContinuous(tree,states,model="OU") 
  fit_eb<-fitContinuous(tree,states,model="EB") 
  fit_wn<-fitContinuous(tree,states,model="white")
  
  aic_bm<-fit_bm$opt$aicc 
  aic_ou<-fit_ou$opt$aicc
  aic_eb<-fit_eb$opt$aicc 
  aic_wn<-fit_wn$opt$aicc
  
  aic_all<-c(aic_bm,aic_ou,aic_eb, aic_wn) 
  names(aic_all)<-c("BM","OU","EB","white")
  
  aicw2<-aicw(aic_all) 
  print(variable) 
  print(aicw2)
  
  fit_lambda<-fitContinuous(tree,states,model="lambda")
  
  geiger_results[[variable]]<-cbind(aicw2, fit_lambda$opt$lambda, fit_ou$opt$alpha)
  
  print(paste("Full parameter estimates for ", variable,':'))
  print(fit_bm) 
  print(fit_ou) 
  print(fit_eb) 
  print(fit_wn) 
  print('lambda')
  
  print(fit_lambda)
  
  #to view variable on phylogeny: 
  obj<-contMap(tree,states,fsize=0.4)
  pdf(paste(variable,".defaultcols.pdf")) 
  #obj<-setMap(obj,c("#0D0887FF","#CC4678FF","#F89441FF","#F0F921FF"))
  plot(obj, fsize= 0.4)
  dev.off()
}

geiger_results

ggplot(geiger_results[[1]], aes(fill=as.factor(rownames(geiger_results[1]$P125.PWY)), y=w, x=names(geiger_results[1]))) + 
  geom_bar(position="fill", stat="identity") + theme_bw() + scale_fill_viridis_d(option = "C") +
  geom_point(aes(y=`fit_lambda$opt$lambda`), stat="identity",
             position="dodge",
             alpha=.8,
             size=5, color = "black") +
  scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Pagel's Lambda [λ]"))


for (i in seq_along(geiger_results)) {
  p1<-ggplot(geiger_results[[i]], aes_string(x="names(geiger_results[i])",y="w", fill=as.factor(rownames(geiger_results[1]$P125.PWY)))) +
    geom_bar(position="fill", stat="identity") + theme_bw() + scale_fill_viridis_c() +
    geom_point(aes(y=`fit_lambda$opt$lambda`), stat="identity",
               position="dodge",
               alpha=.8,
               size=5, color = "black") +
    scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Pagel's Lambda [λ]"))
  print(p1)
}

pdf(file = "comparative_barplots_1020.pdf")
lapply(1:length(geiger_results), function(i){
  ggplot(geiger_results[[i]], aes_string(x="names(geiger_results[i])",y="w", fill=as.factor(rownames(geiger_results[1]$P125.PWY)))) +
    geom_bar(position="fill", stat="identity") + theme_bw() + scale_fill_viridis_d() +
    geom_point(aes(y=`fit_lambda$opt$lambda`), stat="identity",
               position="dodge",
               alpha=.8,
               size=5, color = "black") +
    scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Pagel's Lambda [λ]"))})
dev.off()

library(reshape2)
df.long<-melt(geiger_results,id.vars=c("fit","delta", "w", "fit_lambda$opt$lambda", "fit_ou$opt$alpha"))
df.long
write.csv(x = df.long, file = "Geiger_results_LONG.csv")
df.long<-read.csv(file = "Geiger_results_LONG.csv")

pdf(file = "Barplots_allmodels.pdf", width = 12, height = 4)
ggplot(df.long, aes(x=L1,y=w,fill=as.factor(Model))) +
  geom_bar(position="fill", stat="identity") + theme_bw() + scale_fill_viridis_d() +
  geom_point(aes(y=df.long$fit_lambda.opt.lambda), stat="identity",
             position="dodge",
             alpha=.8,
             size=5, color = "black") +
  scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Pagel's Lambda [λ]")) + theme(axis.text.x = element_text(angle = 45, hjust = 0.95))
dev.off()
