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
