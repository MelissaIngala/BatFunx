####### Random Forests

library("randomForest")
library("plyr") # for the "arrange" function
#install.packages('rfUtilities')
library("rfUtilities") # to test model significance
library("caret") # to get leave-one-out cross-validation accuracies and also contains the nearZeroVar function 

otu_table <- as.data.frame(otu_table(physeq.sub))  
metadata <- physeq.sub@sam_data

dim(otu_table)
#[1] 448 545
dim(metadata)
#[1] 545  10

###### Preprocessing: Remove rare features
otu_nonzero_counts <- apply(otu_table, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")

remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

#Remove Functional OTUs present in less than 10% of samples
otu_table_rare_removed <- remove_rare(table=otu_table, cutoff_pro=0.1)

dim(otu_table_rare_removed)
#[1] 400 545

# Re-normalize so all values are simplex
otu_table_rare_removed_norm <- sweep(otu_table_rare_removed, 2, colSums(otu_table_rare_removed) , '/')*100
colSums(otu_table_rare_removed_norm)

#Transforming the Data to Z score
otu_table_scaled <- scale(otu_table_rare_removed_norm, center = TRUE, scale = TRUE) 

######## Building the RF Model
otu_table_scaled_Binary <- data.frame(t(otu_table_scaled))  
otu_table_scaled_Binary$Binary <- metadata[rownames(otu_table_scaled_Binary), "Binary"]  

set.seed(151)
RF_niche_classify <- randomForest( x=otu_table_scaled_niche[,1:(ncol(otu_table_scaled_niche)-1)] , y=as.factor(otu_table_scaled_niche[ , ncol(otu_table_scaled_niche)]) , ntree=10001, importance=TRUE, proximities=TRUE )
RF_niche_classify
#OOB estimate of  error rate: 13.21%
#Confusion matrix:
#                   Animal_feeder Omnivorous Plant_feeder class.error
#Animal_feeder           386          0           12  0.03015075
#Omnivorous                2          0            4  1.00000000
#Plant_feeder             54          0           87  0.38297872


######## Building the RF Model
otu_table_scaled_niche <- data.frame(t(otu_table_scaled))  
otu_table_scaled_niche$FN <- metadata[rownames(otu_table_scaled_niche), "FeedingNiche"]  

set.seed(152)
RF_FN_classify <- randomForest(x=otu_table_scaled_niche[,1:(ncol(otu_table_scaled_niche)-1)] , y=as.factor(otu_table_scaled_niche[ , ncol(otu_table_scaled_niche)]) , ntree=10001, importance=TRUE, proximities=TRUE )
RF_FN_classify

#OOB estimate of  error rate: 15.6%
#Confusion matrix:
#                Carnivore Frugivore Insectivore Omnivore Piscivore Sanguivore
#Carnivore           0         0           6        0         0          0
#Frugivore           0        92          49        0         0          0
#Insectivore         0        14         346        0         0          0
#Omnivore            0         4           2        0         0          0
#Sanguivore          0         1           8        0         0         22

#class.error
#Carnivore    1.00000000
#Frugivore    0.34751773
#Insectivore  0.03888889
#Omnivore     1.00000000
#Sanguivore   0.29032258

############ Cross-Validate Models by Leave-one-out
require(caret)
fit_control <- trainControl( method = "LOOCV" )    

RF_Binary_classify_loocv <- train( otu_table_scaled_Binary[,1:(ncol(otu_table_scaled_Binary)-1)] , y=otu_table_scaled_Binary[, ncol(otu_table_scaled_Binary)] , method="rf", ntree=501 , tuneGrid=data.frame( mtry=25 ) , trControl=fit_control )
#Resampling: Leave-One-Out Cross-Validation 
#Summary of sample sizes: 544, 544, 544, 544, 544, 544, ... 
#Resampling results:
  
#  Accuracy  Kappa    
#0.866055  0.6258311

#Tuning parameter 'mtry' was held constant at a value of 25

RF_niche_classify_loocv <- train( otu_table_scaled_niche[,1:(ncol(otu_table_scaled_niche)-1)] , y=otu_table_scaled_niche[, ncol(otu_table_scaled_niche)] , method="rf", ntree=501 , tuneGrid=data.frame( mtry=25 ) , trControl=fit_control )

#Resampling results:
  
#  Accuracy   Kappa    
#0.8419118  0.6469638

#########Pull out important variables
#Binary
par(mfrow=c(1,2))
RF_binary_classify_imp <- as.data.frame( RF_niche_classify$importance )
RF_binary_classify_imp$features <- rownames( RF_binary_classify_imp )
RF_binary_classify_imp_sorted <- arrange( RF_binary_classify_imp  , desc(MeanDecreaseAccuracy)  )
barplot(RF_binary_classify_imp_sorted$MeanDecreaseAccuracy, ylab="Mean Decrease in Accuracy (Variable Importance)", main="RF Classification Variable Importance Distribution")

#Discrete
RF_FN_classify_imp <- as.data.frame( RF_FN_classify$importance )
RF_FN_classify_imp$features <- rownames( RF_FN_classify_imp )
RF_FN_classify_imp_sorted <- arrange( RF_FN_classify_imp  , desc(MeanDecreaseAccuracy)  )
barplot(RF_FN_classify_imp_sorted$MeanDecreaseAccuracy, ylab="Mean Decrease in Accuracy (Variable Importance)", main="RF Classification Variable Importance Distribution")


#######Top 10 features
par(mar=c(10,4,4,4))
par(mfrow=c(1,2))
barplot(RF_binary_classify_imp_sorted[1:10,"MeanDecreaseAccuracy"], names.arg=RF_binary_classify_imp_sorted[1:10,"features"] , ylab="Mean Decrease in Accuracy", las=2, ylim=c(0,0.02), cex.names = 0.7, main="Binary Niche RF", col ="#B12A90")  
barplot(RF_FN_classify_imp_sorted[1:10,"MeanDecreaseAccuracy"], names.arg=RF_FN_classify_imp_sorted[1:10,"features"] , ylab="Mean Decrease in Accuracy", las=2, ylim=c(0,0.025), cex.names = 0.7, main="Discrete Niche RF", col ="#FCA636")  

saveRDS( file = "RF_binary_model.rda" , RF_niche_classify )
saveRDS( file = "RF_discrete_model.rda" , RF_FN_classify )

rmarkdown::render("randomForest_niches.R")
rmarkdown::render("randomForest_niches.R", "pdf_document")

######## Building the RF Model-- Family
otu_table_scaled_family <- data.frame(t(otu_table_scaled)) 
metadata<-as.data.frame(metadata)
otu_table_scaled_family$FAM <- metadata$HostFamily  

set.seed(155)
RF_FAM_classify <- randomForest(x=otu_table_scaled_family[,1:(ncol(otu_table_scaled_family)-1)] , y=as.factor(otu_table_scaled_family[ , ncol(otu_table_scaled_family)]) , ntree=10001, importance=TRUE, proximities=TRUE )
RF_FAM_classify

############## Genus
otu_table_scaled_genus <- data.frame(t(otu_table_scaled)) 
otu_table_scaled_genus$GEN <- metadata$HostGenus  

set.seed(156)
RF_GEN_classify <- randomForest(x=otu_table_scaled_genus[,1:(ncol(otu_table_scaled_genus)-1)] , y=as.factor(otu_table_scaled_genus[ , ncol(otu_table_scaled_genus)]) , ntree=10001, importance=TRUE, proximities=TRUE )
RF_GEN_classify
