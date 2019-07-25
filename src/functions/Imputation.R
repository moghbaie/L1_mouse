## Mehrnoosh Oghbaie
## 07/24/2019
## Imputate missing values

## Imputation consists of three stages:
##  1. Remove proteins not identified in both cases and controls
##  2. Impute values for proteins with zero intensity (cases and controls)
##  3. Impute values for proteins that have at least one non-zero intensity


set.seed(123)
#####################################################################################################
##  1. Remove proteins not identified in both cases and controls

Template$set("public","PreImputation", list())

Template$set("public","removeAllZeros", function(){
  self$PreImputation <- self$experiment
  for ( i in names(self$PreImputation)){
    self$PreImputation[[i]] <- self$PreImputation[[i]][rowSums(self$PreImputation[[i]][,-(1:3)], na.rm =T)!=0,]
  }
}
)


#####################################################################################################
### Imputation method 7_1 (* We implement this method)

#Method_1:
# For each replica mean and standard deviation of non-zero proteins intensities are counted. 
# New intensity for each missing value in each replica is sampled from uniform distribution with parameters:
#  start = mu - 3*sd, end = mu - 2*sd
# Intnew=Unif(mean(Intreplica)-3*sd(Intreplica),mean(Intreplica)-2*sd(Intreplica))

#Method_7:
# For outliers (proteins, which have zero values in all but one replica) this method implies one of the methods for imputation of all zero replicas.
# For other proteins uses method6

Template$set("public","experimentImputed", list())

Template$set("public","imputeAll", function(){
  self$experimentImputed <- self$PreImputation
  for(name in names(self$experiment)){
    coln <- colnames(self$experimentImputed[[name]])
    self$experimentImputed[[name]] <- impute(self$experimentImputed[[name]], strsplit(name,".vs.")[[1]], amm = "2", pmm = "6")
  }
}
)




