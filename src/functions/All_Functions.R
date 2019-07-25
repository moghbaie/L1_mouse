# Mehrnoosh Oghbaie
# 07/24/2019
# Repository for all the functions
# This file includes most of the general functions that are going to be used in this project

######################################################################################
# Either download or install the required library from CRAN or bioconductor
#######################################################################################

install.packages.if.necessary <- function(CRAN.packages=c(), bioconductor.packages=c()) {
  #if (length(bioconductor.packages) > 0) {
  #  source("http://bioconductor.org/biocLite.R")
  #}
  
  for (p in bioconductor.packages) {
    if (!require(p, character.only=T)) {
      BiocManager::install(p,version = "3.8") 
    }
    library(p,lib.loc="~/R/win-library/3.5",character.only=T)
  }
  
  for (p in CRAN.packages) {	
    if (!require(p, character.only=T)) { 	
      install.packages(p) 	
    }	
    #library(p, lib.loc="~/R/win-library/3.5",character.only=T) 
    library(p, lib.loc="~/R/win-library/3.5",character.only=T)
  }
}


##################################################################################
### Quality Control functions
##################################################################################
runQC <- function(input){
  txt_folder <-  file.path(dirname(as.character(input[1,"V2"])))
  r = createReport(txt_folder)
  cat(paste0("\nReport generated as '", r$report_file, "'\n\n"))
}

###################################################################################
### Imputation functions
##################################################################################

### function : finding column with minimum zeros
count_zeros <- function(x){
  return(sum(is.na(x)))
}

min_col <- function(x){
  return(names(x[x==min(x)]))
}
### function: generate uniform distribution for given number, mean and std
na_zeros_impute <- function(na_zeros , mu,sd){
  return(matrix(runif(sum(na_zeros), mu-3*sd, mu-2*sd),ncol=1))
}

perseus_zeros_impute <- function(na_zeros , mu,sd){
  return(matrix(rnorm(sum(na_zeros), mean = mu-1.8*sd, sd = 0.3*sd),ncol=1))
}

#####################################################################################################
### Imputation - all zeros

impute_all_zeros <- function(data,con,amm = "2"){
  ## Predicting all missing replicates values
  colnames <- colnames(data)[grepl(con,colnames(data))]
  na_zeros <- rowSums(data[,colnames],na.rm=T)==0
  min_zero_col <- min_col(apply(data[,colnames],2,count_zeros))
  if(amm == "2"){
    sample_min_zeros <- unname(unlist(data[,min_zero_col]))
    mu <- mean(sample_min_zeros, na.rm=T)
    sd <- sd(sample_min_zeros, na.rm=T)
    data[na_zeros,colnames] <- matrix(na_zeros_impute(na_zeros* length(colnames) , mu,sd),ncol= length(colnames) , nrow = sum(na_zeros))
  } 
  if(amm == "1"){
    sample <- data[,colnames]
    stats <- calculate_stats_nonzeros(data,con)
    data[na_zeros,colnames] <- matrix(cbind(na_zeros_impute(na_zeros , stats[1,1],stats[1,2]),
                                            na_zeros_impute(na_zeros , stats[2,1],stats[2,2]),
                                            na_zeros_impute(na_zeros , stats[3,1],stats[3,2])),
                                      ncol= length(colnames),
                                      nrow = sum(na_zeros))
  }
  return(data)
}

##################################################################################################
### Calculate statistics of all non zeros

calculate_stats_nonzeros <- function(data,con){
  ## Predicting all missing replicates values
  colnames <- colnames(data)[grepl(con,colnames(data))]
  sample <- data[,colnames]
  mu1 <- mean(sample[,1], na.rm=T)
  sd1 <- sd(sample[,1], na.rm=T)
  mu2 <- mean(sample[,2], na.rm=T)
  sd2 <- sd(sample[,2], na.rm=T)
  mu3 <- mean(sample[,3], na.rm=T)
  sd3 <- sd(sample[,3], na.rm=T)
  stats <- rbind(cbind(mu1,sd1),
                 cbind(mu2,sd2),
                 cbind(mu3,sd3))
  return(stats)
}
##################################################################################################
### Impute partial zeros

impute_partial_zeros <- function(data, count_na, colnames,condition, nb=0){
  y <- data[count_na>nb,colnames]
  print(paste("There are ",dim(y)[1]," records missing. \n And ", (dim(data)[1]-dim(y)[1]), " records with full values in ",condition))
  if(length(colnames)>2){
    for(i in 1:dim(y)[1]){
      col_NAs <- colnames(y)[is.na(unlist(y[i,]))]
      for(j in col_NAs){
        col_NA <- j
        col_select <- names(y[i,colnames(y)!=col_NA])
        sample <- data[complete.cases(data[,col_select]),col_select]
        delta <- (sample[,1]-sample[,2])/mean(unlist(unname(sample)))
        mu <- mean(unlist(unname(delta)), na.rm=T)
        std <- sd(unlist(unname(delta)),na.rm=T)
        cor <- cor(data[count_na==0,colnames])
        mean_cor <- mean(cor[col_NA,col_select])
        deltanew <- rnorm(1,mu, std*sqrt(2)/mean_cor)
        y[i,col_NA] <- mean(unlist(unname(y[i,col_select])),na.rm=T)*abs(1+deltanew)
      }
    }
  }else{
    na_zeros <- length(y[is.na(y)])
    mu <- mean(y[!is.na(y)])
    sd <- sd(y[!is.na(y)])
    y[is.na(y)] <- na_zeros_impute(na_zeros , mu,sd)
  }
  data[count_na>0,colnames] <- y
  return(data)
}


impute <- function(data,condition, amm = "2", pmm = "6"){
  for(con in condition){
    stats <- calculate_stats_nonzeros(data,con)
    data <- impute_all_zeros(data,con,amm = "2")
    colnames <- colnames(data)[grepl(con,colnames(data))]
    count_na <- apply(data[,colnames],1, count_zeros)
    
    if(pmm == "6"){
      data <- impute_partial_zeros(data, count_na, colnames, condition, nb = 0)
    }
    
    if(pmm == "7"){
      data <- impute_partial_zeros(data, count_na, colnames, condition, nb = 1)
      data_missing <- apply(data, 2, function(x) is.na(x))
      for( i in 1:length(colnames)){
        data[data_missing[,i],i] <- na_zeros_impute(apply(data_missing,2,sum)[i] , stats[i,1], stats[i,2])
      }
    }
  }
  return(data)
}




### Perseus imputation


impute_Perseus <- function(data,condition){
  if(length(condition)>1){
    condition <- paste0(condition, collapse="|")
  }
  colnames <- colnames(data)[grepl(condition,colnames(data))]
  for(i in 1:length(colnames)){
    zero_col <- is.na(data[,colnames[i]])
    sample <- data[,colnames[i]]
    sample <- sample[complete.cases(sample)]
    print(paste("There are ", length(sample), " complete records in ", colnames[i], ".\n There are ",sum(unname(zero_col)), "zero records to fill"))
    #OutVals = boxplot(sample)$out
    #sample <- sample[!sample %in% OutVals]
    print(paste(length(sample), "records are used for imputing zero values in ", colnames[i]))
    mu <- mean(sample, na.rm=T)
    sd <- sd(sample, na.rm=T)
    data[zero_col,colnames[i]] <- perseus_zeros_impute(zero_col, mu,sd)
  }
  return(data)
}


###################################################################################
### Anova functions
##################################################################################
#### Volcano plot

draw_volcanoplot <- function(data, condition){
  case <- data[[condition]]$Case_non_zero
  control <- data[[condition]]$Control_non_zero
  outliers <- case>1
  ds <- data[[condition]]
  ds[["outliers"]] <- outliers
  ds <- ds[complete.cases(ds),]
  ds$uniprotID <- unlist(lapply(ds$uniprotID, function(x) strsplit(x,";")[[1]][1]))
  ds$uniprotID <- unlist(lapply(ds$uniprotID, function(x) strsplit(x,"-")[[1]][1]))
  fold_cutoff = 1
  pvalue_cutoff = 0.05
  ds$ID  <- apply(ds,1, function(x) strsplit(x[["Gene.names"]], ";")[[1]][1])
  ds$ID <- ifelse(is.na(ds$ID),ds$uniprotID,ds$ID)
  ds$ID[ds$uniprotID=="P11260"] <- "L1RE1"
  ds$ID[ds$uniprotID=="P11369"] <- "ORF2"
  # L1RE : P11260
  # Pol : P11369
  cold <- ifelse(ds$ID=="ORF2","purple", ifelse(ds$ID=="L1RE1", "red", ifelse(ds$outliers&(ds$Significant == "Yes")&(ds$Physical_evidence>0),"black", "grey")))
  group <- ifelse(ds$ID=="ORF2","ORF2", ifelse(ds$ID=="L1RE1", "L1RE1",ifelse(ds$outliers&ds$Significant == "Yes"&(ds$Physical_evidence>0), "Significant","Not significant/Outliers")))
  
 
  ds$colour <- ifelse(ds$ID=="ORF2","purple", ifelse(ds$ID=="L1RE1", "red",ifelse((ds$Significant == "Yes"& ds$outliers)&(ds$Physical_evidence>0),"black","grey")))
  lab <- ifelse((ds$Significant=="Yes"& ds$outliers)|ds$ID=="L1RE1",ds$ID,"")
  cols <- c("ORF2"="purple", "L1RE1"="red",  "Significant"="black",  "Not significant/Outliers"="grey")
  
  subds <-  subset(ds, ((ds$outliers)&Significant=="Yes")|(ID=="L1RE1")|(ID=="ORF2"))
  
  #pdf(paste0("../image/volcano_plot/Volcano_plot_",condition,".pdf"),width = 12, height = 18)
  png(paste0("../image/volcano_plot/Volcano_plot_",condition,".png"),width = 1000, height = 700)
  
  p <- ggplot(ds,aes(logfold, -log10(p.adj),label=ID)) +
    geom_point(aes(colour=group),fill = cold, size=2) +
    geom_vline(xintercept = fold_cutoff, col = "blue")+
    geom_vline(xintercept = -fold_cutoff, col = "blue")+
    geom_hline(yintercept = -log10(pvalue_cutoff), col = "green")+
    ggtitle(condition)+theme_minimal()+
    scale_colour_manual(values=cols, aesthetics = c("fill","colour"))+
    geom_text_repel(data  = subds, 
                    colour = ifelse(subds$ID =="ORF2","purple",ifelse(subds$ID=="L1RE1","red","black")), segment.size  = 0.2, 
                    segment.alpha =0.35,
                    box.padding = unit(0.45, "lines"), 
                    point.padding = unit(0.45, "lines"),
                    size=3 )
  # scale_x_continuous(limits = c(-7, 12))+
  #  scale_y_continuous(limits = c(0, 4.5))
  
  print(p)
  dev.off()
}


