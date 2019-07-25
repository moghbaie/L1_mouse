## Mehrnoosh Oghbaie
## 06/21/2019
## Integrate normalized Values

## Normalization consists of two stages:
##  1. Normalize intensities by LORF1 intensity
##  2. Take average of normalized intensity on all cases

Template$set("public","experimentAveragedNormalized", list())
Template$set("public","experimentAveraged", list())
Template$set("public","calculateAverageNormalized", function(run_order){
  average_list <- data.frame(matrix(NA, ncol=0,nrow=dim(self[["input_merged"]])[1]))
  average_list$id <- self[["input_merged"]]$id
  average_list$uniprotID <- unlist(lapply(self[["input_merged"]]$uniprotID,function(x) strsplit(x, ";")[[1]][1]))
  average_list$uniprotID <- unlist(lapply(self[["input_merged"]]$uniprotID,function(x) strsplit(x, "-")[[1]][1]))
  average_list$Gene.names <- unlist(lapply(as.character(self[["input_merged"]]$Gene.names),function(x) strsplit(x, ";")[[1]][1]))
  average_list$Potential.contaminant <- as.character(self[["input_merged"]]$Potential.contaminant)
  average_list$Reverse <- as.character(self[["input_merged"]]$Reverse)
  
  average_normalized_list <- average_list
  
  for(name in  unique(c(run_order[,1],run_order[,2]))){
    if(!grepl("ORF1",name)){
      norm_int <- sum(apply(self$input_merged[grepl("ORF2",self$input_merged$uniprotID)& self$input_merged$uniprotID %in% self[["Significant_list"]]$uniprotID,colnames(self$input_merged)[grepl(paste0("LFQ.intensity.",name), colnames(self$input_merged))]],1, function(x) mean(ifelse(x==0,0, log(x)))))
      if(norm_int==0){
        norm_int <- max(apply(self$input_merged[,colnames(self$input_merged)[grepl(paste0("LFQ.intensity.",name), colnames(self$input_merged))]],1, function(x) mean(ifelse(x==0,0, log(x)))))
        
      }
    } else {
      norm_int <- sum(apply(self$input_merged[(grepl("ORF1",self$input_merged$uniprotID)& self$input_merged$uniprotID %in% self[["Significant_list"]]$uniprotID)|self$input_merged$uniprotID=="Q9UN81",
                                              colnames(self$input_merged)[grepl(paste0("LFQ.intensity.",name), colnames(self$input_merged))]],1, function(x) mean(ifelse(x==0,0, log(x)))))
      
    }
    print(name)
    print(norm_int)
    average_normalized_list[[name]] <- apply(self$input_merged[,colnames(self$input_merged)[grepl(paste0("LFQ.intensity.",name), colnames(self$input_merged))]],1, function(x) mean(ifelse(x==0,0, log(x)))/norm_int)
  }
  self[["experimentAveragedNormalized"]] <- average_normalized_list
}
)

####################################################################################################################
### Draw venn diagram
Template$set("public","drawCommonVenndiagram", function(){
  venn_list <- self[["Significant_proteins"]][-(1:3)]
  venn_list[venn_list!=0] <- 1
  png(filename="../Image/Integrated_plot/VennDiagram.png",width = 1000, height = 800)
  a <- vennCounts(venn_list)
  vennDiagram(a,circle.col = c("darkmagenta", "darkblue",  "orange"), main="Comparison between expressed proteins")
  dev.off()
})

###################################################################################################################
### Draw heatmap of average value for significant ones
### Draw heatmap of average value for significant ones
Template$set("public","drawHeatmap", function(){
  
  da <- melt(self$Significant_proteins, id.vars=c( "id", "uniprotID","Gene.names"))
  db <- melt(self[["experimentAveragedNormalized"]], id.vars=c( "id", "uniprotID","Gene.names"))
  #db <- melt(average_unnormalized_list, id.vars=c( "id", "uniprotID","Gene.names"))
  dl <- cbind(da,db[-c(1:4)])
  dl$Gene.names <- unlist(lapply(as.character(dl$Gene.names),function(x) strsplit(x, ";")[[1]][1]))
  dl$Gene.names <- ifelse(is.na(dl$Gene.names),dl$uniprotID,dl$Gene.names)
  
  dl2 <- dl[(dl$uniprotID %in% self[["Significant_list"]]$uniprotID)|grepl("ORF",dl$uniprotID),]
  a <-ifelse(unique(dl2[["Gene.names"]])[order(unique(dl2[["Gene.names"]]))] %in% unique(eLife_list$geneName),"red","black")
  sequence_length = length(unique(dl2$uniprotID))
  second_sequence = c((sequence_length%/%2+1):sequence_length)
  first_sequence = c(1:(sequence_length%/%2)) 
  first_angles = c(90 - 180/length(first_sequence) * first_sequence)
  second_angles = c(-90 - 180/length(second_sequence) * second_sequence)
  colnames(dl2) <- c("id", "uniprotID", "Gene.names", "condition","Significance","Expression")
  dl2$Significance <- ifelse((grepl("ORF",dl2$uniprotID)&dl2$Expression>0)|dl2$Significance==1,1,0) 
  ## polar plot with only colon significant
  png(filename="../Image/Integrated_plot/heatmap_significant_eLife.png",width = 1200, height = 1200)
  q <- ggplot(dl2, aes(condition, Gene.names)) + 
    geom_tile(aes(fill = Expression, alpha=Significance),color="gray35")+
    scale_alpha(range = c(0.5, 1))+
    geom_tile(alpha=0.6-dl2$Significance*(3/5), fill="gray")+
    #scale_fill_gradient( limits=c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),  low = "blue", high = "red3")+
    scale_fill_gradient2(mid = "green",low = "black",  high = "red3")+
    scale_y_discrete(position = "right")+
    coord_polar(theta = "y")+
    theme_bw()+
    theme(axis.text.x = element_text(angle= c(first_angles,second_angles),size=10, colour = a),
          axis.text.y = element_text(size=14),
          axis.title=element_text(size=14,face="bold"))
  
  print(q)
  dev.off()
  
  png(filename="../Image/Integrated_plot/heatmap_significant_all.png",width = 600, height = 1600)
  qq <- ggplot(dl2, aes(condition, Gene.names)) + 
    geom_tile(aes(fill = Expression, alpha=Significance),color="gray35")+
    scale_alpha(range = c(0.5, 1))+
    geom_tile(alpha=0.6-dl2$Significance*(3/5), fill="gray")+
    #scale_fill_gradient( limits=c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),  low = "blue", high = "red3")+
    scale_fill_gradient2(mid = "green",low = "black",  high = "red3")+
    scale_y_discrete(position = "right")+
    #coord_polar(theta = "y")+
    theme_bw()+
    theme(
      axis.text.x = element_text(
        #angle= c(first_angles,second_angles),
        angle =90, size=10, colour = a),
      axis.text.y = element_text(size=12),
      axis.title=element_text(size=14,face="bold"))
  
  print(qq)
  dev.off()
  
})



###################################################################################################################
### MDS plot
Template$set("public","drawMDSplot", function(){
  average_list <- self$experimentAveragedNormalized
  dz <- average_list%>%filter(uniprotID%in% self$Significant_list$uniprotID)
  rownames(dz) <- dz$id
  d4 = dist(dz[-(1:2)])
  
  fit4 = cmdscale(d4, eig=TRUE, k=3) #k is number of dimensions
  dim1 = fit4$points[,1]
  dim2 = fit4$points[,2]
  dim3 = fit4$points[,3]
  
  colfunc <- colorRampPalette(c("green", "red3"))
  color <- data.frame(seq(10)/10,colfunc(10))
  color <- color$colfunc.10.[match(round(dz$Average.163T_ORF1_6,1),color$seq.10..10)]
  
  gene_protein <- read.delim("protein.tab")
  gene_protein <- gene_protein %>% dplyr::mutate(Gene.names = as.character(Gene.names),
                                                 Entry = as.character(Entry))
  gene_protein$Gene.names  <- apply(gene_protein,1, function(x) strsplit(x[["Gene.names"]], " ")[[1]][1])
  
  name <- gene_protein$Gene.names[match(dz$uniprotID, gene_protein$Entry)]
  
  dz[["GeneName"]] <- ifelse(!is.na(name), name , dz$uniprotID)
  dz[["GeneName"]][dz[["GeneName"]]=="Q9UN81"] <- "L1RE1"
  plot3d(dim1, dim2, dim3, type="s", size =1, lwd=4,col=color)
  text3d(dim1, dim2, dim3+0.08, ifelse(!is.na(color),dz[["GeneName"]],""), fontweight="bold", cex= 0.7, col="black")
  #create a spinning object
  s = spin3d(axis = c(0,0,1), rpm=2)
  #play 3d
  play3d(s, duration=33)
  ### You need to have Imagemagick installed and give the path to convert.exe to it
  imconvertstring<-"\"c:\\Program Files\\ImageMagick-7.0.8-Q16\\convert.exe\" -delay 1x%d %s*.png %s.%s"
  movie3d(s, duration=33, dir="../Image/Integrated_plot/", clean=T, convert = imconvertstring, type = "gif")
  
})
#######################################################################

