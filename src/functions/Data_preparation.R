## Mehrnoosh Oghbaie
## 07/24/2019
## Preparing data from MaxQuant or ...

## Data preparation consists of four stages:
##  1. Remove contaminants and reverse proteins
##  2. Log transformation
##  3. Separate different experiment


Template <- R6Class("Template",
                    private = list(
                      input.info.dir = "../Input.info"
                    ), 
                    public = list(df = NA,
                                  df_log=NA),
                    active = list(
                      input_dir = function() {
                        base::dirname(private$input.info.dir)
                      }, 
                      input = function(){
                        read.delim(private$input.info.dir, header = FALSE)
                      },
                      data_type = function(){
                        self$input %>%
                          dplyr::filter(V1=="data_type") %>%
                          dplyr::mutate(V2=as.character(V2))%>% 
                          dplyr::select(V2) %>% 
                          .$V2
                      }, 
                      
                      input_merged = function(){
                        file_name <- self$input %>%
                          dplyr::filter(V1=="MaxQuant_Output") %>%
                          dplyr::mutate(V2=as.character(V2))%>% 
                          dplyr::select(V2) %>% 
                          .$V2
                        
                        tbl <- read.delim(file_name)
                        tbl$uniprotID <- apply(tbl, 1, function(x) strsplit(as.character(x[["Protein.IDs"]]), ";")[[1]][1])
                        
                        
                        return(tbl)
                      }
                      ,
                      contaminant_list =function(){
                        contaminants_tbl <- read.table(paste(self$input_dir,
                                                             self$input %>%
                                                               dplyr::filter(V1=="Contaminants") %>%
                                                               dplyr::mutate(V2=as.character(V2))%>% 
                                                               dplyr::select(V2) %>% 
                                                               .$V2, sep="/"),
                                                       sep=";", quote="\"")
                        contaminant_list <- unlist(lapply(as.character(contaminants_tbl[,1])[grepl(">", as.character(contaminants_tbl[,1]))], function(x) strsplit(x," |>")[[1]][2]))
                        return(contaminant_list)
                      } 
                    )
)

Template$set("public","non_unique_ORF", list())


Template$set("public","removeContaminant", function(){
  cols <- colnames(self$input_merged)
  dl <-  self[["input_merged"]]
  dl <- dl %>% dplyr::mutate(Peptide.counts..unique. = as.character(Peptide.counts..unique.),
                             Protein.IDs = as.character(Protein.IDs),
                             Potential.contaminan = as.character(Potential.contaminant))
  dl$Peptide.counts..unique. <-  apply(dl,1,function(x) strsplit(x[["Peptide.counts..unique."]],";")[[1]][1])
  self[["non_unique_ORF"]] <- dl  %>% dplyr::filter(Peptide.counts..unique.=="0" & grepl("ORF1",Protein.IDs)) %>% .$uniprotID
  self$df <- dl[!(dl$Potential.contaminant=="+"|dl$Reverse=="+") ,c("id","uniprotID","Gene.names",cols[grepl(self$data_type,cols)])]    
})


Template$set("public","logTransformation", function(){
  
  self$df[,colnames(self$df)[grepl(self$data_type,colnames(self$df))]] <- apply(self$df[,colnames(self$df)[grepl(self$data_type,colnames(self$df))]],2, as.numeric)
  tp <- lapply(self$df, class)
  df_log <- self$df
  for(i in colnames(self$df)[tp =="numeric"]){
    df_log[i] <- log(unlist(unname(df_log[i])))
    df_log[is.infinite(unname(unlist(df_log[i]))),i] <- NA
  }
  self$df_log <- df_log
  
}
)

Template$set("public","experiment", list())

Template$set("public","separatedList", function(run_order){
  selected_col <- colnames(self$df_log)
  
  for(i in 1:dim(run_order)[1]){
    self$experiment[[paste(run_order[i,], collapse=".vs.")]] <-  self$df_log[,c("id","uniprotID","Gene.names",selected_col[grepl(paste(paste0(self$data_type,".",run_order[i,]),collapse="|"),selected_col)])]
  }
}
)




