library(readr)
library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

#Genes
meta_file <- read_tsv(file=args[1], col_names = T, skip_empty_rows = T)
sample_name_dict <- paste0(meta_file$sample_id, "_", meta_file$condition)
names(sample_name_dict) <- meta_file$sample_id

#dict
if(args[2] == "hs"){
  dict <- read_tsv("/groups/zuber/USERS/florian.andersch/data/transcripts/human/mRNAs_human2021.txt") 
}
if(args[2] == "mm"){
  dict <- read_tsv("/groups/zuber/USERS/florian.andersch/data/transcripts/mouse/mRNAs_mouse2021_mm10.txt")
}

dict$entrez_id <- as.character(dict$entrez_id)
dict_stripped <- dict %>%
  select(entrez_id, symbol) %>%
  distinct

setwd(args[3])

files <- list.files(pattern = "*_tcount_collapsed.csv")
final_file<-NULL
for(i in 1:length(files)){
  file_buff <- read_tsv(file=files[i], col_names = T)
  file_buff$nonTcReadCount <- file_buff$readCount - file_buff$tcReadCount
  file_buff$RPMu <- (file_buff$tcReadCount / sum(file_buff$nonTcReadCount)) * 10^6
  file_buff$gene_name <- as.character(file_buff$gene_name)
  file_buff <- file_buff %>%
    left_join(dict_stripped, by=c("gene_name"="entrez_id")) %>%
    select(entrez=gene_name, starts_with("symbol"), starts_with("entrez"), length, everything()) %>%
    arrange(symbol)
  
  sample_id <- files[i] %>% str_split(pattern="_") %>% unlist 
  sample_id<-sample_id[1]
  print(sample_id)
  
  sample_name <- sample_name_dict[sample_id]
  
  if(!is.na(sample_name)){
    colnames(file_buff)<-c(colnames(file_buff)[1:3], paste0(sample_name, "_", colnames(file_buff)[4:length(colnames(file_buff))]))
    
    
    
    if(is.null(final_file)){
      final_file <- file_buff
      
      final_counts <- file_buff %>%
        select(entrez, symbol, length, ends_with("_readCount"))
      colnames(final_counts)[4]<-sample_name
      
      final_tcCount <- file_buff %>%
        select(entrez, symbol, length, ends_with("_tcReadCount"))
      colnames(final_tcCount)[4]<-sample_name
      
      final_RPM <- file_buff %>%
        select(entrez, symbol, length, ends_with("_readsCPM"))
      colnames(final_RPM)[4]<-sample_name
      
      final_RPMu <- file_buff %>%
        select(entrez, symbol, length, ends_with("_RPMu"))
      colnames(final_RPMu)[4]<-sample_name
      
    }else{
      final_file <- final_file %>%
        full_join(file_buff)
      
      final_counts <- final_counts %>%
        full_join(file_buff %>%
                    select(entrez, symbol, length, ends_with("_readCount")))
      colnames(final_counts)[ncol(final_counts)]<-sample_name
      
      final_tcCount <- final_tcCount %>%
        full_join(file_buff %>%
                    select(entrez, symbol, length, ends_with("_tcReadCount")))
      colnames(final_tcCount)[ncol(final_tcCount)]<-sample_name
      
      final_RPM <- final_RPM %>%
        full_join(file_buff %>%
                    select(entrez, symbol, length, ends_with("_readsCPM")))
      colnames(final_RPM)[ncol(final_RPM)]<-sample_name
      
      final_RPMu <- final_RPMu %>%
        full_join(file_buff %>%
                    select(entrez, symbol, length, ends_with("_RPMu")))
      colnames(final_RPMu)[ncol(final_RPMu)]<-sample_name
      
    }
  }
}


write_tsv(final_file, file="merged_files_all_samples.txt", col_names = T)
write_tsv(final_counts, file="counts.txt", col_names = T)
write_tsv(final_tcCount, file="tcCounts.txt", col_names = T)
write_tsv(final_RPM, file="RPMs.txt", col_names = T)
write_tsv(final_RPMu, file="RPMus.txt", col_names = T)



