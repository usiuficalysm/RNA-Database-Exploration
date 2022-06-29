library(tidyverse)
library(data.table)
library(stringi)

columns <- c('chr', 'source', 'feature', 'start',
              'end', 'score', 'strand',
              'phase', 'attributes')
genome <- fread("gencode.v40.annotation.gtf")
setnames(genome, colnames(genome), columns)


get_gene_df <- function(gtf_file = genome, gene_filter){
  dt_gene <- genome[grepl(gene_filter, genome$attributes)]
  test_strings <- dt_gene$attributes %>% str_split(";")
  
  column_names <- c()
  for (i in 1:length(test_strings)){
    pre <- unlist(test_strings[i])
    for (j in 1:length(pre)){
      pre_str <- pre[j] %>% str_split('/"') %>% unlist() %>% stri_remove_empty() %>% str_trim()
      if(length(pre_str) == 2){
        column_names <- column_names %>% append(pre_str[1]) %>% unique()
      }
    }
  }
  
  my_df <- data.frame(matrix(nrow = 1, ncol = length(column_names)))
  colnames(my_df) <- column_names
  
  for (i in 1:length(test_strings)){
    pre <- unlist(test_strings[i])
    tags <- ""
    for (j in 1:length(pre)){
      pre_str <- pre[j] %>% str_split('/"') %>% unlist() %>% stri_remove_empty() %>% str_trim()
      if(length(pre_str) == 2){
        if(pre_str[1] == 'tag'){
          tags <- tags %>% paste(pre_str[2], sep = "; ")
        }else{
          my_df[i, eval(pre_str[1])] <- pre_str[2]
        }
      }
    }
    tags <- tags %>% substring(3)
    my_df$tag[i] <- tags
    
  }
  dt_gene <- cbind(dt_gene, my_df)
  
  return(dt_gene)
}

samd <- get_gene_df(genome, gene_filter = 'SAMD')
linc <- get_gene_df(genome, gene_filter = 'LINC')

## ====== No more trash DBs ======

library(readxl)
basepath <- "C:/Users/Malcolm/OneDrive - University of Cambridge/Desktop/## Cambridge Academics/NST IB/Computational Things/RNA Screening/20220627 RNA Information"
file_list <- dir(basepath)

### ====  Parsing Disease_association as db2. ====
db2 <- fread(paste0(basepath, "/", file_list[4]))

# Extracts RNAs based on association score threshold and disease name matching.
extract_rnas <- function(df, disease_name, disease_colname = 'Disease name', threshold = 0.50){
  
  # Renaming - Adds anchor and end for exact match.
  disease_colname <- paste0('^', disease_colname, '$')
  names(df)[which(grepl(disease_colname, names(df)))] <- 'Disease name'
  
  if(!('Association score' %in% colnames(df))){
    return(df[tolower(`Disease name`) %like% tolower(disease_name)])
  }else{
    return(df[tolower(`Disease name`) %like% tolower(disease_name)
              & `Association score` >= threshold])
  }
}
test_infd <- extract_rnas(db2, 'Infectious disease')

### ====  Parsing Var_Disease_association as db3 and all_ncRNA as db1. ====

db3 <- fread(paste0(basepath, "/", file_list[5]))
test_immune <- extract_rnas(db3, 'immune', threshold = 0.2)
test_cancer <- extract_rnas(db3, 'cancer', threshold = 0.5)

db1 <- fread(paste0(basepath, "/", file_list[1]))
test_cancer <- extract_rnas(db1, disease_name = 'cancer', disease_colname = 'Disease Name')

### ==== Sorting by Disease Type - Generating Graphs ====

gene_types <- unique(db1$`Disease Name`)

for (type in gene_types){
  ## Name Cleaning
  type <- type %>% str_replace_all('/', ' or ')
  
  ## Filter by Disease types.
  first_filter <- db1[`Disease Name` == eval(type)]
  
  #### Operation 1: Graph Dysfunction Patterns.
  patterns <- unique(first_filter$`Dysfunction Pattern`)
  le_count <- c()
  for (i in patterns){
    le_count <- le_count %>% append(sum(first_filter$`Dysfunction Pattern` == i))
  }
  
  ## Calculate N/A % and remove N/A values - who needs them anyways?
  first_filter_dyspattern <- data.table(`Pattern Type` = patterns, 'Frequency' = le_count)
  na_percentage <- round(sum(first_filter$`Dysfunction Pattern` == 'N/A') / length(first_filter$`Dysfunction Pattern`), 3)
  first_filter_dyspattern <- first_filter_dyspattern[`Pattern Type` != 'N/A']
  
  ## All N/A catch; Commence Plotting
  if (dim(first_filter_dyspattern)[1] == 0){
    test_plot <- ggplot() + annotate("text", x = 10, y = 10, label = "All N/A Values.", size = 18) + 
          theme_void() + ggtitle(paste0("Pattern Types - ", type)) +
           theme(plot.title = element_text(size = 24, face = 'bold', hjust = 0.5))
  }else{
    test_plot <- ggplot(first_filter_dyspattern, aes(x = `Pattern Type`, y = Frequency, fill = `Pattern Type`)) +
      geom_bar(width = 0.3, stat = 'identity') +
      scale_fill_brewer(palette = 'RdPu') + 
      # geom_text(x = 0.2*length(first_filter_dyspattern$`Pattern Type`),
      #           y = 0.5*max(first_filter_dyspattern$Frequency),
      #           label = paste0('N/A Proportion = ', na_percentage),
      #           size = 5) +
      ggtitle(paste0("Pattern Types - ", type, "\n (N/A Proportion = ", na_percentage, ')')) +
      theme_minimal() +
      theme(legend.position = 'none',
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 14, face = 'bold'),
            plot.title = element_text(size = 24, face = 'bold', hjust = 0.5)
      ) + 
      coord_flip()
  }
  
  ## Commence Printing.
  mywd <- getwd()
  dir.create(file.path(mywd, "Pattern Types and Frequencies"), showWarnings = F)
  setwd(file.path(mywd, "Pattern Types and Frequencies"))
  
  jpeg(filename = paste0(Sys.Date(), '_', type, '.jpeg') ,width = 11, height = 8, units = 'in', res = 300)
  print(test_plot)
  print(type)
  dev.off()
  
  setwd('..')
}
