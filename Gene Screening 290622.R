# library(tidyverse) # Cannot install for some reason.
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)
library(stringi)
library(tidyr)

### ==== (Ignore this) Was attempting to parse attribute strings in .gtf file ====

# columns <- c('chr', 'source', 'feature', 'start',
#               'end', 'score', 'strand',
#               'phase', 'attributes')
# genome <- fread("gencode.v40.annotation.gtf")
# setnames(genome, colnames(genome), columns)
# 
# 
# get_gene_df <- function(gtf_file = genome, gene_filter){
#   dt_gene <- genome[grepl(gene_filter, genome$attributes)]
#   test_strings <- dt_gene$attributes %>% str_split(";")
#   
#   column_names <- c()
#   for (i in 1:length(test_strings)){
#     pre <- unlist(test_strings[i])
#     for (j in 1:length(pre)){
#       pre_str <- pre[j] %>% str_split('/"') %>% unlist() %>% stri_remove_empty() %>% str_trim()
#       if(length(pre_str) == 2){
#         column_names <- column_names %>% append(pre_str[1]) %>% unique()
#       }
#     }
#   }
#   
#   my_df <- data.frame(matrix(nrow = 1, ncol = length(column_names)))
#   colnames(my_df) <- column_names
#   
#   for (i in 1:length(test_strings)){
#     pre <- unlist(test_strings[i])
#     tags <- ""
#     for (j in 1:length(pre)){
#       pre_str <- pre[j] %>% str_split('/"') %>% unlist() %>% stri_remove_empty() %>% str_trim()
#       if(length(pre_str) == 2){
#         if(pre_str[1] == 'tag'){
#           tags <- tags %>% paste(pre_str[2], sep = "; ")
#         }else{
#           my_df[i, eval(pre_str[1])] <- pre_str[2]
#         }
#       }
#     }
#     tags <- tags %>% substring(3)
#     my_df$tag[i] <- tags
#     
#   }
#   dt_gene <- cbind(dt_gene, my_df)
#   
#   return(dt_gene)
# }
# 
# samd <- get_gene_df(genome, gene_filter = 'SAMD')
# linc <- get_gene_df(genome, gene_filter = 'LINC')

## ====== No more weird DBs ======

library(readxl)
basepath <- file.path(getwd(), 'RNA Databases')
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

### ==== Sorting by Disease Type - Generating Graphs (Pattern Types) ====

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

rm(first_filter, first_filter_dyspattern, test_plot, gene_types, i, le_count, mywd,
   na_percentage, patterns, type)

#### ==== Tackling db2/3: Association Score ====
## Possible to set a association score filter
### db2
extract_rna_assoc <- function(df, threshold = 0.5){
  return(df[`Association score` >= 0.5])
}
db2_essential <- extract_rna_assoc(db2)
db2_assoc_percentage <- round(length(db2_essential$`Disease name`) / length(db2$`Disease name`) * 100, 2)
print(paste0('% of diseases with at least 1 RNA with >=0.5 association score: ', db2_assoc_percentage, '.'))
# We will need to supplement diseases without sufficiently associated RNA.
hist(db2_essential$`Association score`, breaks = 50)

### db3
names(db3)
db3_essential <- extract_rna_assoc(db3)
# % of diseases with sufficiently associated RNA.
print(paste0('% of diseases with at least 1 RNA with >=0.5 association score: ', round((nrow(db3_essential) / nrow(db3) * 100), 1), '.'))
hist(db3_essential$`Association score`, breaks = 50)

### ==== Heatmapping the diseases (db1) with Regulatory Filters ====
db1_clean <- db1[`Dysfunction Pattern` != 'N/A']
regulation_cat <- db1_clean$`Dysfunction Pattern` %>% str_remove(' \\[.*$')
regulation_cat[which(grepl('up-regulated|down-regulated', db1_clean$`Dysfunction Pattern`))] <- db1_clean$`Dysfunction Pattern`[which(grepl('up-regulated|down-regulated', db1_clean$`Dysfunction Pattern`))]
db1_clean[, `Regulation Category` := regulation_cat]

## Further data preparation for heatmapping.
combined_test <- paste(db1_clean$`Disease Name`, db1_clean$`Regulation Category`, sep = '+')
heatmappable <- as.data.table(table(combined_test))
heatmappable <- heatmappable[,`:=` (`Disease Name` = str_to_title(str_remove(heatmappable$combined_test, '\\+.*$')),
                                    DysPattern = str_remove(heatmappable$combined_test, '^.*\\+'),
                                    LogN = log(heatmappable$N))]

## Heatmapping.
disease_dys_heatmap <- ggplot(heatmappable, aes(y = `Disease Name`, x = DysPattern)) +
                        geom_tile(aes(fill = LogN), colour = 'white') +
                        scale_fill_gradient(low = 'white', high = 'red') +
                        ggtitle("Heatmap of Diseases and Dysfunction Patterns") + xlab("Dysfunction Pattern") +
                        theme_minimal() + theme(axis.title.x = element_text(size = 12),
                                                axis.title.y = element_text(size = 12),
                                                axis.text.x = element_text(size = 8, angle = 45,
                                                                           hjust = 1, vjust = 1),
                                                axis.text.y = element_text(size = 8),
                                                title = element_text(size = 16),
                                                legend.key.size = unit(12, 'points'),
                                                legend.title = element_text(size = 12))

jpeg(filename = 'DiseaseDysHeatmap.jpeg', width = 11, height = 66, units = 'in', res = 300)
print(disease_dys_heatmap)
dev.off()

### ==== Creating BarChart of Frequency of RNA to Dyspattern ====
top_ten <- heatmappable[order(heatmappable$N, decreasing = T),]
top_ten <- top_ten[grepl("Regulation \\[", DysPattern)]

plot_bar_regulation <- function(df, keyword = 'up', color = '#70a8e0'){
  reg_bar <- ggplot(top_ten[grepl(paste0("Regulation \\[", keyword, "-regulated\\]"), DysPattern)][1:20,]) +
    geom_bar(aes(y = reorder(`Disease Name`, N), x = LogN), stat = 'identity', width = 0.5, fill = color) +
    theme_minimal() +
    ggtitle(paste0("Log(nRNA) of Top 20 Disease Types for ", stri_trans_totitle(keyword), "-Regulation")) +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size= 12),
          axis.title.y = element_text(size= 12),
          plot.title = element_text(size = 16, hjust = 0.5))
  return(reg_bar)
}

upreg_bar <- plot_bar_regulation(top_ten)
downreg_bar <- plot_bar_regulation(top_ten, keyword = "down", color = "#df5454")

jpeg(filename = 'Top20-Upregulated.jpeg', width = 11, height = 8, units = 'in', res = 300)
print(upreg_bar)
dev.off()

jpeg(filename = 'Top20-Downregulated.jpeg', width = 11, height = 8, units = 'in', res = 300)
print(downreg_bar)
dev.off()

### ==== Finding Name Linkages ====
db2_sims <- db2_essential[db2_essential$ncRNA %in% toupper(db1_clean$`ncRNA Symbol`)]
db2_sims <- db2_essential[db2_essential$ncRNA %in% db3_essential$ncRNA]

db1_sims <- db1_clean[toupper(db1_clean$`ncRNA Symbol`) %in% db2_essential$ncRNA]
db1_sims <- db1_clean[toupper(db1_clean$`ncRNA Symbol`) %in% db3_essential$ncRNA]

commonRNAlist <- c(toupper(unique(db1_sims$`ncRNA Symbol`)), db2_sims$ncRNA, db3_sims$ncRNA) %>% unique()

col_names <- unique(c(names(filter_1[c(1, 4:8)]), names(filter_2[c(5,3,8)]), names(filter_3[c(2,6,7,4)])))
compiled_rnas <- data.frame(matrix(nrow = 0, ncol = length(col_names)))
colnames(compiled_rnas) <- col_names

for (RNA in commonRNAlist){
  filter_1 <- db1_sims[`ncRNA Symbol` == RNA][,c(1,4:8)]
  names(filter_1)[1] <- 'ncRNA'
  filter_2 <- db2_sims[ncRNA == RNA][,c(5,3,8)]
  filter_3 <- db3_sims[ncRNA == RNA][,c(2,6,7,4)]
  huge_mess <- rbind(filter_1, filter_2, filter_3, fill = T)
  compiled_rnas <- compiled_rnas %>% rbind(huge_mess)
  
  rm(filter_1, filter_2, filter_3, huge_mess, RNA)
}

write.csv(compiled_rnas, "20220713 commonRNAs.csv", row.names = T, na = '')

print(paste0("Total unique common and significant RNAs across all 3 databases: ", length(unique(compiled_rnas$ncRNA))))
