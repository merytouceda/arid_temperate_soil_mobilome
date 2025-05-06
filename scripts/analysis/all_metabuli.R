# metabuli processing for all contigs of base project


library(tidyverse)
library(car)
library(lme4)
library(car)
library(performance)
library(tidyverse)


############
### join taxonomy to contig id
############

##
# 1. Create a table of contig to ID from the classifications
##
# load wd
setwd("/Volumes/BunnyBike/mge_urban/base/metabuli/metabuli/classifications/")

# initialize empty table
classifications_all <- data.frame(matrix(NA, nrow = 0, ncol = 5))  #after reviewing this I think the "wide" format is most useful since it keep the LOD information
colnames(classifications_all) <- c("classified_yesno", "contig", "tax_id" ,"notsure" ,"level")


files <- list.files(pattern="*.tsv", full.names=F, recursive=FALSE)

# Put them all together: 
for (file in files){
  
  file <- tryCatch(read_tsv(file), error=function(e) NULL)
  
  colnames(file) <- c("classified_yesno", "contig", "tax_id" ,"notsure" ,"level")
  
  classifications_all <- rbind(classifications_all, file)
  
}

classifications_all_final <-  classifications_all %>%
  dplyr::select(c("contig", "tax_id")) %>%
  filter(tax_id != 0)



##
# 2. Create a table with all the taxonomy and their ids
##
# load wd
setwd("/Volumes/BunnyBike/mge_urban/base/metabuli/metabuli/reports/")

# initialize empty table
reports_all <- data.frame(matrix(NA, nrow = 0, ncol = 6))  #after reviewing this I think the "wide" format is most useful since it keep the LOD information
colnames(reports_all) <- c("noidea1", "noidea2", "noidea3" ,"rank" ,"tax_id", "taxon")

files <- list.files(pattern="*.tsv", full.names=F, recursive=FALSE)

# Put them all together: 
for (file in files){
  
  file <- tryCatch(read_tsv(file), error=function(e) NULL)
  
  colnames(file) <- c("noidea1", "noidea2", "noidea3" ,"rank" ,"tax_id", "taxon")
  
  reports_all <- rbind(reports_all, file)
  
}

reports_all_final <- reports_all %>%
  dplyr::select(-c("noidea1", "noidea2", "noidea3" ,"rank")) %>%
  unique()


##
# 3. Join them
##
taxon_id_table <- classifications_all_final %>%
  left_join(reports_all_final, by = "tax_id")

taxon_distrib <- taxon_id_table %>%
  group_by(taxon) %>%
  tally() %>%
  arrange(-n) %>%
  slice_head(n = 20)

ggplot(taxon_distrib, aes(x= reorder(taxon, n), y = n))+
  geom_col()+
  xlab("Taxonomic classification")+
  ylab("Number of plasmids")+
  coord_flip()+
  theme_bw()
ggsave("/Volumes/BunnyBike/mge_urban/local/figures/mges_taxons_20.pdf", device = "pdf", width = 5, height = 6 , units = "in")


