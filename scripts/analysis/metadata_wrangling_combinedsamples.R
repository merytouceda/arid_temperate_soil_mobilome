# metadata base wrangling for 549 samples


library(tidyverse)
library(ggbiplot)
library(cowplot)

md1 <- read_csv("/Volumes/BunnyBike/mge_urban/base/package_metadata_bpa_23b1dae2_20241007T2157_amdb-metagenomics-novaseq.csv")
md2 <- read_csv("/Volumes/BunnyBike/mge_urban/base/package_metadata_bpa_23b1dae2_20241007T2157_base-metagenomics.csv")

# remove columns that are all NAs
md1_clean <- md1 %>%
  dplyr::select(where(~ !all(is.na(.)))) %>%
  separate(title, into = c("title", "sample_id"), sep = "\\/")

md2_clean <- md2 %>%
  dplyr::select(where(~ !all(is.na(.)))) %>%
  separate(title, into = c("bla", "bla2", "sample_id"), sep = " ")# %>%
  #separate(title, into = c("sample_id", "bla3"), sep = "_")


# Make sure the samples are the ones I am processing

# load list of samples being processed
samples_combined <- read.table("/Volumes/BunnyBike/mge_urban/base/samples_combined_list.txt", header = F)
colnames(samples_combined) <- c("sample_id")

samples_combined <- samples_combined %>%
  mutate(sample_id = gsub(sample_id, pattern = "_combined", replace = ""))
         

samples_in_md2 <- samples_combined %>%
  filter(sample_id %in% md2_clean$sample_id)


samples_not_in_md2 <- samples_combined %>%
  filter(!sample_id %in% md2_clean$sample_id)


samples_in_md1 <- samples_combined %>%
  filter(sample_id %in% md1_clean$sample_id)

samples_in_md1$sample_id == samples_not_in_md2$sample_id


# check the vegetation type and environmet of the samples I am processing


md2_clean_hiseq <- md2_clean %>%
  filter(!str_detect(notes, pattern = "18S"))

md2_clean$vegetation_type

ggplot(md1_clean, aes(x = general_env_feature )) + 
  geom_bar() +
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(md1_clean, aes(x = vegetation_type )) + 
  geom_bar() +
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(md2_clean_hiseq, aes(x = vegetation_type )) + 
  geom_bar() +
  facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))


md2_clean_hiseq_temparid <- md2_clean_hiseq %>%
  filter(general_env_feature %in% c("Arid", "Temperate")) %>%
  filter(vegetation_type %in% c("Forest", "Grassland", "Shrubland", "Woodland"))

ggplot(md2_clean_hiseq_temparid, aes(x = vegetation_type )) + 
  geom_bar(stat = "count") +
  stat_count(geom = "text", colour = "white", size = 3.5,
             aes(label = ..count..),position=position_stack(vjust=0.5)) +
  facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))

md2_clean_hiseq_temparid_samples <- md2_clean_hiseq_temparid %>%
  separate(sample_id, into = c("sample_id", "nada"), sep = "_") %>%
  mutate(sample_id = paste(sample_id, "combined", sep = "_")) %>%
  select(sample_id) %>%
  unique()

# save md
#md_save <-  md2_clean_hiseq_temparid %>%
#  distinct() %>%
#  mutate(Sample = gsub(sample_id, pattern = "_[0-9]*$", replacement = ""))
#write_csv(md_save , "/Volumes/BunnyBike/mge_urban/base/md_base_march25.csv")

# save sample names to send to HPC for filtering samples for processing
#write.csv(md2_clean_hiseq_temparid_samples, "/Volumes/BunnyBike/mge_urban/base/samples_to_process_120524.csv", quote = F)

moved_files <- read.table("/Volumes/BunnyBike/mge_urban/base/moved_files.txt")

missing_files <- md2_clean_hiseq_temparid_samples %>%
  mutate(sample_id = gsub(sample_id, pattern = "_combined", replacement = "")) %>%
  filter(!sample_id %in% moved_files$V1)





# ---------------------------------------------------------------------------------------------------------
# MD visualizations
# use the md final from ptu analyses, this md will only include the samples that are in the md and have been processed
# they should be 152

# testing the relationships between variables
md_variables <- md_final %>%
  select(sample_id, where(is.numeric)) %>%
  filter(!sample_id == "12891_1") %>% # this sample had a bunch of NA values, so it makes plotting very difficult
  select(-sample_id) %>%
  select(-c("depth", "depth_lower", "depth_upper", "slope_aspect", "slope_gradient", 
            "water_content", "sand", "silt", "fine_sand", "clay", "coarse_sand", 
            "ph_solid_h2o", "latitude", "longitude")) # these variables had a lot of NAs too
  


######
# Corplot
######
corrplot::corrplot(cor(md_variables, use = "complete.obs"), 
                   type = "lower", 
                   tl.col = "black", 
                   order = "hclust", 
                   diag = FALSE)


######
# PCA
######
# Standardize the data (PCA works best with standardized data)
md_variables_scaled <- scale(md_variables)

# run PCA
md_variables_result <- prcomp(md_variables_scaled , center = TRUE, scale. = TRUE)
# see summary results
summary(md_variables_result)

# make sure they are in the same order
all(rownames(md_final) == rownames(as.data.frame(md_variables_result$x)))
all(rownames(md_final) == rownames(as.data.frame(pca_scores)))

# Extract scores and add groupings
pca_scores <- as.data.frame(md_variables_result$x)
pca_scores$general_env_feature <-md_final$general_env_feature
pca_scores$vegetation_type <-md_final$vegetation_type

# Extract PCA loadings (arrows for variables)
pca_loadings <- as.data.frame(md_variables_result$rotation)
pca_loadings$Variable <- rownames(pca_loadings)  # Variable names

# Scaling the loadings for better visualization
arrow_scale <- 20 # Adjust scaling factor as needed
pca_loadings <- pca_loadings %>%
  mutate(PC1 = PC1 * arrow_scale, PC2 = PC2 * arrow_scale)


# Create PCA Biplot with ggplot2
ggplot() +
  # Plot PCA scores (points for samples)
  geom_point(data = pca_scores, aes(x = PC1, y = PC2, 
                                    color = vegetation_type, shape = general_env_feature),
             size = 4) +
  
  # Add PCA loadings (arrows for variables)
  geom_segment(data = pca_loadings, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black", linewidth = 1.2) +
  
  # Add labels to arrows
  geom_text(data = pca_loadings, aes(x = PC1, y = PC2, label = Variable), 
            vjust = 1, hjust = 1, size = 5, color = "black") +
  
  # Titles and themes
  labs(title =,x = "PC1", y = "PC2") +
  theme_bw() +
  scale_color_manual(values = pal1)+
  scale_shape_manual(values=c(16, 13))+
  theme(legend.position = "bottom")
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/pca_samples.pdf", device = "pdf", width = 7, height = 8 , units = "in")






######
# Variables loadings
######
md_variables_result$rotation

# Convert the loadings to a data frame for easier plotting
loadings_df <- as.data.frame(md_variables_result$rotation)

# Create a tidy version of the data (long format)
loadings_long <- loadings_df %>%
  rownames_to_column("Variable") %>%
  gather(PC, Loading, -Variable) %>%
  filter(PC %in% c("PC1", "PC2"))

# Sort variables by the absolute value of their loading for the first principal component (PC1)
loadings_long <- loadings_long %>%
  mutate(Variable = factor(Variable, 
                           levels = unique(loadings_long$Variable[order(abs(loadings_long$Loading), decreasing = TRUE)])))

# Plot the loadings for the first two principal components
ggplot(loadings_long, aes(x = reorder(Variable, Loading), y = Loading, fill = PC)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ PC, scales = "free_x") +
  coord_flip() +
  labs(title = "Variable Contributions to Principal Components",
       x = "Variable", y = "Loading") +
  theme_minimal()+
  theme(legend.position = "none")


####
## Variance explained
#####
# Calculate the squared loadings
squared_loadings <- md_variables_result$rotation^2

# Calculate the proportion of variance explained by each variable for each PC
proportion_variance <- apply(squared_loadings, 2, function(x) x / sum(x))
print(proportion_variance)


# plot
squared_loadings_df <- as.data.frame(md_variables_result$rotation^2)

# Create a tidy version of the data (long format)
squared_loadings_long <- squared_loadings_df  %>%
  rownames_to_column("Variable") %>%
  gather(PC, Variance, -Variable) %>%
  filter(PC %in% c("PC1", "PC2"))


# Create two separate plots and combine them
p1 <- squared_loadings_long %>%
  filter(PC == "PC1") %>%
  ggplot(aes(x = reorder(Variable, Variance), y = Variance, fill = PC)) +
  geom_bar(stat = "identity") +
  labs(title = "PC1 Contributions", x = "Variable", y = "Proportion of variance") +
  coord_flip() +
  theme_minimal()+
  theme(legend.position = "none")

p2 <- squared_loadings_long %>%
  filter(PC == "PC2") %>%
  ggplot(aes(x = reorder(Variable, Variance), y = Variance, fill = PC)) +
  geom_bar(stat = "identity", fill = "#619CFF") +
  labs(title = "PC2 Contributions", x = "Variable", y = "Proportion of variance") +
  coord_flip() +
  theme_minimal()+
  theme(legend.position = "none")

# Combine the plots using patchwork or gridExtra

combined_plot <- plot_grid(p1, p2, ncol = 2)
# Display the combined plot
combined_plot
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/pca_samples_varianceexplained.pdf", device = "pdf", width = 7, height = 5 , units = "in")


# Expand the data to all PCs
squared_loadings_long <- squared_loadings_df  %>%
  rownames_to_column("Variable") %>%
  gather(PC, Variance, -Variable) 

squared_loadings_long$PC <- as.factor(squared_loadings_long$PC)
squared_loadings_long$PC <- factor(squared_loadings_long$PC, levels = c("PC1",   "PC2" , "PC3",  "PC4",  "PC5",  
                                                                        "PC6",  "PC7" , "PC8" , "PC9",
                                                                        "PC10", "PC11", "PC12", "PC13" ,"PC14", "PC15",
                                                                        "PC16", "PC17", "PC18", "PC19","PC20",
                                                                        "PC21","PC22", "PC23" ))

ggplot(squared_loadings_long , aes(x = Variable, y = Variance, fill = PC)) +
  geom_bar(stat = "identity") +
  labs(x = "Variable", y = "Proportion of variance explained") +
  coord_flip() +
  theme_minimal()+
  facet_wrap(~PC)+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/pca_samples_varianceexplained_all.pdf", device = "pdf", width = 15, height = 17 , units = "in")


# ---------------------------------------------------------------------------------------------------------
# count number of samples
md_final_count <- md_final %>%
  count(c("general_env_feature", "vegetation_type"))

colnames(md_final)

ggplot(md_final_count, aes(x = vegetation_type, y = freq))+
  geom_col(aes(fill = vegetation_type))+
  geom_text(aes(label = freq), vjust = -0.5)+
  facet_wrap(~general_env_feature)+
  xlab("")+
  ylab("Number of samples")+
  scale_fill_manual(values = pal1)+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/samples_distribution.pdf", device = "pdf", width = 5, height = 5 , units = "in")



# ---------------------------------------------------------------------------------------------------------
# Make MAP of samples
library(maps)
library(ggrepel)
# Create a basic map of Australia
australia_map <- map_data("world")
australia_map <- australia_map %>% filter(region == "Australia")

md_final_map <- md_final %>%
  rownames_to_column(var = "Sample")
# Plot map of Australia
ggplot() +
  # Add the map
  geom_polygon(data = australia_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
  # Add points (latitude and longitude)
  geom_point(data = md_final_map , aes(x = longitude, y = latitude, color = vegetation_type, shape = general_env_feature),
             size = 5, alpha = 0.8) +
  # Add labels for the points
  #geom_text_repel(data = md_final_map , aes(x = longitude, y = latitude, label = Sample)) +
  # Customize the theme
  theme_void() +
  scale_color_manual(values = pal1)+
  scale_shape_manual(values=c(16, 13))+
  theme(legend.position = "bottom")
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/map_samples.pdf", device = "pdf", width = 10, height = 10 , units = "in")


# ---------------------------------------------------------------------------------------------------------
# Plot sequencing depth and number of samples

total_reads <- read.table("/Volumes/BunnyBike/mge_urban/base/number_reads_sample.txt", sep = " ", header = F)
colnames(total_reads) <- c("Sample", "tot_reads")

total_reads <- total_reads %>%
  mutate(Sample = gsub(Sample, pattern = "_clean_1.fastq:", replacement = "")) %>%
  mutate(Sample = gsub(Sample, pattern = "\\./", replacement = "")) 

md_final_totreads <- md_final %>%
  rownames_to_column(var = "Sample") %>%
  left_join(.,total_reads, by = "Sample")

ggplot(md_final_totreads, aes(x = vegetation_type, y = tot_reads))+
  geom_boxplot(aes(color = vegetation_type, fill = vegetation_type), alpha = 0.4, outlier.shape = NA)+
  geom_jitter(aes(color = vegetation_type))+
  facet_wrap(~general_env_feature)+
  scale_color_manual(values = pal1)+
  scale_fill_manual(values = pal1)+
  xlab("")+
  ylab("Total reads")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/total_reads_distribution.pdf", device = "pdf", width = 5, height = 5 , units = "in")


write_csv(md_final_totreads, "/Volumes/BunnyBike/mge_urban/base/md_final_final_march25.csv")


# ---------------------------------------------------------------------------------------------------------
# Plot season of collection

# Read the final final metadata (the one with total reads)
md_final_final <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_final_march25.csv") %>%
  column_to_rownames(var = "Sample") %>%
  # Add season
  mutate(
    # First convert the string dates to proper date objects
    date_clean = as_date(str_remove(collection_date, " UTC")),
    # Extract month number
    month = month(date_clean),
    # Create season based on meteorological seasons (Southern Hemisphere)
    season = case_when(
      month %in% c(12, 1, 2) ~ "Summer",
      month %in% c(3, 4, 5) ~ "Fall",
      month %in% c(6, 7, 8) ~ "Winter",
      month %in% c(9, 10, 11) ~ "Spring",
      TRUE ~ NA_character_
    )
  )

md_final_final_plot <- md_final_final %>%
  count(c("general_env_feature", "vegetation_type", "season")) %>%
  mutate(season = ifelse(is.na(season), "Unknown", season))

md_final_final_plot$season <- factor(md_final_final_plot$season, levels = c("Winter", "Spring", "Summer", "Fall", "Unknown"))
ggplot(md_final_final_plot, aes(x = season, y = freq, fill = vegetation_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~general_env_feature) +
  labs(
    x = "",
    y = "Number of samples",
    fill = "Vegetation Type"
  ) +
  scale_fill_manual(values = pal1)+
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/seasons_distribution.pdf", device = "pdf", width = 5, height = 5 , units = "in")


