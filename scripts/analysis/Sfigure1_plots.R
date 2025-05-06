


pal1 <- c("#D81B60", "#997404", "#004D40")
pal2 <- c("lightpink", "#D81B60","#FFC107", "#997404", "green4","#004D40")

md <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_final_votus.csv")

ggplot(md, aes(x = vegetation_type, fill = vegetation_type)) + 
  geom_bar(stat = "count") +
  stat_count(geom = "text", colour = "white", size = 3.5,
             aes(label = ..count..),position=position_stack(vjust=0.5)) +
  facet_wrap(~general_env_feature)+
  scale_fill_manual(values = c("#D81B60", "#997404","#004D40"))+
  xlab("")+
  ylab("Number of samples")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/samples_number.pdf", device = "pdf", width = 5, height = 5 , units = "in")


# testing the relationships between variables
md_variables <- md %>%
  dplyr::select(sample_id, where(is.numeric)) %>%
  dplyr::filter(!sample_id == "12891_1") %>% # this sample had a bunch of NA values, so it makes plotting very difficult
  dplyr::select(-sample_id) %>%
  dplyr::select(-c("depth", "depth_lower", "depth_upper", "slope_aspect", "slope_gradient", 
            "water_content", "sand", "silt", "fine_sand", "clay", "coarse_sand", 
            "ph_solid_h2o", "latitude", "longitude", "plasmid_rich", "plasmid_shannon", 
            "virus_rich", "virus_shannon", "Axis02", "plasmids_relab", "virus_av_len", "Axis01",
            "month", "av_len", "ratio_vt", "Sample", "tot_reads", "potassium_colwell")) # these variables had a lot of NAs too



######
# Corplot
######
corrplot::corrplot(cor(md_variables, use = "complete.obs"), 
                   type = "lower", 
                   tl.col = "black", 
                   order = "hclust", 
                   diag = FALSE)


#####
# PCA
######
# Standardize the data (PCA works best with standardized data)
md_variables_scaled <- scale(md_variables)

# run PCA
md_variables_result <- prcomp(md_variables_scaled , center = TRUE, scale. = TRUE)
# see summary results
summary(md_variables_result)

# make sure they are in the same order
all(rownames(md) == rownames(as.data.frame(md_variables_result$x)))
all(rownames(md) == rownames(as.data.frame(pca_scores)))

# Extract scores and add groupings
pca_scores <- as.data.frame(md_variables_result$x)
pca_scores$general_env_feature <-md$general_env_feature
pca_scores$vegetation_type <-md$vegetation_type

# Extract PCA loadings (arrows for variables)
pca_loadings <- as.data.frame(md_variables_result$rotation)
pca_loadings$Variable <- rownames(pca_loadings)  # Variable names

# Scaling the loadings for better visualization
arrow_scale <- 20 # Adjust scaling factor as needed
pca_loadings <- pca_loadings %>%
  mutate(PC1 = PC1 * arrow_scale, PC2 = PC2 * arrow_scale)


pca_scores$combined_var <- paste(pca_scores$general_env_feature, pca_scores$vegetation_type)
pca_scores$combined_var <- factor(pca_scores$combined_var, levels = c("Arid Grassland", "Temperate Grassland",
                                                                              "Arid Shrubland", "Temperate Shrubland", 
                                                                              "Arid Woodland", "Temperate Woodland"))

# Create PCA Biplot with ggplot2
ggplot() +
  # Plot PCA scores (points for samples)
  geom_point(data = pca_scores, aes(x = PC1, y = PC2, 
                                    color = combined_var, shape = general_env_feature),
             size = 5) +
  
  # Add PCA loadings (arrows for variables)
  geom_segment(data = pca_loadings, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black", linewidth = .8) +
  
  # Add labels to arrows
  geom_text(data = pca_loadings, aes(x = PC1, y = PC2, label = Variable), 
            vjust = 1, hjust = 1, size = 5, color = "black") +
  
  # Titles and themes
  labs(title =,x = "PC1", y = "PC2") +
  theme_bw() +
  scale_color_manual(values = pal2)+
  scale_shape_manual(values=c(18, 16))+
  theme(legend.position = "none", text = element_text(size = 16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/pca_samples.pdf", device = "pdf", width = 8, height = 8 , units = "in")


# # Extract variance explained
variance_explained <- md_variables_result$sdev^2
prop_variance <- variance_explained / sum(variance_explained)
cumulative_variance <- cumsum(prop_variance)

# Create a data frame with this information
variance_df <- data.frame(
  PC = 1:length(variance_explained),
  Variance = variance_explained,
  PropVariance = prop_variance,
  CumulativeVariance = cumulative_variance
)

# View the results
print(variance_df)



# Make MAP of samples
library(maps)
library(ggrepel)
# Create a basic map of Australia
australia_map <- map_data("world")
australia_map <- australia_map %>% filter(region == "Australia")

md_final_map <- md 
# Plot map of Australia
ggplot() +
  # Add the map
  geom_polygon(data = australia_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
  # Add points (latitude and longitude)
  geom_point(data = md_final_map , aes(x = longitude, y = latitude, color = combined_var, shape = general_env_feature),
             size = 8, alpha = 0.8) +
  # Add labels for the points
  #geom_text_repel(data = md_final_map , aes(x = longitude, y = latitude, label = Sample)) +
  # Customize the theme
  theme_void() +
  scale_color_manual(values = pal2)+
  scale_shape_manual(values=c(18, 16))+
  theme(legend.position = "none")
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/map_samples.pdf", device = "pdf", width = 10, height = 10 , units = "in")
