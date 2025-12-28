# Package loading
library(tidyverse)
library(glue)
library(vegan)

# Read table needed for analysis
taxpasta_original <- read_tsv("test_metemgee/data/taxprofiler/taxpasta/bracken_B_standard_16gb.tsv")

# Analysis was done with the following source: https://www.youtube.com/watch?v=oLf0EpMJ4yA&ab_channel=RiffomonasProject
# TRANSPOSED
# Wide matrix, transposed
taxpasta_matrix <- taxpasta_original %>%
  select(name, starts_with("JLa")) %>% 
  pivot_longer(
    cols = -name,
    names_to = "Sample",
    values_to = "Reads") %>%
  group_by(Sample, name) %>% # Correct error "Can't convert `fill`<double> to <list>" with post: "https://stackoverflow.com/questions/63365516/tidyr-pivot-wider-error-cant-convert-double-to-list"
  summarise(Reads = sum(Reads), .groups = "drop") %>% # group-by and summarise to combine duplicates (duplicate row issue needs to be aggregated) before they are pivoted
  pivot_wider(
    names_from = name,
    values_from = Reads,
    values_fill = 0
  ) # New transposed matrix alphabetically and numerically (sample) ordered

# Path to table_metrics containing metadata
table_metrics <- "/mnt/lustre/groups/maier/maina479/projects/detaxizer/code/table_metrics.tsv"

# Change metadata to wanted labeling
metadata <- read.table(
  table_metrics,
  header = TRUE,
  sep = "\t") %>%
  rename(Sample = ID) %>% # Column with sample names has to be named correctly
  mutate(
    SubjectLabel = case_when(
      Subject == "Com20" & Group == "DMSO control" ~ "Control Com20",
      Subject == "Com20" ~ "Com20",
      grepl("^MSB", Subject) & Group == "DMSO control" ~ paste("Control", Subject),
      grepl("^MSB", Subject) ~ Subject,
      .default = Subject), # To avoid NA for other data sets; everything else goes here
    Concentration = case_when(
      Group == "DMSO control" ~ "1% DMSO", # Label for control
      .default = str_extract(Group, "\\d+(\\.\\d*)?")), # Extract concentration number from drug
    DrugLabel = case_when(
      Group == "DMSO control" ~ "Control",
      .default = str_trim(str_remove(Group, "\\d+(\\.\\d*)?")), # Remove concentration to extract drug name
    )
  )

# Matrix for Com20
edit_matrix_com20 <- taxpasta_matrix %>%
  mutate(Sample = sub("_R1.*", "", Sample)) %>% # ID is extracted from sample name
  left_join(metadata, by = "Sample") %>% # Metadata info is added
  filter(grepl("Com20", SubjectLabel)) %>%  # Com20 filter
  filter(Group != "blank") # No blanks

# Matrix for MSB
edit_matrix_msb <- taxpasta_matrix %>%
  mutate(Sample = sub("_R1.*", "", Sample)) %>% # ID is extracted from sample name
  left_join(metadata, by = "Sample") %>% # Metadata info is added
  filter(grepl("MSB", SubjectLabel)) %>% # MSB filter
  filter(Group != "blank") # No blanks

# Com20 manipulate data to filter out samples greater equal 1M seq
sample_counts_com20 <- edit_matrix_com20 %>%
  select(Sample, where(is.numeric)) %>%
  mutate(across(where(is.numeric), as.integer)) %>%
  mutate(N = rowSums(across(-Sample))) %>%  # All columns except sample are summed
  filter(N >= 1000000) %>% # Samples >= 1M are filtered
  select(-N)

# MSB manipulate data to filter out samples greater equal 1M seq
sample_counts_msb <- edit_matrix_msb %>%
  select(Sample, where(is.numeric)) %>%
  mutate(across(where(is.numeric), as.integer)) %>%
  mutate(N = rowSums(across(-Sample))) %>%  # All columns except sample are summed
  filter(N >= 1000000) %>% # Samples >= 1M are filtered
  select(-N)

# Distance matrix Com20
# Convert to a matrix format
# avgdist() function needs taxa as columns and samples as row names
tibble_matrix_com20 <- sample_counts_com20 %>% # for getting distances
  column_to_rownames("Sample")

# To have max comparability and reproducibility, the same seed is used for Com20 and MSB
set.seed(1)
# Use bray–curtis dissimilarity method with rarefaction
dist_com20 <- avgdist( # To add rarefaction to distance calculation
  x = tibble_matrix_com20,
  dmethod = "bray",
  sample = 100000,  # Sample = 1170895 is min after 1M, only 100,000 needed for good results
  iterations = 1) # iterations = 1 Number of iterations using the avgdist() function is by default 100 and can be set to a 1000 for example, but that take much more time and does not make a big difference in the result.

# Convert the bray–curtis dissimilarity object to a matrix format for further analysis
dist_matrix_com20 <- dist_com20 %>%
  as.matrix()

# # Save redable, tabs and file included
# write.table(dist_matrix_com20,
#             file = "/mnt/lustre/groups/maier/maina479/projects/beta_div/output/dist_matrix_com20.tsv", sep = "\t", quote = FALSE, col.names = NA)

# Distance Matrix MSB
tibble_matrix_msb <- sample_counts_msb %>% # Data frame instead of tibl for input for avgdist desk??ersae? for getting distances
  select(Sample, where(is.numeric)) %>%
  column_to_rownames("Sample")

# To have max comparability and reproducibility, the same seed is used for Com20 and MSB
set.seed(1)
# Use bray–curtis dissimilarity method with rarefaction
dist_msb <- avgdist( # To add rarefaction to distance calculation
                           x = tibble_matrix_msb,
                           dmethod = "bray",
                           sample = 100000,  # Sample = 1170895 is min after 1M, only 100,000 needed for good results
                           iterations = 1) # iterations = 1 Number of iterations using the avgdist() function is by default 100 and can be set to a 1000 for example, but that take much more time and does not make a big difference in the result.

# Convert the bray–curtis dissimilarity object to a matrix format for further analysis
dist_matrix_msb <- dist_msb %>%
  as.matrix()

# # Save redable, tabs and file included
# write.table(dist_matrix_msb,
#             file = "/mnt/lustre/groups/maier/maina479/projects/beta_div/output/dist_matrix_msb.tsv", sep = "\t", quote = FALSE, col.names = NA)

# PCoA Com20 classical multidimensional scaling
bray_pcoa_obj_com20 <- cmdscale(dist_matrix_com20, k = 4, eig = TRUE, add = TRUE)

# PCoA Com20 classical multidimensional scaling
bray_pcoa_obj_msb <- cmdscale(dist_matrix_msb, k = 4, eig = TRUE, add = TRUE)

# Single data frame with sample and pcoa data is created
bray_pcoa_com20 <- bray_pcoa_obj_com20$points %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  left_join(metadata, by = "Sample")  # Merges pcoa coordinates with the sample metadata by comparing with sample ID

# # Save Com20 pcoa data frame
# write_tsv(bray_pcoa_com20, "/mnt/lustre/groups/maier/maina479/projects/beta_div/output/bray_pcoa_com20.tsv")

# Single data frame with sample and pcoa data is created
bray_pcoa_msb <- bray_pcoa_obj_msb$points %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  left_join(metadata, by = "Sample")  # Merges pcoa coordinates with the sample metadata by comparing with sample ID

# # Save MSB pcoa data frame
# write_tsv(bray_pcoa_msb, "/mnt/lustre/groups/maier/maina479/projects/beta_div/output/bray_pcoa_msb.tsv")

# Extract proportion of variance of each principal component
pc1_com20 <- round(bray_pcoa_obj_com20$eig[1] / sum(bray_pcoa_obj_com20$eig), 4) * 100
pc2_com20 <- round(bray_pcoa_obj_com20$eig[2] / sum(bray_pcoa_obj_com20$eig), 4) * 100

# Labels for axes of the plot
labels <- c(glue("PCo1 ({pc1_com20}%)"), # Label: https://bioinfo.cd-genomics.com/principal-co-ordinates-analysis.html#:~:text=Among%20them%2C%20PCoA%20Axis%201,represents%20the%20distance%20between%20samples.
            glue("PCo2 ({pc2_com20}%)"))

c25 <- c( # Color palette with more differences according to: "https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes"
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1", # 6 and 7
  "skyblue2", "#FB9A99", # pink
  "#7eec7e", #palegreen2
  "#CAB2D6", # purple #11
  "#FDBF6F", # orange
  "gray70", "khaki2", #13 14
  "maroon", "orchid1", "deeppink1", "blue1", "#457aa7", #19 steelblue4
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown" #24
)

# Select 14 distinct colors to differentiate treatments
drug_palette <- c25[c(1,2,3,4,7,5,8,9,15,10,16,17,18,24,20,25)]

# Shape for Com20 data
subject_shape_com20 <- c("Com20" = 21, "Control Com20" = 1)

# Pcoa plot for Com20 and control
pcoa_com20 <- ggplot() +
  # Com20 data points
  geom_point(data = bray_pcoa_com20 %>% filter(SubjectLabel == "Com20"),
             aes(x = V1, y = V2, color = Concentration, fill = Concentration, shape = SubjectLabel),
             alpha = 0.55,
             size = 3,
             show.legend = c(fill = TRUE, shape = TRUE, color = FALSE)) +  # Crucial for displaying correct fill and color in legend
  # Com20 Control data points
  geom_point(data = bray_pcoa_com20 %>% filter(SubjectLabel == "Control Com20") %>%
               select(-DrugLabel), # add controls in each Drug plot according to: "https://stackoverflow.com/questions/62656776/best-way-to-add-points-to-all-facets-in-ggplot2-in-r"
             aes(x = V1, y = V2, shape = SubjectLabel),
             color = "black",
             fill = "black",
             alpha = 0.6,
             size = 3,
             show.legend = c(shape = TRUE, color = TRUE, fill = TRUE) ) + # Crucial for displaying correct fill and color in legend
  # Windows per drug
  facet_wrap(~ DrugLabel, ncol = 5) +
  # X- and Y-axis labels
  labs(x = labels[1], y = labels[2]) + #,
      # title = "PCoA - PCo1 vs PCo2: Com20 β-Diversity across 14 Drugs, Concentrations, and Controls") +
  # Plot outline colors (no legend)
  scale_color_manual(values = drug_palette, na.value = "grey50", breaks = c("1.25","5","20","80","320")) + 
  # Plot fill colors for drug and the legend 
  scale_fill_manual(values = drug_palette, na.value = "grey50", breaks = c("1.25","5","20","80","320"), name = "Concentration") +   
  # Plot shape for subject and the legend
  scale_shape_manual(values = c("Com20" = 21, "Control Com20" = 0), name = "Subject") + 
  # Design of plot
  theme(panel.background = element_rect(fill = 'white', color = 'black')) +
  theme(text = element_text(size = 14),
        legend.position = "right")

# Save the plot
ggsave("/mnt/lustre/groups/maier/maina479/projects/beta_div/output/pcoa_com20.png",
       plot = pcoa_com20, width = 12, height = 8, dpi = 300)


# BELOW PCOA MSB
# # Subjects must have a shaped defined
subject_shape_msb <- c(
  "Control MSB 110" = 1,
  "Control MSB 352" = 0,
  "Control MSB 488" = 5,
  "Control MSB 558" = 2,
  "Control MSB 633" = 6,
  "MSB 110" = 21,
  "MSB 352" = 22,
  "MSB 488" = 23,
  "MSB 558" = 24,
  "MSB 633" = 25
)

# Pcoa plot for each drug
pcoa_msb <- ggplot() +
  # MSB subjects plotting
  geom_point(data = bray_pcoa_msb %>% filter(SubjectLabel == Subject),
             aes(x = V1, y = V2, color = Concentration, fill = Concentration, shape = SubjectLabel),
             alpha = 0.55,
             size = 3,
             show.legend = c(fill = TRUE, shape = TRUE, color = FALSE)) +  # Crucial for displaying only fill shape in legend
  # MSB Control plotting
  geom_point(data = bray_pcoa_msb %>% filter(grepl("Control", SubjectLabel)) %>% # Each control from each donor should have its own shape
               select(-DrugLabel), # add controls in each Drug plot according to: "https://stackoverflow.com/questions/62656776/best-way-to-add-points-to-all-facets-in-ggplot2-in-r"
             aes(x = V1, y = V2, shape = SubjectLabel),
             color = "black",
             fill = NA,
             alpha = 0.6,
             size = 3,
             show.legend = c(shape = TRUE, color = FALSE, fill = FALSE)) + # Crucial for displaying only fill shape in legend
  # Windows per drug
  facet_wrap(~ DrugLabel, ncol = 5) +
  # X- and Y-axis labels
  labs(x = labels[1], y = labels[2]) + #,
       # title = "PCoA - PCo1 vs PCo2: MSB β-diversity across 14 Drugs, Concentrations, and Donor Controls") +
  # Plot outline colors (no legend)
  scale_color_manual(values = drug_palette, na.value = "grey50", breaks = c("1.25","5","20","80","320")) +  
  # Plot fill colors for drug and the legend
  scale_fill_manual(values = drug_palette, na.value = "grey50", breaks = c("1.25","5","20","80","320"), name = "Concentration") +
  # Plot shape for subject and show in the legend
  scale_shape_manual(values = subject_shape_msb, name = "Subject") + 
  # Design of plot
  theme(panel.background = element_rect(fill = 'white', color ='black')) +
  theme(text = element_text(size = 14),
        legend.position = "right") +
  guides( # Correct fill and shape in legend according to: "https://stackoverflow.com/questions/43266636/incorrect-shape-and-fill-of-ggplot-legend"
    fill = guide_legend( # Legend for Drug is now with colors and not in black thanks to "https://stackoverflow.com/questions/74368241/ggplot-points-have-color-but-legend-guide-points-are-all-black"
      override.aes = list(shape = 21, color = NA)),  # Concentration in legend is shown as circles with filled colors
    shape = guide_legend(
      override.aes = list(
        fill = c(NA, NA, NA, NA, NA, "black", "black", "black", "black", "black"), # MSB samples are filled and controls not
        color = "black"
      )
    )
  )


# Save the plot
ggsave("/mnt/lustre/groups/maier/maina479/projects/beta_div/output/pcoa_msb.png",
       plot = pcoa_msb, width = 12, height = 8, dpi = 300)

# Shapes for subjects and control
subject_shape_donor <- c(
  "Control MSB 110" = 0,
  "Control MSB 352" = 0,
  "Control MSB 488" = 0,
  "Control MSB 558" = 0,
  "Control MSB 633" = 0,
  "MSB 110" = 21,
  "MSB 352" = 21,
  "MSB 488" = 21,
  "MSB 558" = 21,
  "MSB 633" = 21
)

# MSB PER DONOR
# Both controls and drugs must share the same ID (DonorID)
bray_pcoa_msb <- bray_pcoa_msb %>%
  mutate(DonorID = gsub("Control ", "", SubjectLabel))   # Only MSB with ID

# Shapes for control and drug samples
subject_shape_msb <- c("Control" = 0, "Drug" = 21)

# MSB pcoa plot per donor
pcoa_msb_donor <- ggplot(data = bray_pcoa_msb, aes(x = V1, y = V2)) +
  # Samples are circles colored by drug
  geom_point(data = bray_pcoa_msb %>% filter(!grepl("Control", SubjectLabel)),
             aes(color = DrugLabel, fill = DrugLabel, shape = "Drug"),
             alpha = 0.75, size = 3,
             show.legend = c(fill = TRUE, shape = FALSE, color = FALSE)) + # Drugs are filled in legend
  # Controls are hollow squares
  geom_point(data = bray_pcoa_msb %>% filter(grepl("Control", SubjectLabel)),
             aes(shape = "Control"),  
             color = "black", fill = NA,
             alpha = 0.9, size = 3,
             show.legend = c(shape = TRUE, color = TRUE, fill = TRUE)) + # Label control hollow
  # Windows per donor
  facet_wrap(~ DonorID, ncol = 5) +   # 5 donor windows
  # X- and Y-axis labels
  labs(x = labels[1], y = labels[2]) +
  # Plot outline colors
  scale_color_manual(values = drug_palette, na.value = "grey50") +
  # Plot fill colors for drug and the legend
  scale_fill_manual(values = drug_palette, na.value = "grey50", name = "Drug") +
  # Plot shape for subject and show in the legend
  scale_shape_manual(
    values = c("Control" = 0, "Drug" = 21),
    name = "Experiment Group") +
  # Design plot
  theme(
    panel.background = element_rect(fill = 'white', color = 'black'),
    text = element_text(size = 14),
    legend.position = "right"
  ) +
  # Correct legend with fill, color and shape
  guides(
    fill = guide_legend(override.aes = list(shape = 21, color = NA)),
    shape = guide_legend(
      override.aes = list(fill = "black", color = "black")
    )
  )

# Save the plot
ggsave("/mnt/lustre/groups/maier/maina479/projects/beta_div/output/pcoa_per_Donor.png",
       plot = pcoa_msb_donor, width = 12, height = 8, dpi = 300)

