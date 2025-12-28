# Package loading
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(Hmisc)
library(scales)
library(ggrepel)
library(forcats)

# Quality analysis of data
# Load table with information for quality analysis
table_path <- "/mnt/lustre/groups/maier/maina479/projects/sequencing_check/scripts/table_metrics.tsv"

# Reading table and converting to R data frame
table_original <- read.table(table_path, header = TRUE, sep = "\t")

# Classifying differently table_original, new data frame for plot labeling
table_metrics <- table_original %>%
  mutate(
    SubjectLabel = case_when(
      grepl("^MSB", Subject) ~ "MSB", # If columns starts with "MSB..." name it MSB
      .default = Subject), # Otherwise name it Com20
    GroupLabel = case_when( 
      Group == "blank" ~ paste("Blank (", SubjectLabel, ")", sep = ""), # "Blank (MSB)" or "Blank (Com20)"
      Group == "DMSO control"  ~ paste("Control (", SubjectLabel, ")", sep = ""), # "Control (MSB)" or "Control (Com20)"
      .default = paste("Drug (", SubjectLabel, ")", sep = ""), # "Drug (MSB)" or "Drug (Com20)"
    )
  )

# Out of metadata create correct name for samples
metadata <- table_original %>% 
  mutate(
    SubjectLabel = case_when(
      Subject == "Com20" & Group == "DMSO control" ~ "Control Com20",
      Subject == "Com20" ~ "Com20",
      grepl("^MSB", Subject) & Group == "DMSO control" ~ paste("Control", Subject),
      grepl("^MSB", Subject) ~ Subject,
      .default = Subject), # To avoid NA for other data sets; everything else goes here
    Concentration = case_when(
      Group == "DMSO control" ~ "1% DMSO", # Control label for plot
      .default = str_extract(Group, "\\d+(\\.\\d*)?")), # Extracts concentration number from drug
    GroupLabel = case_when(
      Group == "DMSO control" ~ "Control",
      .default = str_trim(str_remove(Group, "\\d+(\\.\\d*)?"))), # Remove concentration to extract drug name
  )

# Chosen order of x-axis labels for the plots
table_metrics$GroupLabel <- factor(
  table_metrics$GroupLabel,
  levels = c("Blank (Com20)", "Blank (MSB)", "Control (Com20)", "Control (MSB)", "Drug (Com20)", "Drug (MSB)")
)

# Com20 and also MSB box plot for total sequences per sample group
total_seq_box_plot <- ggplot(table_metrics, aes(x = GroupLabel, y = Total.Sequences)) +
  geom_boxplot(fill = "#7eec7e", outlier.size = 2, outlier.alpha = 0.55) + # size and oppacity manipulation through: "https://stackoverflow.com/questions/5677885/ignore-outliers-in-ggplot2-boxplot" 
  labs(
    #title = "Total Sequences in Com20 and MSB: Blanks, Controls, and Drug Treatments",
    x = "Sample Group",
    y = "Total Sequences"
  ) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

ggsave("/mnt/lustre/groups/maier/maina479/projects/sequencing_check/output/total_seq_box_plot.png",
       plot = total_seq_box_plot, width = 8, height = 6, dpi = 300)

# Com20 and MSB total sequences median and IQR table from box plot
totseq_boxplot_table <- table_metrics %>%
  group_by(GroupLabel) %>%
  summarise(
    Median = median(Total.Sequences, na.rm = TRUE), # calc. of data: "#  https://stackoverflow.com/questions/68400943/aggregate-multiple-columns-in-a-data-frame-at-once-calculating-different-statist"
    Q1 = quantile(Total.Sequences, 0.25, na.rm = TRUE),
    Q3 = quantile(Total.Sequences, 0.75, na.rm = TRUE),
    IQR = IQR(Total.Sequences, na.rm = TRUE),
    .groups = "drop"
  )

# Save table
write_csv(totseq_boxplot_table, "/mnt/lustre/groups/maier/maina479/projects/sequencing_check/output/totseq_boxplot_table.csv")


# Com20 and MSB coverage of all sample groups
coverage_all <- ggplot(table_metrics, aes(x = GroupLabel, y = `C...`)) +
  geom_boxplot(fill = "#7eec7e", outlier.size = 2, outlier.alpha = 0.55) + # size and oppacity manipulation through: "https://stackoverflow.com/questions/5677885/ignore-outliers-in-ggplot2-boxplot" 
  labs(#title = "Coverage in Com20 and MSB: Blanks, Controls, and Drug Treatments",
       x = "Sample Group", y = "Coverage (%)") +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

ggsave("/mnt/lustre/groups/maier/maina479/projects/sequencing_check/output/coverage_all.png",
       plot = coverage_all, width = 8, height = 6, dpi = 300)

# Com20 and MSB coverage median and IQR table from box plot
cov_boxplot_table <- table_metrics %>%
  group_by(GroupLabel) %>%
  summarise(
    Median = median(`C...`, na.rm = TRUE), # calc. of data: "#  https://stackoverflow.com/questions/68400943/aggregate-multiple-columns-in-a-data-frame-at-once-calculating-different-statist"
    Q1 = quantile(`C...`, 0.25, na.rm = TRUE),
    Q3 = quantile(`C...`, 0.75, na.rm = TRUE),
    IQR = IQR(`C...`, na.rm = TRUE),
    .groups = "drop"
  )

# Save table
write_csv(cov_boxplot_table, "/mnt/lustre/groups/maier/maina479/projects/sequencing_check/output/cov_boxplot_table.csv")

# COVERAGE VS. TOTAL SEQUENCES GRAPHIC
# Filter for Com20 and no blanks
filtered_metrics_com20 <- table_metrics %>%
  filter(Group != "blank", SubjectLabel == "Com20")

# Filter for MSB and no blanks
filtered_metrics_msb <- table_metrics %>%
  filter(Group != "blank", SubjectLabel == "MSB")

# Plot Com20 Coverage vs Total Sequences
cov_tot_seq_plot_com20 <- ggplot() +
  # Com20 sample points
  geom_point(data = filtered_metrics_com20 %>% filter(SubjectLabel == "Com20"),
             aes(x = Total.Sequences, y = `C...`, color = SubjectLabel),
             alpha = 0.80,
             size = 2) +
  # Control Com20 sample points
  geom_point(data = metadata %>% filter(SubjectLabel == "Control Com20"),
             aes(x = Total.Sequences, y = `C...`, color = SubjectLabel),
             alpha = 0.75,
             size = 2) +
  # Add regression line
  geom_smooth(data = filtered_metrics_com20,
              aes(x = Total.Sequences, y = `C...`),
              method = "lm", se = FALSE, color = "black") + # one regression line according to: "https://stackoverflow.com/questions/15633714/adding-a-regression-line-on-a-ggplot"
  # Add equation of regression line
  # regression line values and R^2 indicated according to "https://rpkgs.datanovia.com/ggpubr/reference/stat_regline_equation.html"
  stat_regline_equation(data = filtered_metrics_com20,
                   aes(x = Total.Sequences, y = `C...`,
                       label = paste(..eq.label.., ..rr.label.., sep = "~~~~")),
                   formula = y ~ x,
                   label.x.npc = "center",
                   label.y.npc = "middle") +
  labs(# title = "Com20: Relationship between Coverage (%) and Total Sequences",
       x = "Total Sequences",
       y = "Coverage (%)") +
  # Legend labeling and color
  scale_color_manual(values = c("Com20" = "#7eec7e", "Control Com20" = "black"), name = "Subject") +
  # Design of plot
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

# Save plot
ggsave("/mnt/lustre/groups/maier/maina479/projects/sequencing_check/output/coverage_totseq_com20.png",
       plot = cov_tot_seq_plot_com20, width = 8, height = 6, dpi = 300)


# Only com20 filter
com20_input <- filtered_metrics_com20 %>% filter(SubjectLabel == "Com20")

# Calculating Pearson's correlation coefficient
pearson_coef <- cor(com20_input$Total.Sequences, com20_input$`C...`, method = "pearson") # according to "https://statsandr.com/blog/correlation-coefficient-and-correlation-test-in-r/"

# Print answer
print(pearson_coef)

# MSB coverage vs total sequences
# Define donor colors (shared by donors and controls)
donor_colors <- c("MSB 110" = "dodgerblue2",
                  "MSB 352" = "#E31A1C",
                  "MSB 488" = "green4",
                  "MSB 558" = "#6A3D9A",
                  "MSB 633" = "gold1",
                  "Control MSB 110" = "dodgerblue2",
                  "Control MSB 352" = "#E31A1C",
                  "Control MSB 488" = "green4",
                  "Control MSB 558" = "#6A3D9A",
                  "Control MSB 633" = "gold1")

# Shapes
subject_shapes <- c("MSB 110" = 21, "MSB 352" = 21, "MSB 488" = 21, "MSB 558" = 21, "MSB 633" = 21,
                    "Control MSB 110" = 2, "Control MSB 352" = 2, "Control MSB 488" = 2, 
                    "Control MSB 558" = 2, "Control MSB 633" = 2)

# FILLED SEPRATELY
# Plot
cov_tot_seq_plot_msb <- ggplot() +
  # Donors have filled circles
  geom_point(data = metadata %>% filter(Group != "blank", grepl("^MSB", Subject) & !grepl("Control", SubjectLabel)),
             aes(x = Total.Sequences, y = `C...`,  color = SubjectLabel, shape = SubjectLabel, fill = SubjectLabel),
             size = 3, alpha = 0.8, show.legend = FALSE) +
  
  # Controls have hollow triangles
  geom_point(data = metadata %>% filter(grepl("Control MSB", SubjectLabel)),
             aes(x = Total.Sequences, y = `C...`, shape = SubjectLabel, color = SubjectLabel),
             size = 3, alpha = 1, show.legend = TRUE) +
  # Regression line
  # one regression line according to: "https://stackoverflow.com/questions/15633714/adding-a-regression-line-on-a-ggplot"
  geom_smooth(data = filtered_metrics_msb,
              aes(x = Total.Sequences, y = `C...`),
              method = "lm", se = FALSE, color = "black") + 
  # Regression equation for plot
  # regression line values and R^2 indicated according to "https://rpkgs.datanovia.com/ggpubr/reference/stat_regline_equation.html"
  stat_regline_equation(data = filtered_metrics_msb,
                   aes(x = Total.Sequences, y = `C...`,
                       label = paste(..eq.label.., ..rr.label.., sep = "~~~~")),
                   formula = y ~ x,
                   label.x.npc = "center",
                   label.y.npc = "middle") +
  labs(#title = "MSB: Relationship between Coverage (%) and Total Sequences",
       x = "Total Sequences", y = "Coverage (%)") +
  # Legend naming and colors
  scale_fill_manual(values = donor_colors, name = "Subject") +
  scale_color_manual(values = donor_colors, name = "Subject") +
  scale_shape_manual(values = subject_shapes, name = "Subject") +
  # Override legend to have donors filled and controls hollow
  guides(
    shape = guide_legend(
      override.aes = list(
                          shape = c(rep(2, 5), rep(21, 5)),  # 5 donors with filled circles and 5 controls with triangles
                          fill  = c(donor_colors[1:5], donor_colors[1:5]),
                          color = c(donor_colors[1:5], donor_colors[1:5])
      )
    ),
    fill = "none"
  ) +
  # Plot design
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"))

# Save plot
ggsave("/mnt/lustre/groups/maier/maina479/projects/sequencing_check/output/coverage_totseq_msb.png",
       plot = cov_tot_seq_plot_msb, width = 8, height = 6, dpi = 300)
