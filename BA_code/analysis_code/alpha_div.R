# Package loading
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(Hmisc)
library(scales)
library(ggrepel)
library(forcats)

# # Load table with information for quality analysis
table_path <- "/mnt/lustre/groups/maier/maina479/projects/alpha_div/scripts/table_metrics.tsv"

# Reading table and converting to R data frame
table_original <- read.table(table_path, header = TRUE, sep = "\t")

# Classifying differently table_original, new data frame for plot labeling
table_metrics <- table_original %>%
  mutate(
    SubjectLabel = case_when(
      grepl("^MSB", Subject) ~ "MSB", # If columns starts with "MSB..." name it MSB
      .default = Subject), # Otherwise name it Com20
    GroupLabel = case_when( # ERASEEE NOT USED?
      Group == "blank" ~ paste("Blank (", SubjectLabel, ")", sep = ""), # "Blank (MSB)" or "Blank (Com20)"
      Group == "DMSO control"  ~ paste("Control (", SubjectLabel, ")", sep = ""), # "Control (MSB)" or "Control (Com20)"
      .default = paste("Drug (", SubjectLabel, ")", sep = "") # "Drug (MSB)" or "Drug (Com20)"
    )
  )

# BOX PLOT DIVERSITY DRUGS
# Rename metadata for plots
table_diversity <- table_metrics %>%
  filter(Group != "blank") %>% # Exclude blanks at this point
  mutate(
    SubjectLabel = case_when( # DO I HAVE TO INCLUDE THIS? REPEATED IN TABLE_METRICS
      grepl("^MSB", Subject) ~ "MSB", 
      .default = "Com20"
    ),
    Concentration = case_when(
      Group == "DMSO control" ~ "1% DMSO", # Label for control
      .default = str_extract(Group, "\\d+(\\.\\d*)?") # Extract the concentration number from drug
    ),
    DrugLabel = case_when( 
      Group == "DMSO control" ~ "Control",
      .default = str_trim(str_remove(Group, "\\d+(\\.\\d*)?")) # Remove concentration to extract drug name
    ), 
    Treatment = case_when(
      SubjectLabel == "MSB" ~ paste(DrugLabel, Concentration), # Drug lable and concentration at the same time
      .default = paste(DrugLabel, Concentration)
    )
  )

# Com20 alpha diversity plot drugs ordered alphabetically and concentration in ascending order
com20_facet_order <- table_diversity %>%
  filter(SubjectLabel == "Com20") %>%
  mutate(
    # DrugLabel ordering by plotting the control in the first window and then alphabetically the drugs
    DrugLabel = factor(
                       DrugLabel, # To order the windows in the plot according to "https://stackoverflow.com/questions/15116081/controlling-order-of-facet-grid-facet-wrap-in-ggplot2"
                       levels = c("Control", sort(unique(DrugLabel[DrugLabel != "Control"])))), # Alphabetically ordered according to: "https://stackoverflow.com/questions/36936613/unique-sorted-rows-single-column-from-r-data-table/36951561#36951561"
    # Ascending concentrations and "1% DMSO" for the control
    Concentration = factor(
      Concentration,
      levels = c(
        sort(unique(as.numeric(Concentration[Concentration != "1% DMSO"]))),
        "1% DMSO"
      )
    )
  )

# Calculate control median diversity of Com20 for including a ashed line in alpha div boxplots
median_control_com20 <- table_diversity %>%
  filter(SubjectLabel == "Com20", DrugLabel == "Control") %>%
  pull(diversity) %>% # use a specific column: "https://www.rdocumentation.org/packages/lplyr/versions/0.1.6/topics/pull"
  median(na.rm = TRUE) # Removes missing values of data set

# Table Com20 every drug and its concentration median and IQR table
boxplot_com20_table <- com20_facet_order %>%
  group_by(DrugLabel, Concentration) %>%
  summarise(
    Median = median(diversity, na.rm = TRUE), # calc. of data: "#  https://stackoverflow.com/questions/68400943/aggregate-multiple-columns-in-a-data-frame-at-once-calculating-different-statist"
    Q1 = quantile(diversity, 0.25, na.rm = TRUE),
    Q3 = quantile(diversity, 0.75, na.rm = TRUE),
    IQR = IQR(diversity, na.rm = TRUE),
    .groups = "drop"
  )

# Save table
write_csv(boxplot_com20_table, "/mnt/lustre/groups/maier/maina479/projects/alpha_div/output/boxplot_com20_table.csv")

# Com20 alpha div boxplot
com20_alpha_div_boxplot <- ggplot(com20_facet_order, aes(x = Concentration, y = diversity)) +
  # Box plots plotting
  geom_boxplot(fill = "#7eec7e", outlier.size = 0.7, outlier.alpha = 0.55) + # size and oppacity manipulation (NO OUTLIERS HERE THOUGH!!) through: "https://stackoverflow.com/questions/5677885/ignore-outliers-in-ggplot2-boxplot"
  # Windows by drug label
  facet_wrap(~ DrugLabel, ncol = 5, scales = "free_x") +
  # y-axis set limits
  coord_cartesian(ylim = c(14, 18)) + # limits from: https://stackoverflow.com/questions/57000066/clip-only-1-axis-using-coord-cartesian
  # Control median dashed line
  geom_hline(yintercept = median_control_com20, linetype = "dashed", color = "red") + # dashed line from: "https://stackoverflow.com/questions/57177608/how-to-add-dashed-horizontal-line-with-label-in-ggplot"
  labs(
    #title = expression("Com20 Nonpareil " * N[d] * " Sequence Diversity by Drug and its Concentrations"),
    x = "Concentration (µM)",
    y = expression(bold("Metagenome " * N[d]))) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))

ggsave("/mnt/lustre/groups/maier/maina479/projects/alpha_div/output/alpha_div_com20_boxplot.png", plot = com20_alpha_div_boxplot, width = 12, height = 8, dpi = 300)

# MSB PLOT NOW
# MSB alpha diversity plot drugs ordered alphabetically and concentration in ascending order
msb_facet_order <- table_diversity %>%
  filter(SubjectLabel == "MSB") %>%
  mutate(
    # DrugLabel ordering by plotting the control in the first window and then alphabetically the drugs
    DrugLabel = factor(
                       DrugLabel, # to order the windows in the plot according to "https://stackoverflow.com/questions/15116081/controlling-order-of-facet-grid-facet-wrap-in-ggplot2"
                       levels = c("Control", sort(unique(DrugLabel[DrugLabel != "Control"])))), # alphabetically ordered according to: "https://stackoverflow.com/questions/36936613/unique-sorted-rows-single-column-from-r-data-table/36951561#36951561"
    # Ascending concentrations and "1% DMSO" for the control
    Concentration = factor(
      Concentration,
      levels = c(
        sort(unique(as.numeric(Concentration[Concentration != "1% DMSO"]))),
        "1% DMSO"
      )
    )
  )


# MSB Calculate Control median diversity for including a ashed line in alpha div boxplots
median_control_msb <- table_diversity %>%
  filter(SubjectLabel == "MSB", DrugLabel == "Control") %>%
  pull(diversity) %>%
  median(na.rm = TRUE)

# Table MSB every drug and its concentration median and IQR table
boxplot_msb_table <- msb_facet_order %>%
  group_by(DrugLabel, Concentration) %>%
  summarise(
    Median = median(diversity, na.rm = TRUE), # calc. of data: "#  https://stackoverflow.com/questions/68400943/aggregate-multiple-columns-in-a-data-frame-at-once-calculating-different-statist"
    Q1 = quantile(diversity, 0.25, na.rm = TRUE),
    Q3 = quantile(diversity, 0.75, na.rm = TRUE),
    IQR = IQR(diversity, na.rm = TRUE),
    .groups = "drop"
  )

# Save table
write_csv(boxplot_msb_table, "/mnt/lustre/groups/maier/maina479/projects/alpha_div/output/boxplot_msb_table.csv")


# MSB alpha div boxplot by drug and concentration
msb_alpha_div_boxplot <- ggplot(msb_facet_order, aes(x = Concentration, y = diversity)) +
  # Plotting of boxplots and outliers
  geom_boxplot(fill = "#7eec7e", outlier.size = 0.7, outlier.alpha = 0.55) + # size and oppacity manipulation through: "https://stackoverflow.com/questions/5677885/ignore-outliers-in-ggplot2-boxplot"
  # Windows per drug
  facet_wrap(~ DrugLabel, ncol = 5, scales = "free_x") + 
  # y-axis limits set
  coord_cartesian(ylim = c(14, 18)) + # limits from: https://stackoverflow.com/questions/57000066/clip-only-1-axis-using-coord-cartesian
  # Control median in dashed line over the whole plot
  geom_hline(yintercept = median_control_msb, linetype = "dashed", color = "red") + # dashed line from: "https://stackoverflow.com/questions/57177608/how-to-add-dashed-horizontal-line-with-label-in-ggplot"
  labs(
    # title = expression("MSB Nonpareil " * N[d] * " Sequence Diversity by Drug and its Concentrations"),
    x = "Concentration (µM)",
    y = expression(bold("Metagenome " * N[d]))) +
  # Design of plot
  theme(
    strip.text = element_text(face = "bold"),
    panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))

ggsave("/mnt/lustre/groups/maier/maina479/projects/alpha_div/output/alpha_div_msb_boxplot.png", 
       plot = msb_alpha_div_boxplot, width = 12, height = 8, dpi = 300)

# RATIO DIVERISTY PLOT FOR MSB BIG
# Calculating the average of diversity per subject
msb_avg_table <- table_diversity %>%
  filter(Subject != "Com20") %>%
  group_by(Subject, DrugLabel, Concentration) %>%
  summarise(
    DivMean = mean(diversity, na.rm = TRUE),
    .groups = "drop"
  )

# Save table
write_csv(msb_avg_table, "/mnt/lustre/groups/maier/maina479/projects/alpha_div/output/msb_avg_table.csv")

# # Calculating the diversity average ratio compared to control (drug/control)

# Calculating the average of the control
control_values <- msb_avg_table %>%
  filter(DrugLabel == "Control") %>%
  select(Subject, Concentration, Control_DivMean = DivMean)

# Joining msb_avg_table to do the ratio calculation of drug diverity by control diversty
ratio_div_table <- msb_avg_table %>%
  left_join(control_values %>% select(Subject, Control_DivMean), by = "Subject") %>%
  mutate(ratio = DivMean / Control_DivMean,
    Treatment = paste(DrugLabel, Concentration)
  )

# Save table
write_csv(ratio_div_table, "/mnt/lustre/groups/maier/maina479/projects/alpha_div/output/msb_ratio_div_table.csv")

# X-axis is now reordered by descending median value
ratio_div_table <- ratio_div_table %>%
  mutate(
    DrugConcentration = interaction(DrugLabel, Concentration, sep = " "),
    DrugConcentration = fct_reorder(DrugConcentration, ratio, .fun = median, .desc = TRUE)
  )

# MSB ratio diversity median and IQR table
boxplot_msb_ratio_table <- ratio_div_table %>%
  group_by(DrugConcentration) %>%
  summarise(
    Median = median(ratio, na.rm = TRUE), # calc. of data: "#  https://stackoverflow.com/questions/68400943/aggregate-multiple-columns-in-a-data-frame-at-once-calculating-different-statist"
    Q1 = quantile(ratio, 0.25, na.rm = TRUE),
    Q3 = quantile(ratio, 0.75, na.rm = TRUE),
    IQR = IQR(ratio, na.rm = TRUE),
    .groups = "drop"
  )

# Save table
write_csv(boxplot_msb_ratio_table, "/mnt/lustre/groups/maier/maina479/projects/alpha_div/output/ratio_table_msb_boxplot.csv")

# For labeling data in the zoomed plots, only drugs below the 0.95 line
label_data <- ratio_div_table %>%
  group_by(Treatment, DrugLabel, Concentration) %>%
  summarise(med_ratio = median(ratio, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(med_ratio)) %>% # This keeps non-missing values
  filter(med_ratio < 0.95)

## BIG MSB DIV RATIO
div_ratio_plot_msb <- ggplot(ratio_div_table, aes(x = DrugConcentration, y = ratio)) + 
  # Boxplot plotting with outliers
  geom_boxplot(fill = "#7eec7e", outlier.size = 2, outlier.alpha = 0.55) +
  # Add reference line at 0.95
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(
    # title = "MSB α-Diversity Ratio (Drug/Control) per Treatment",
    x = "Treatment",
    y = "α-Diversity Ratio",
    color = "Drug"
  ) +
  # Design of plot
  theme(
    panel.background = element_rect(fill = 'white', colour = 'black'), 
    axis.text.x = element_blank(),   # Removes all x-labels
    axis.ticks.x = element_blank(),  # Removes all x-marks
    axis.title.x = element_text(face = "bold"),  # For title of x-axis
    axis.title.y = element_text(face = "bold") # For y-axis title
  ) #+
  #coord_cartesian(ylim = c(0.85, NA)) # Comment out to do bigger supplement Graph with outlier included

# Save the plot
ggsave("/mnt/lustre/groups/maier/maina479/projects/alpha_div/output/div_ratio_msb_plot_big.png", 
      plot = div_ratio_plot_msb, width = 12, height = 8, dpi = 300)

# Filter the main table to only include DrugConcentration with med_ratio < 0.95
zoomed_drugs <- label_data$Treatment

# Order the x-axis per descending median value
ratio_div_table <- ratio_div_table %>%
  mutate(
    Treatment = fct_reorder(Treatment, ratio, .fun = median, .desc = TRUE)
  )

div_ratio_plot_zoom <- ggplot(ratio_div_table %>% filter(Treatment %in% zoomed_drugs),
                              aes(x = Treatment, y = ratio)) +
  geom_boxplot(fill = "#7eec7e", outlier.size = 2, outlier.alpha = 0.55) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(
    # title = "Diversity Ratio (Drug / Control) per Subject (Zoomed: <0.95)",
    x = "Treatment",
    y = "α-Diversity Ratio",
    color = "Drug"
  ) +
  # Plot design
  theme(
    panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.text.x = element_text(
      angle = 0,          # Horizontal
      vjust = 1,          # Adjusting vertical position
      hjust = 0.5,        # Centering horizontally
      size = 10,
      color = "black"
    ),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) +
  # Start y-axis at 0.85 for visuals
  coord_cartesian(ylim = c(0.85, NA))

# Save the zoomed plot
ggsave("/mnt/lustre/groups/maier/maina479/projects/alpha_div/output/div_ratio_plot_zoom_msb.png",
       plot = div_ratio_plot_zoom, width = 12, height = 8, dpi = 300)


# RATIO DIVERSITY ALPHA COM20
# Table for Com20 diversity
table_com20 <- table_diversity %>%
  filter(SubjectLabel == "Com20") %>%
  select(Subject, DrugLabel, Concentration, diversity) %>%
  ungroup()

# Save table
write_csv(table_com20, "/mnt/lustre/groups/maier/maina479/projects/alpha_div/output/table_div_com20.csv")

# Calculating the median of the control
control_values_com20 <- table_com20 %>%
  filter(DrugLabel == "Control") %>%
  group_by(Subject) %>%
  summarise(control_diversity = median(diversity, na.rm = TRUE), .groups = "drop")

# Joining table_com20 to do the ratio calculation of drug diverity by control diversty
ratio_div_table_com20 <- table_com20 %>%
  left_join(control_values_com20, by = "Subject") %>%
  mutate(ratio = diversity / control_diversity)

# Save table
write_csv(ratio_div_table_com20, "/mnt/lustre/groups/maier/maina479/projects/alpha_div/output/ratio_div_table_com20.csv")

# X-axis is now reordered by descending median value
ratio_div_table_com20 <- ratio_div_table_com20 %>%
  mutate(
    Treatment = paste(DrugLabel, Concentration),
    Treatment = fct_reorder(Treatment, ratio, .fun = median, .desc = TRUE)
  )

# Com20 ratio diversity mean and IQR table
com20_mean_ratio_div <- ratio_div_table_com20 %>%
  group_by(DrugLabel, Concentration) %>%
  summarise(
    MeanRatio = mean(ratio, na.rm = TRUE),
    SD = sd(ratio, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

# Save table
write_csv(com20_mean_ratio_div, "/mnt/lustre/groups/maier/maina479/projects/alpha_div/output/com20_mean_ratio_div.csv")

# For labeling data in the zoomed plots, only drugs below the 0.95 line
label_data_com20 <- ratio_div_table_com20 %>%
  group_by(Treatment, DrugLabel, Concentration) %>%
  summarise(mean_ratio = mean(ratio, na.rm = TRUE),
            sd_ratio = sd(ratio, na.rm = TRUE),
            .groups = "drop") %>%
  filter(mean_ratio < 0.95)

# Labels for big plot
all_mean_data_com20 <- ratio_div_table_com20 %>%
  group_by(Treatment, DrugLabel, Concentration) %>%
  summarise(mean_ratio = mean(ratio, na.rm = TRUE),
            sd_ratio = sd(ratio, na.rm = TRUE),
            .groups = "drop")

# Com20 all data ratio diversity
bar_plot_com20 <- ggplot(all_mean_data_com20, aes(x = Treatment, y = mean_ratio)) +
  # Bar plot of mean ratio
  geom_bar(stat = "identity", fill = "grey", color = "black") +
  # Add SD to plot
  geom_errorbar(aes(ymin = mean_ratio - sd_ratio, ymax = mean_ratio + sd_ratio), width = 0.2) +
  # Add 0.95 dashed line
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(
    #title = "Mean Diversity Ratio (Drug / Control) per Drug and Concentration",
    x = "Treatment",
    y = "Mean α-Diversity Ratio") +
  # Plot design
  theme(
    panel.background = element_rect(fill = 'white', colour = 'black'), 
    axis.text.x = element_blank(),   # Removes all x-labels
    axis.ticks.x = element_blank(),  # Removes all x-marks
    axis.title.x = element_text(face = "bold"),  # For title of x-axis
    axis.title.y = element_text(face = "bold") # For y-axis title
  ) +
  # Y-axis starts at 0.8
  coord_cartesian(ylim = c(0.8, NA))

# Save the plot
ggsave("/mnt/lustre/groups/maier/maina479/projects/alpha_div/output/div_ratio_barplot_com20.png",
       plot = bar_plot_com20, width = 12, height = 8, dpi = 300)

# Same plot but zoomed in for the drugs with ratio below 0.95
div_ratio_plot_zoom_com20 <- ggplot(label_data_com20 %>% filter(Treatment %in% zoomed_drugs),
  aes(x = Treatment, y = mean_ratio)
) +
  # Add mean bars
  geom_bar(stat = "identity", fill = "#7eec7e", color = "black") +
  # Add SD
  geom_errorbar(aes(ymin = mean_ratio - sd_ratio, ymax = mean_ratio + sd_ratio), width = 0.2) +
  # Add dashed 0.95 line
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(
    # title = "Diversity Ratio (Drug / Control) per Subject (Zoomed: <0.95)",
    x = "Treatment",
    y = "Mean α-Diversity Ratio",
    color = "Drug"
  ) +
  # Design of plot
  theme(
    panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.text.x = element_text(
      angle = 0,          # Horizontal
      vjust = 1,          # Adjusting vertical position
      hjust = 0.5,        # Centering horizontally
      size = 10,
      color = "black"
    ),
    axis.ticks.x = element_blank(),  # Removes all x-marks
    axis.title.x = element_text(face = "bold"),  # For title of x-axis
    axis.title.y = element_text(face = "bold") # For y-axis title
  ) +
  # Y-axis starts at 0.8
  coord_cartesian(ylim = c(0.8, NA))



# Save the zoomed plot
ggsave("/mnt/lustre/groups/maier/maina479/projects/alpha_div/output/div_ratio_barplot_zoom_com20.png",
       plot = div_ratio_plot_zoom_com20, width = 12, height = 8, dpi = 300)