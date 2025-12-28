# Make a table with important information for the analysis
# Read nonpareil_all_samples.tsv + transform in table
# Path to file
nonpareil_path <- "/mnt/lustre/groups/maier/maina479/projects/test_metemgee/data/taxprofiler/nonpareil/nonpareil_all_samples.tsv"

# Read file and draw sample names (ID)
# The sample names are stored as row names in the .tsv. read.delim assumes HEADER and tabulator
nonpareil <- read.delim(nonpareil_path,
                  row.names = 1) # 1. column treated as a row name

# Row names into new column for correct reading
nonpareil$sample_name <- rownames(nonpareil)

# Extract ID
nonpareil$ID <- sub("_R1.*", "", nonpareil$sample_name)

# C (coverage) in percentage
nonpareil$`C(%)` <- nonpareil$C * 100

# Total bases sequenced in Gbp
nonpareil$`LR(Gbp)` <- nonpareil$LR / 1e9

# Read the full TSV file (multiqc_fastqc_1.txt) and convert into table (for extracting "Total sequences")
multiqc_path <- "/mnt/lustre/groups/maier/maina479/projects/test_metemgee/data/taxprofiler/multiqc/multiqc_data/multiqc_fastqc_1.txt"

multiqc <- read.table(
  file = multiqc_path,
  header = TRUE,
  sep = "\t"
)

# Extract ID
multiqc$ID <- sub("_R1", "", multiqc$Sample)

metadata_all = "/mnt/lustre/groups/maier/maina479/projects/sequencing_check/scripts/metadata_all.csv"

metadata <- read.table(
  file = metadata_all,
  header = TRUE,
  sep = ";",
  colClasses = c("character", "character", "character"))

metadata$ID <- (metadata$samples)

# Combine the three data frames
combine <- merge(merge(nonpareil, multiqc, by = "ID"), metadata, by = "ID")

# Change column order as wished:
table <- combine[, c("ID", "Subject", "Group", "C", "C(%)", "LR(Gbp)", "Total.Sequences", "LRstar", "diversity")]

output_path <- "/mnt/lustre/groups/maier/maina479/projects/sequencing_check/scripts/table_metrics.tsv"

# Save table
write.table(table, output_path, sep = "\t", row.names = FALSE, quote = FALSE)

