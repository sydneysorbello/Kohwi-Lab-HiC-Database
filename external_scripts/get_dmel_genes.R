library(biomaRt)
library(dplyr)

# Connect to Ensembl for Drosophila melanogaster
ensembl <- useEnsembl(biomart = "genes", dataset = "dmelanogaster_gene_ensembl")

# Define the attributes
attributes <- c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position")

# Query the database
gene_table <- getBM(attributes = attributes, mart = ensembl)

# Filter out mitochondrial chromosome
gene_table <- gene_table %>%
  arrange(start_position) %>%
  mutate(gene_id = row_number()) %>%
  select(gene_id, external_gene_name, chromosome_name, start_position, end_position)

# Rename columns for clarity
colnames(gene_table) <- c("gene_id", "gene_name", "chromosome", "start", "end")

# Preview
head(gene_table)


write.csv(gene_table, "dmel_gene_table.csv", row.names = FALSE)

main_chromosomes <- c("X", "2L", "2R", "3L", "3R", "4", "Y")

gene_table_filtered <- gene_table %>%
  filter(chromosome %in% main_chromosomes)

write.csv(gene_table_filtered, "gene_table_filtered.csv", row.names = FALSE, quote = FALSE)



gene_table_filtered <- gene_table_filtered %>%
  mutate(gene_length = end - start + 1)  # +1 for inclusive range

# Find the largest gene
largest_gene <- gene_table_filtered %>%
  filter(gene_length == max(gene_length))

# Count genes larger than 25kb
num_large_genes <- gene_table_filtered %>%
  filter(gene_length > 25000) %>%
  nrow()

# Show results
print(largest_gene)
cat("Number of genes larger than 25kb:", num_large_genes, "\n")
