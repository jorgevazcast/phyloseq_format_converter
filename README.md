# Phyloseq Format Converter

Utility functions for converting multiple data formats into phyloseq objects.

## Requirements
- [phyloseq](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html)

## Handling Unknown Taxa
Unknown taxa inherit the taxonomy of the closest assigned higher taxonomic level.

For example:

k_Bacteria;p_Actinobacteriota;c_Coriobacteriia;o_Coriobacteriales;f_Eggerthellaceae;g_  
**uc_f_Eggerthellaceae**

## Usage

Load required libraries and functions:
```r
# Set seed for reproducibility
set.seed(12345)

# Load necessary libraries
library(phyloseq)

# Define the main directory path
path_dir <- "./"

# Source helper functions for MetaPhlAn and HUMAnN 4 preprocessing
source(paste0(path_dir, "/Functions/MetaPhlAn_HUMAnN_4_preprocessing_functions.R"))
```

---

## Convert MetaPhlAn output to phyloseq object
```r
# SGB-level
SGB <- read_infile_MetaPhlAn(
  in.file = "./Example_data/Example_mpa.tsv",
  tax_level = "SGB"
)

SGB.phylo <- phyloseq_format_MetaPhlAn(
  in.data = SGB,
  database = "mpa"
)

# Species-level
species <- read_infile_MetaPhlAn(
  in.file = "./Example_data/Example_mpa.tsv",
  tax_level = "species"
)


species.phylo <- phyloseq_format_MetaPhlAn(
  in.data = species,
  database = "mpa"
)


# Genus-level
genus <- read_infile_MetaPhlAn(
  in.file = "./Example_data/Example_mpa.tsv",
  tax_level = "genus"
)

genus.phylo <- phyloseq_format_MetaPhlAn(
  in.data = genus,
  database = "mpa"
)

# Assign genus names as taxa identifiers
taxa_names(genus.phylo) <- tax_table(genus.phylo)[, 6]
```

---

## Convert MetaPhlAn GTDB_r207 output to phyloseq object
```r
# Species-level
species <- read_infile_MetaPhlAn_GTDB(
  in.file = "./Example_data/Example_mpa.GTDB_r207.tsv",
  tax_level = "species"
)

species.phylo.GTDB_r207 <- phyloseq_format_MetaPhlAn(
  in.data = species,
  database = "GTDB_r207"
)

# Genus-level
genus <- read_infile_MetaPhlAn_GTDB(
  in.file = "./Example_data/Example_mpa.GTDB_r207.tsv",
  tax_level = "genus"
)

genus.phylo.GTDB_r207 <- phyloseq_format_MetaPhlAn(
  in.data = genus,
  database = "GTDB_r207"
)

# Assign genus names as taxa identifiers
taxa_names(genus.phylo.GTDB_r207) <- tax_table(genus.phylo.GTDB_r207)[, 6]
```

---

## Convert HUMAnN 4 output to phyloseq object

The `read_infile_HUMAnN_4()` function reads HUMAnN 4 output files and returns a list with three data frames:

| Element | Description |
|---------|-------------|
| `All` | Complete data including all features and taxa stratification |
| `Features` | Only unstratified features (e.g., KOs, pathways) |
| `Taxa_stratified` | Only taxa-stratified features (contains `\|` in row names) |

### Example: KO abundance per taxa
```r
# Read HUMAnN 4 output file
KO_taxa_list <- read_infile_HUMAnN_4("./Example_data/HUMAnN.KO_Taxa.tsv")

# Convert each element to phyloseq object individually
All.phylo <- phyloseq_format_HUMAnN(in.data = KO_taxa_list$All)
Features.phylo <- phyloseq_format_HUMAnN(in.data = KO_taxa_list$Features)
Taxa_stratified.phylo <- phyloseq_format_HUMAnN(in.data = KO_taxa_list$Taxa_stratified)
```

### Using the wrapper function

The `create_phyloseq_files()` function is a convenient wrapper that processes all elements at once and saves them as `.rds` files:
```r
# Create and save all phyloseq objects with a custom prefix
create_phyloseq_files(Abund.list = KO_taxa_list, prefix = "Out")

# This generates:
# - Out_All.rds
# - Out_Features.rds
# - Out_Taxa_stratified.rds
```

---

## Functions Reference

| Function | Description |
|----------|-------------|
| `read_infile_MetaPhlAn_GTDB()` | Reads MetaPhlAn output with GTDB taxonomy |
| `phyloseq_format_MetaPhlAn()` | Converts MetaPhlAn data to phyloseq object |
| `read_infile_HUMAnN_4()` | Reads HUMAnN 4 output (gene families, pathways, KOs) |
| `phyloseq_format_HUMAnN()` | Converts HUMAnN data to phyloseq object |
| `create_phyloseq_files()` | Wrapper to batch convert and save phyloseq objects |
