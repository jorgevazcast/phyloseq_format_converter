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

# Source helper functions for MetaPhlAn GTDB preprocessing
source(paste0(path_dir, "/Functions/MetaPhlAn_HUMAnN_4_preprocessing_functions.R"))
```

Convert the MPA GTDB_r207 outfiles to a phyloseq object

```
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

# Assign genus names
taxa_names(genus.phylo.GTDB_r207) <- tax_table(genus.phylo.GTDB_r207)[, 6]
```
			
