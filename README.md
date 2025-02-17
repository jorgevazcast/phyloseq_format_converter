# Phyloseq Format Converter

## Requirements
- [ ] [phyloseq] (https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html)

## Add your files
Functions for converting different dataformats to a phyloseq object.

## Handling Unknown Taxa
Unknown taxa will inherit the taxonomy of the previously assigned taxonomic level.

For example:
k_Bacteria;p_Actinobacteriota;c_Coriobacteriia;o_Coriobacteriales;f_Eggerthellaceae;**g_** 
"uc_f_Eggerthellaceae"


## Basic usages (phyloseq objects)
Remove the first line of the metaphlan 4 output, its it start with an '#'

```
# Standar output
tail -n +2 mpa_v4.1.1._gtdb_vdp5yl.tsv  > mpa_v4.1.1._gtdb_vdp5yl.trim.tsv
# GTDB output
tail -n +2 mpa_v4.1.1._gtdb_vdp5yl.tsv > mpa_v4.1.1._gtdb_vdp5yl.trim.tsv

```

Convert the outfiles to a phyloseq object

```
>r
# Set seed for reproducibility
set.seed(12345)

# Load necessary libraries
library(phyloseq)

# Define the main directory path
path_dir = "./"

# Source the custom function for taxonomic assignment
source(paste0(path_dir, "/Functions/MetaPhlAn_4_preprocessing_functions.R"))

# Read the data.frame (GTDB_r207)
genus <- read_infile_MetaPhlAn_GTDB( in.file="/home/luna.kuleuven.be/u0141268/Postdoc_Raes/Projects/SEQuentialREcolonizaTion/1_infiles/FGFP5yrl/mpa_v4.1.1._gtdb_vdp5yl.trim.tsv",tax_level="genus")	

# Create the phyloseq object
genus.phylo.GTDB_r207 <- phyloseq_format_MetaPhlAn(in.data = genus , database="GTDB_r207")

genus.phylo.GTDB_r207


```
