# Phyloseq Format Converter

## Requirements
- [ ] [phyloseq] (https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html)

## Add your files
Functions for converting different dataformats to a phyloseq object.

## Handling Unknown Taxa
Unknown taxa will inherit the taxonomy of the previously assigned taxonomic level.

For example:
k_Bacteria;p_Actinobacteriota;c_Coriobacteriia;o_Coriobacteriales;f_Eggerthellaceae;**g_** 
**"uc_f_Eggerthellaceae"**


## Basic usages (phyloseq objects)
Remove the first line of the metaphlan 4 output, its it start with an '#'

```
# GTDB output
tail -n +2 ./Example_data/mpa_v4.1.1._gtdb_vdp5yl.tsv > mpa_v4.1.1._gtdb_vdp5yl.trim.tsv

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
genus <- read_infile_MetaPhlAn_GTDB( in.file="mpa_v4.1.1._gtdb_vdp5yl.trim.tsv",tax_level="genus")	

# Create the phyloseq object
genus.phylo.GTDB_r207 <- phyloseq_format_MetaPhlAn(in.data = genus , database="GTDB_r207")

genus.phylo.GTDB_r207
			
```

## Verify the compatibility with the enterotype package

To verify the compatibility with the enterotype package, I am using the taxonomy of the GTDB_VDP5YL cohort as an example to predict its enterotypes. **This is not optimal, as the predictor was trained with the FGFP 3000 dataset using the GTDB_r86 database**. This is just to check if the taxonomy format is compatible with the enterotype pipeline

```
library(Enterotypes)
Path_db = "/home/luna.kuleuven.be/u0141268/github_projects/Enterotype_data/GTDB_r86/"

ent_list <- Enterotype(input.obj=genus.phylo.GTDB_r207, prefix_cluster_files="TEST", 
			Enterotype_detection_method = "predict_ent",background_population_dir = Path_db)
			
```			
			
			
			
