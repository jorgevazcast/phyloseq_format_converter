set.seed(12345)
library(phyloseq)

#####################################################################
########################      MetaPhlAn 4    ########################
rar_phylo <- function(phylo.obj,Nreads){
	cat("N reads: ", Nreads, "\n")
	phylo.obj.rar = rarefy_even_depth( phylo.obj , replace=FALSE,  sample.size = Nreads )
	taxSums <- taxa_sums(phylo.obj.rar)
	phylo.obj.rar <- prune_taxa(names(taxSums[taxSums > 0]),phylo.obj.rar)

	Otable <- data.frame(t(otu_table(phylo.obj.rar)))
	(raremax <- min(rowSums(Otable)))	
	pdf( paste0("rarecurve.VDP.mOTUs.",Nreads,".pdf"),width=10)
	rarecurve(Otable, step = 20, sample = raremax, col = "blue", cex = 0.6, label = F)
	dev.off()
	
	return(phylo.obj.rar)
	
}

# MetaPhlAn reports abundances at the MOST SPECIFIC level possible.
# If reads map to a specific SGB, they're reported there, NOT at genus level.
# This means filtering for rows ending in "g__" gives incomplete data (92-100%).
#
# This function:
# 1. Takes only the deepest level (SGBs with |t__) - avoids double counting
# 2. Extracts genus from each SGB's taxonomic path
# 3. Aggregates (sums) all SGBs belonging to the same genus
#
# Result: Complete genus-level abundances that sum to ~100%

aggregate_to_genus_simple <- function(in_df) {
	in_df$clade_name <- as.character(in_df$clade_name)
	in_df <- in_df[grepl("\\|t__", in_df$clade_name), ]
	
	in_df$genus <- sapply(in_df$clade_name, function(x) {
		parts <- unlist(strsplit(x, "\\|"))
		genus_idx <- which(grepl("^g__", parts))
		paste(parts[1:genus_idx], collapse = "|")
	})

	numeric_cols <- sapply(in_df, is.numeric)
	genus_agg <- aggregate(in_df[, numeric_cols], by = list(genus = in_df$genus), FUN = sum)
	
	rownames(genus_agg) <- genus_agg$genus
	genus_agg$clade_name <- genus_agg$genus
	genus_agg$genus <- NULL
	
	return(genus_agg)
}

read_infile_MetaPhlAn <- function(in.file = "", tax_level = "genus", skip.rows = 1, unclassified = TRUE, Eukaryota = T, use_agg_genus_simple = T){
	
	# Validar inputs
	if(in.file == "") stop("Please provide input file path")
	if(!tax_level %in% c("SGB", "species", "genus")) {stop("tax_level must be one of: 'SGB', 'species', 'genus'")}
	
	# Leer archivo
	Taxa_df <- utils::read.delim(in.file, sep = "\t", check.names = FALSE, skip = skip.rows)
	Taxa_df$clade_name <- as.character(Taxa_df$clade_name)
	
	# Extraer UNCLASSIFIED
	UNCLASSIFIED <- Taxa_df[Taxa_df$clade_name == "UNCLASSIFIED", ]
	
	if(nrow(UNCLASSIFIED) > 0) {
		UNCLASSIFIED$clade_name <- NULL
		rownames(UNCLASSIFIED) <- "UNCLASSIFIED"
		message("Found UNCLASSIFIED row")
	} else {
		message("No UNCLASSIFIED row found")
		UNCLASSIFIED <- NULL
	}
	
	# Filtrar por nivel taxonómico
	if (tax_level == "SGB") {
		Taxa_df <- Taxa_df[grepl("[|]t__", Taxa_df$clade_name), ]
		Select.Taxa.level <- sapply(Taxa_df$clade_name, function(x) {
			vect <- unlist(strsplit(x, "[|]"))
			grepl("t__", vect[length(vect)])
		})
	}
	
	if (tax_level == "species") {
		Taxa_df <- Taxa_df[grepl("[|]s__", Taxa_df$clade_name), ]
		Select.Taxa.level <- sapply(Taxa_df$clade_name, function(x) {
			vect <- unlist(strsplit(x, "[|]"))
			grepl("s__", vect[length(vect)])
		})
	}
	
	if (tax_level == "genus") {
		Taxa_df <- Taxa_df[grepl("[|]g__", Taxa_df$clade_name), ]
		Select.Taxa.level <- sapply(Taxa_df$clade_name, function(x) {
			vect <- unlist(strsplit(x, "[|]"))
			grepl("g__", vect[length(vect)])
		})
	}


	### Use the function for the genus agglomeration method not ust the filtering ####		
	if(use_agg_genus_simple == T & tax_level == "genus" ){
		cat("\nuse_agg_genus_simple function\n")
		Taxa_df <- aggregate_to_genus_simple(in_df=Taxa_df)
	
	}else{
		rownames(Taxa_df) <- Taxa_df$clade_name
		Taxa_df <- Taxa_df[Select.Taxa.level == TRUE, ]

	}

	### Remove Eukaryota ####
	if(Eukaryota == F){ Taxa_df <- Taxa_df[!grepl("k__Eukaryota", Taxa_df$clade_name), ]}

	Taxa_df$clade_name <- NULL
				
	# Add the unclassifed ones
	if(unclassified == TRUE && !is.null(UNCLASSIFIED)) {  Taxa_df <- rbind(UNCLASSIFIED, Taxa_df)    }
	
	# Clean empty samples and taxa
	Taxa_df <- Taxa_df[, colSums(Taxa_df) != 0, drop = FALSE]
	Taxa_df <- Taxa_df[rowSums(Taxa_df) != 0, , drop = FALSE]
	
	# Summary
	message(sprintf("Loaded %d taxa x %d samples", nrow(Taxa_df), ncol(Taxa_df)))
	message(sprintf("Column sums range: %.2f - %.2f", min(colSums(Taxa_df)), max(colSums(Taxa_df))))
	
	return(Taxa_df)
}


read_infile_MetaPhlAn_GTDB <- function(in.file = "", tax_level = "", skip.rows = 1){  # Levels:  species, genus

	Taxa_df <- utils::read.delim(in.file, sep = "\t", check.names = F, , skip = skip.rows)
	Taxa_df$clade_name <- as.character(Taxa_df$clade_name)

		
	if (tax_level == "species") {
		Taxa_df <- Taxa_df[grepl(";s__",Taxa_df$clade_name),]
		Select.Taxa.level <- sapply(Taxa_df$clade_name, function(x){
				vect <- unlist(strsplit(x, ";"))
				TaxLevel <- grepl("^s__",vect[length(vect)])	
				return(TaxLevel)})	
	}
	if (tax_level == "genus") {
		Taxa_df <- Taxa_df[grepl(";g__",Taxa_df$clade_name),]
		Select.Taxa.level <- sapply(Taxa_df$clade_name, function(x){
				vect <- unlist(strsplit(x, ";"))
				TaxLevel <- grepl("^g__",vect[length(vect)])	
				return(TaxLevel)})			
	}
	
	row.names(Taxa_df) = Taxa_df$clade_name
	print(all(names(Select.Taxa.level) == rownames(Taxa_df)))
	Taxa_df <- Taxa_df[ Select.Taxa.level == T, ]
	Taxa_df$clade_name = NULL
	
	Taxa_df = Taxa_df[, colSums(Taxa_df) != 0]
	Taxa_df = Taxa_df[rowSums(Taxa_df) != 0, ]
	return(Taxa_df)
}

MetaPhlAn_names_databases <- function(namesTaxa = c(), DB=c("mpa","GTDB_r207")  ){

	taxa_names <- namesTaxa
	if(DB == "mpa"){
		### Convert the MetaPhlAn taxonomy to the one that of the dada2 format
		taxa_names <- gsub("k__","k_", taxa_names)
		taxa_names <- gsub("p__","p_", taxa_names)
		taxa_names <- gsub("c__","c_", taxa_names)
		taxa_names <- gsub("o__","o_", taxa_names)	
		taxa_names <- gsub("f__","f_", taxa_names)			
		taxa_names <- gsub("g__","g_", taxa_names)			
		taxa_names <- gsub("s__","s_", taxa_names)			
		taxa_names <- gsub("t__","t_", taxa_names)					
	}
	if(DB == "GTDB_r207"){
		### Convert the MetaPhlAn taxonomy to the one that of the dada2 format
		taxa_names <- gsub("d__","k_", taxa_names)
		taxa_names <- gsub("p__","p_", taxa_names)
		taxa_names <- gsub("c__","c_", taxa_names)
		taxa_names <- gsub("o__","o_", taxa_names)	
		taxa_names <- gsub("f__","f_", taxa_names)			
		taxa_names <- gsub("g__","g_", taxa_names)			
		taxa_names <- gsub("s__","s_", taxa_names)			
		taxa_names <- gsub("t__","t_", taxa_names)					
	}
	return(taxa_names)
}

phyloseq_format_MetaPhlAn <- function(in.data, database=c("mpa","GTDB_r207"), scale1M = F){

	cat("\nCheck that the taxa is in the rownames\nrownames: \n", head(rownames(in.data)),"\n\n")
	cat("\nDB taxonomy: ", database,"\n\n")

	taxa_names <- MetaPhlAn_names_databases( namesTaxa= rownames(in.data), DB=database )
	
	rownames(in.data) <- taxa_names

	if(database == "mpa"){ split_data <- stringr::str_split(rownames(in.data), "\\|") }
	if(database == "GTDB_r207"){ split_data <- stringr::str_split(rownames(in.data), ";") }
		
	max_length <- max(lengths(split_data))

	# Fill in missing elements with NA
	filled_data <- lapply(split_data, function(x) c(x, rep(NA, max_length - length(x))))
	
	# Convert the list to a data frame
	Tax_df <- as.data.frame(do.call(rbind, filled_data))	
	rownames(Tax_df) <- taxa_names	

	Tax_df[Tax_df == "t_"] <- NA
	Tax_df[Tax_df == "s_"] <- NA	
	Tax_df[Tax_df == "g_"] <- NA
	Tax_df[Tax_df == "f_"] <- NA
	Tax_df[Tax_df == "o_"] <- NA
	Tax_df[Tax_df == "c_"] <- NA
	Tax_df[Tax_df == "p_"] <- NA
	Tax_df[Tax_df == "k_"] <- NA			
		
	### Complete the NA taxonomy 
	for( i in 1:nrow(Tax_df) ){
		# i <- 927
		temp.line <- as.character(Tax_df[i,])
		if( any(is.na(temp.line)) ){
			
			na.ind  <- is.na(temp.line)			
			na.lev <- min( grep("TRUE",as.character( na.ind )) ) ### Locates the first NA taxonomic label

			if(na.lev == 1 & temp.line[na.lev] == "unassigned"){
				Tax_df[i,] <- "unassigned"
				
			}else{
				### Takes the last-known taxonomic category
				last.tax <- tail(as.character(temp.line[!na.ind]),n=1)
				 ### Add an uc_LAST_KNOWN_TAX_LEVEL to all the unkwon taxonomic levels, ie: 
				#temp.line[,(na.lev+1):length(temp.line)] <- paste0("uc_",last.tax)
				temp.line[is.na(temp.line)] <- paste0("uc_",last.tax)
				### Replace the taxonomy
				#print(temp.line)
				Tax_df[i,]<-temp.line
			}	
		}
	}

	TaxColumns<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","SGB")		
	
	colnames(Tax_df) <- TaxColumns[1:length(colnames(Tax_df))]
	
	if(scale1M == T){
		in.data <-  in.data * 10000 
	}
	
	OTUs <- phyloseq::otu_table(in.data,taxa_are_rows=T)
	TAXA <- phyloseq::tax_table(as(Tax_df,"matrix"))
	
	out.phylo <- phyloseq::phyloseq(OTUs,TAXA)
	return(out.phylo)
}

#####################################################################
########################       HUMAnN 4      ########################

read_infile_HUMAnN_4 <- function(in.file){

	Abund_df <- utils::read.delim(in.file, sep = "\t", check.names = F)

	colnames(Abund_df) <- gsub("Gene Family","names",colnames(Abund_df))
	colnames(Abund_df) <- gsub("Pathway","names",colnames(Abund_df))
	row.names(Abund_df) = Abund_df$names
	Abund_df$names = NULL

	Abund_df = Abund_df[, colSums(Abund_df) != 0] 
	Abund_df = Abund_df[rowSums(Abund_df) != 0, ]

	All <-  Abund_df
	Features <-  Abund_df[!grepl("[|]",rownames(Abund_df)),]
	Taxa_stratified <-  Abund_df[grepl("[|]",rownames(Abund_df)),]	
	
	list_ret <- list(All,Taxa_stratified,Features)
	names(list_ret) <- c("All","Taxa_stratified","Features")
	return(list_ret)	
}

phyloseq_format_HUMAnN <- function(in.data){
	cat("\nCheck that the functions are in the rownames\nrownames: \n", 
	head(rownames(in.data)), "\n\n")
	
	Tax_df <- data.frame(Function = rownames(in.data) )
	rownames(Tax_df) <- rownames(in.data)
	OTUs <- phyloseq::otu_table(in.data, taxa_are_rows = T)


	OTUs <- phyloseq::otu_table(in.data, taxa_are_rows = T)
	TAXA <- phyloseq::tax_table(as(Tax_df, "matrix"))
	out.phylo <- phyloseq::phyloseq(OTUs, TAXA)
	
	return(out.phylo)
}

create_phyloseq_files <- function(Abund.list, prefix = "Out"){

	All.pyhlo <- phyloseq_format_HUMAnN(in.data =  Abund.list$All )
	Features.pyhlo <- phyloseq_format_HUMAnN(in.data =  Abund.list$Features )
	Taxa_stratified.phylo <- phyloseq_format_HUMAnN(in.data =  Abund.list$Taxa_stratified )
	
	saveRDS(file = paste0(prefix,".All.phylo.rds"), All.pyhlo )
	saveRDS(file = paste0(prefix,".Features.phylo.rds"), Features.pyhlo )
	saveRDS(file = paste0(prefix,".Taxa_stratified.phylo.rds"), Taxa_stratified.phylo )			
}




