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

read_infile_MetaPhlAn <- function(in.file = "", tax_level = "", skip.rows = 1){  # Levels: SGB, species, genus

	Taxa_df <- utils::read.delim(in.file, sep = "\t", check.names = F, skip = skip.rows)
	Taxa_df$clade_name <- as.character(Taxa_df$clade_name)

	if (tax_level == "SGB") {
		Taxa_df <- Taxa_df[grepl("[|]t__",Taxa_df$clade_name),]
		Select.Taxa.level <- sapply(Taxa_df$clade_name, function(x){
				vect <- unlist(strsplit(x, "[|]"))
				TaxLevel <- grepl("t__",vect[length(vect)])	
				return(TaxLevel)})	
	}		
	if (tax_level == "species") {
		Taxa_df <- Taxa_df[grepl("[|]s__",Taxa_df$clade_name),]
		Select.Taxa.level <- sapply(Taxa_df$clade_name, function(x){
				vect <- unlist(strsplit(x, "[|]"))
				TaxLevel <- grepl("s__",vect[length(vect)])	
				return(TaxLevel)})	
	}
	if (tax_level == "genus") {
		Taxa_df <- Taxa_df[grepl("[|]g__",Taxa_df$clade_name),]
		Select.Taxa.level <- sapply(Taxa_df$clade_name, function(x){
				vect <- unlist(strsplit(x, "[|]"))
				TaxLevel <- grepl("g__",vect[length(vect)])	
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




