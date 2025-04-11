set.seed(12345)
library("phyloseq"); 


standardize_taxonomy <- function(phylo.obj){

	OTU = otu_table(phylo.obj)
	tax = tax_table(phylo.obj)	
	tax <- as(tax,"matrix")

	for(i in rownames(tax)){
		temp.line<-tax[i,]
		temp.line[1]<-paste("k",temp.line[1],sep="_")
		temp.line[2]<-paste("p",temp.line[2],sep="_")
		temp.line[3]<-paste("c",temp.line[3],sep="_")
		temp.line[4]<-paste("o",temp.line[4],sep="_")
		temp.line[5]<-paste("f",temp.line[5],sep="_")
		temp.line[6]<-paste("g",temp.line[6],sep="_")
		#temp.line[7]<-paste("s",temp.line[7],sep="_")
		#print(temp.line)
		#temp.line<- tax[1321,]
		#if( length(grep("T", is.na(temp.line))) > 0 ){ #### Looks for the NA (missing taxonomic labels)
		if(length(grep("_NA",  temp.line )) > 0 ){ #### Looks for the NA (missing taxonomic labels)
			na.ind<-grep("_NA",  temp.line )  ### Locates in the row where are such values
			na.lev<-min( na.ind )                   ### Locates the first NA taxonomic label
			last.tax<-temp.line[na.lev-1]           ### Takes the last-known taxonomic category
			###last.tax<-as.matrix(last.lev)[1]         
			temp.line[ na.ind ] <- rep( paste( "uc",last.tax, sep="_" ),length(temp.line[ na.ind ]) ) ### Add an uc_LAST_KNOWN_TAX_LEVEL to all the unkwon taxonomic levels, ie: 
			tax[i,]<-temp.line     ### Change the line for the one whith the new taxonomic info
		} else{
			tax[i,]<-temp.line     ### Change the line for the one whith the new taxonomic info
		}
			
	}

	OTU = otu_table(phylo.obj,taxa_are_rows = T)
	TAX = tax_table(tax[  match( rownames(OTU) , rownames(tax) ) ,])
#	SD = sample_data(phylo.obj)
#	SD = SD[  match( colnames(OTU) , rownames(SD) ) ,]
		
#	physeq = phyloseq(OTU, TAX, SD)
	physeq = phyloseq(OTU, TAX)
	return(physeq)

}


standardize_taxonomy_ASV <- function(phylo.obj){

	OTU = otu_table(phylo.obj)
	tax = tax_table(phylo.obj)	
	tax <- as(tax,"matrix")

	for(i in rownames(tax)){
		temp.line<-tax[i,]
		temp.line[1]<-paste("k",temp.line[1],sep="_")
		temp.line[2]<-paste("p",temp.line[2],sep="_")
		temp.line[3]<-paste("c",temp.line[3],sep="_")
		temp.line[4]<-paste("o",temp.line[4],sep="_")
		temp.line[5]<-paste("f",temp.line[5],sep="_")
		temp.line[6]<-paste("g",temp.line[6],sep="_")
		temp.line[7]<-paste("s",temp.line[7],sep="_")
		temp.line[8]<-paste("asv",temp.line[8],sep="_")						
		#temp.line[7]<-paste("s",temp.line[7],sep="_")
		#print(temp.line)
		#temp.line<- tax[1321,]
		#if( length(grep("T", is.na(temp.line))) > 0 ){ #### Looks for the NA (missing taxonomic labels)
		if(length(grep("_NA",  temp.line )) > 0 ){ #### Looks for the NA (missing taxonomic labels)
			na.ind<-grep("_NA",  temp.line )  ### Locates in the row where are such values
			na.lev<-min( na.ind )                   ### Locates the first NA taxonomic label
			last.tax<-temp.line[na.lev-1]           ### Takes the last-known taxonomic category
			###last.tax<-as.matrix(last.lev)[1]         
			temp.line[ na.ind ] <- rep( paste( "uc",last.tax, sep="_" ),length(temp.line[ na.ind ]) ) ### Add an uc_LAST_KNOWN_TAX_LEVEL to all the unkwon taxonomic levels, ie: 
			tax[i,]<-temp.line     ### Change the line for the one whith the new taxonomic info
		} else{
			tax[i,]<-temp.line     ### Change the line for the one whith the new taxonomic info
		}
			
	}

	tax[,8] <- paste0( c(tax[,8]), "_", 1:nrow(tax))
	OTU = otu_table(phylo.obj,taxa_are_rows = T)
	TAX = tax_table(tax[  match( rownames(OTU) , rownames(tax) ) ,])
#	SD = sample_data(phylo.obj)
#	SD = SD[  match( colnames(OTU) , rownames(SD) ) ,]
		
#	physeq = phyloseq(OTU, TAX, SD)
	physeq = phyloseq(OTU, TAX)
	return(physeq)

}

