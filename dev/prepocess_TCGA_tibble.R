setwd("~/PhD/deconvolution/TCGA_all_cancers")
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")

library(tidyverse)

load("tcga.ex.RData")

# Collect annotation of genes
conversion_ids = biomaRt::getBM(
	attributes = c(
		'ensembl_gene_id', 
		'entrezgene','hgnc_symbol', 
		'chromosome_name',
		'start_position',
		'end_position',
		'strand',
		'transcript_length',
		'percentage_gene_gc_content'
	),
	mart = biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
) %>%
	as_tibble() %>%
	mutate(ens = ensembl_gene_id, symbol=hgnc_symbol, entrez = entrezgene) 

# Annotation for correct ID
convert_name_df = 
	read_delim("TCGA_sample_annotation.tab", "\t") %>%
	separate(cases_0_samples_0_portions_0_analytes_0_aliquots_0_submitter_id, c("t1", "t2", "t3", "t4"), sep = "-") %>%
	unite("sample_armonised", c("t1", "t2", "t3", "t4"), sep = "-") 

# Annotation for samples
sample_info = 
	read_delim("TCGA_sample_annotation.tab", "\t") %>%
	separate(cases_0_samples_0_portions_0_analytes_0_aliquots_0_submitter_id, c("t1", "t2", "t3", "t4"), sep = "-") %>%
	unite("sample_armonised", c("t1", "t2", "t3", "t4"), sep = "-") %>%
	mutate(sample = file_id)

# Fformat data set
tcga.ex %>% 
	apply(2, as.integer) %>% 
	as_tibble() %>%
	mutate(ens_iso = rownames(tcga.ex)) %>% 
	gather(sample, `read count`, -ens_iso) %>%
	
	# Add info
	left_join( sample_info %>% select("sample_armonised", "sample")) %>%
	left_join( convert_name_df %>% select("cases_0_project_disease_type", "cases_0_samples_0_sample_type", "sample_armonised") %>% distinct()) %>%
	select(-sample) %>%
	rename(sample = sample_armonised) %>%
	
	# Write to file
	group_by(cases_0_project_disease_type, cases_0_samples_0_sample_type) %>% 
	do(
		{
			data_set__name = (.)[1,] %>% 
				select(cases_0_project_disease_type, cases_0_samples_0_sample_type) %>% 
				as.character()  %>% 
				paste(collapse=" ") %>%
				gsub(" ", "_", .)
			
			write_csv((.), 	sprintf("tibble_TCGA_files/%s.csv", data_set__name))
			
			tibble( data_set = data_set__name, saved = "yes" )
		}
	)


