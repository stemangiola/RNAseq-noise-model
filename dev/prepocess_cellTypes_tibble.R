library(tidyverse)
library(magrittr)
library(foreach)
library(doParallel)
registerDoParallel()

check_df_FANTOM5 = function(){
	# FANTOM
	
	library(readr)
	library(GSEABase)
	library(edgeR)
	library(reshape2)
	
	
	read_delim(
		"~/PhD/deconvolution/FANTOM5/hg19.cage_peak_phase1and2combined_counts_ann.osc.txt.uncommented", 
		"\t", 
		escape_double = FALSE, 
		trim_ws = TRUE
	) %>%
	filter(grepl("p1@", short_description)) %>%
	gather(sample, `read count`, -c(1:7)) %>%
	separate(sample,  c("dummy", "info", "cnhs", "onto_link"), sep="\\.", remove = F) %>%
	separate(info, sprintf("info_%s", 1:20), sep="%20|%2c", remove = F) %>%	
		
	# filter for non induced
	filter(
		!grepl("day[0-9]+|[0-9]+hr", info) |
			grepl("00hr00min|day00", info)
	) %>%
	
	# Get gene symbol
	separate(short_description, c("dummy", "symbol"), remove = F) %>%
		
	# Attach ontology
	{
	
		onto = ontologyIndex::get_ontology(
			"~/PhD/deconvolution/FANTOM5/ff-phase2-140729.obo.txt", 
			propagate_relationships = "is_a", 
			extract_tags = "minimal"
		)

		(.) %>% left_join(
			
			(.) %>%
				distinct(onto_link) %>%
				rowwise() %>%
				do(
					ontologyIndex::get_term_property(
						ontology=onto, 
						property="ancestors", 
						term=sprintf("FF:%s", .$onto_link), 
						as_names=TRUE
					) %>% 
						as.data.frame() %>% 
						t() %>% 
						as_tibble() %>%
						mutate(onto_link =  names(.)[length(names(.))] %>% gsub("FF:", "", .) ) %>%
						setNames(
							c(
								names(.)[-c(length(names(.)), length(names(.))-1)], 
								c("FF:main_info", "onto_link")
								)
						) %>%
							gather(onto_category, onto_value, -onto_link)
				)  %>%
				
				# Filter onyl human samples
				right_join((.) %>% filter(onto_value == "human sample") %>% distinct(onto_link)  ) %>%
				
				# Establish cell types
				mutate(`Cell type` = onto_value) %>%
				mutate(`Cell type formatted` = NA) %>%
				mutate(`Cell type formatted` = ifelse(grepl("fibrobl", onto_value, ignore.case = T), "fibroblast", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("epithelial", onto_value, ignore.case = T), "epithelial", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("human CD14-positive CD16-negative Monocytes sample", onto_value, ignore.case = T), "monocyte", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("blood", onto_value, ignore.case = T), "epithelial", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("natural_killer", onto_value, ignore.case = T), "natural_killer", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("endothelial", onto_value, ignore.case = T) & !grepl("progenitor", onto_value, ignore.case = T), "endothelial", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("mature adipocyte", onto_value, ignore.case = T), "adipocyte", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("adipocyte", onto_value, ignore.case = T), "adipocyte", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("smooth muscle", onto_value, ignore.case = T), "smooth_muscle", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("brain", onto_value, ignore.case = T), "neural", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("astrocyte", onto_value, ignore.case = T), "neural", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("neuron", onto_value, ignore.case = T), "neural", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("neural", onto_value, ignore.case = T), "neural", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("osteoblast", onto_value, ignore.case = T), "osteoblast", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("skeletal muscle", onto_value, ignore.case = T), "skeletal_muscle", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("mesenchymal", onto_value, ignore.case = T), "stem_cell", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("stem", onto_value, ignore.case = T), "stem_cell", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("CD19-positive B cell", onto_value, ignore.case = T), "b_cell", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("CD8-positive T cell", onto_value, ignore.case = T), "t_CD8", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("CD4-positive T cell", onto_value, ignore.case = T), "t_CD4", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("neutrophil", onto_value, ignore.case = T), "neutrophil", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("human eosinophil sample", onto_value, ignore.case = T), "eosinophil", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("human natural killer", onto_value, ignore.case = T), "natural_killer", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("naive conventional T cells", onto_value, ignore.case = T), "t_cell_naive", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("memory conventional T cells", onto_value, ignore.case = T), "t_memory", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("Dendritic Cells - monocyte immature derived", onto_value, ignore.case = T), "dendritic", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("stem", onto_value, ignore.case = T), "stem_cell", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("stem", onto_value, ignore.case = T), "stem_cell", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("stem", onto_value, ignore.case = T), "stem_cell", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("stem", onto_value, ignore.case = T), "stem_cell", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("stem", onto_value, ignore.case = T), "stem_cell", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("stem", onto_value, ignore.case = T), "stem_cell", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("stem", onto_value, ignore.case = T), "stem_cell", `Cell type formatted`)) %>%
				mutate(`Cell type formatted` = ifelse(grepl("stem", onto_value, ignore.case = T), "stem_cell", `Cell type formatted`)) %>%
				
				# filter only recognised cell types
				filter(!is.na(`Cell type formatted`)) %>%
				distinct(onto_link, `Cell type`, `Cell type formatted`)
		
		)	

	}  %>% 
 	dplyr::distinct(sample, symbol, `Cell type`, `Cell type formatted`, `read count`, entrezgene_id) %>%
		
	mutate(`Data base` = "FANTOM5")

}

check_df_BLUEPRINT = function(){
	# Parse BLUEPRINT data base
	

	read_delim(
		"/wehisan/home/allstaff/m/mangiola.s/PhD/deconvolution/BLUEPRINT_db/blueprint_files.tsv", 
		"\t", 
		escape_double = FALSE, 
		trim_ws = TRUE
	) %>%
	filter(`File type` == "Transcription quantification (Genes)") %>%
	filter(!is.na(`Cell type`)) %>%
	
	# Get data from URL
	separate(URL, sep="/|\\.", sprintf("URL_%s", 1:30)) %>%
	mutate(`Cell type` = gsub("_", " ", URL_15)) %>%
	mutate(sample = URL_22) %>%
	
	# Chenage cell type names
	mutate(`Cell type formatted` = NA) %>%
	mutate(`Cell type formatted` = ifelse(grepl("plasma cell", `Cell type`, ignore.case=T), "plasma_cell", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("plasma cell", `Cell type`, ignore.case=T), "plasma_cell", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("band form neutrophil", `Cell type`, ignore.case=T), "neutrophil", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("mature neutrophil", `Cell type`, ignore.case=T), "neutrophil", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("segmented neutrophil of bone marrow", `Cell type`, ignore.case=T), "neutrophil", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("myeloid cell", `Cell type`, ignore.case=T), "myeloid", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("lymphocyte of B lineage", `Cell type`, ignore.case=T), "b_cell", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("CD14-positive, CD16-negative classical monocyte", `Cell type`, ignore.case=T), "monocyte", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("common lymphoid progenitor", `Cell type`, ignore.case=T), "lymphoid", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("hematopoietic multipotent progenitor cell", `Cell type`, ignore.case=T), "stem_cell", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("hematopoietic stem cell", `Cell type`, ignore.case=T), "stem_cell", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("CD38-negative naive B cell", `Cell type`, ignore.case=T), "b_cell_naive", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("CD8-positive, alpha-beta T cell", `Cell type`, ignore.case=T), "t_CD8", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("cytotoxic CD56-dim natural killer cell", `Cell type`, ignore.case=T), "natural_killer", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("common myeloid progenitor", `Cell type`, ignore.case=T), "myeloid", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("CD4-positive, alpha-beta T cell", `Cell type`, ignore.case=T), "t_CD4", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("inflammatory macrophage", `Cell type`, ignore.case=T), "macrophage_M1", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("macrophage", `Cell type`, ignore.case=T), "macrophage", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("endothelial", `Cell type`, ignore.case=T), "endothelial", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("alternatively activated macrophage", `Cell type`, ignore.case=T), "macrophage_M2", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("conventional dendritic cell", `Cell type`, ignore.case=T), "dendritic", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("germinal center B cell", `Cell type`, ignore.case=T), "b_cell_germinal_center", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("naive B cell", `Cell type`, ignore.case=T), "b_cell_naive", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("immature conventional dendritic cell", `Cell type`, ignore.case=T), "dendritic_immature", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("mature conventional dendritic cell", `Cell type`, ignore.case=T), "dendritic", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("osteoclast", `Cell type`, ignore.case=T), "osteoclast", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("central memory CD4-positive, alpha-beta T cell", `Cell type`, ignore.case=T), "t_memory_central", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("effector memory CD4-positive, alpha-beta T cell", `Cell type`, ignore.case=T), "t_memory_effector", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("regulatory T cell", `Cell type`, ignore.case=T), "t_reg", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("central memory CD8-positive, alpha-beta T cell", `Cell type`, ignore.case=T), "t_CD8_memory_central", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("effector memory CD8-positive, alpha-beta T cell", `Cell type`, ignore.case=T), "t_CD8_memory_effector", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("class switched memory B cell", `Cell type`, ignore.case=T), "b_cell_memory", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("memory B cell", `Cell type`, ignore.case=T), "b_cell_memory", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("monocyte", `Cell type`, ignore.case=T), "monocyte", `Cell type formatted`)) %>%	
	mutate(`Cell type formatted` = ifelse(grepl("peripheral blood mononuclear cell", `Cell type`, ignore.case=T), "mono_derived", `Cell type formatted`)) %>%	

	# Filter 
	filter(!is.na(`Cell type formatted`)) %>%
	
	# Select info
	select(sample, `Cell type`, `Cell type formatted`) %>%
	
	# Add expression
	left_join(
		
		foreach(
			my_file = dir(
				path="/wehisan/home/allstaff/m/mangiola.s/PhD/deconvolution/BLUEPRINT_db/", 
				pattern="results", 
				full.names=T
			), 
			.combine = bind_rows
		) %dopar% {
			read_delim(
				my_file,
				"\t", 
				escape_double = FALSE, 
				trim_ws = TRUE
			) %>% 
					mutate(sample = my_file %>% basename())
		} %>%
		separate(sample, sep="\\.", sprintf("sample_%s", 1:6)) %>%
		mutate(sample = sample_5) %>%
		select(gene_id, expected_count, sample)
		
	) %>%
		
		rename(`read count` =  expected_count) %>%
		mutate(`Data base` = "BLUEPRINT")
	

}

check_df_bloodRNA = function(){
	# Parse BLUEPRINT data base
	
	source("~/PhD/deconvolution/ARMET_dev/lib/ARMET_utilities.R")
	library(readr)
	
	load("~/PhD/deconvolution/test_ARMET/myTest_TME_first_run_pure_populations_TME/humanRNASeqCnts_CdG160630.rda")

	bl.count = counts$counts
	
	groups = (seq(ncol(bl.count))-1) %/% 2
	temp =                     data.frame(ct = groups, t(bl.count), row.names=NULL)
	df.mean =                  do.call("cbind", by(temp, temp$ct, function(df) colMedians(as.matrix(df[,-1]), na.rm=T)))
	rownames(df.mean) =        rownames(bl.count)
	colnames(df.mean) = 			colnames(bl.count)[unique((groups+1)*2)]
	bl.count = df.mean
	
	info = as_tibble(counts$samples)
	
	unformatted_cell_types = colnames(bl.count)
	
	formatted_cell_types = c( "dendritic"  ,          "dendritic"   ,         "dendritic", "t_CD4"     ,
														"t_CD8"    ,         "eosinophil",     "dendritic"  ,          
														 "b_cell_memory",          "monocyte",       "b_cell_naive"        ,   "neutrophil",  
														"natural_killer",         "dendritic_plasmocytoid"        ,     "dendritic",
														 "t_CD4"    ,         "t_CD8"       ,      "eosinophil"  ,   "mDCs"      ,   
														"b_cell_memory"  ,        "monocyte"   ,    "b_cell_naive"          ,
														 "natural_killer" ,        "t_CD4"        ,     "t_CD8"       ,      "Mem"        ,  
														"monocyte",       "b_cell_naive"        ,   "neutrophil"    ,
														 "natural_killer"  ,       "t_CD4"         ,    "t_CD8"       ,      "monocyte"   , 
														"b_cell_naive"     ,      "natural_killer"       ,  "t_CD4"            ,
													 "t_CD8"         ,    "mDC"            , "monocyte"   ,    "b_cell_naive"          ,
														"neutrophil"  ,   "natural_killer"        , "dendritic_plasmocytoid"  )
		
	
	colnames(bl.count) = formatted_cell_types
	
	library('biomaRt')
	#mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))
	#save(mart, file=sprintf("%s/PhD/deconvolution/mart_ENSAMBL_to_SYMBOL.RData", Sys.getenv( "HOME" )))
	load(sprintf("%s/PhD/deconvolution/mart_ENSAMBL_to_SYMBOL.RData", Sys.getenv( "HOME" )))
	G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=rownames(bl.count),mart= mart)
	G_list = G_list[G_list$hgnc_symbol!="",]
	bl.count = bl.count[G_list$ensembl_gene_id,]
	rownames(bl.count) = G_list$hgnc_symbol
	bl.count = bl.count[-duplicated(rownames(bl.count)),]
	
	return(bl.count)
}

check_df_ENCODE = function(){


	metadata <- read_delim("/wehisan/home/allstaff/m/mangiola.s/PhD/deconvolution/ENCODE/metadata.tsv",  "\t", escape_double = FALSE, trim_ws = TRUE)
	metadata = metadata %>% filter(`Output type` == "gene quantifications") %>% filter(Assembly == "hg19")
	metadata$`Cell type formatted` = NA
	metadata$`Cell type formatted`[grep("epithelial", metadata$`Biosample term name`, ignore.case=T)] = "epithelial"
	metadata$`Cell type formatted`[grep("stem cell", metadata$`Biosample term name`, fixed=T) ] = "stem_cell"
	metadata$`Cell type formatted`[grep("fibroblast", metadata$`Biosample term name`, fixed=T) ] = "fibroblast"
	metadata$`Cell type formatted`[grep("endothelial", metadata$`Biosample term name`, fixed=T) ] = "endothelial"
	metadata$`Cell type formatted`[grep("smooth muscle", metadata$`Biosample term name`, fixed=T) ] = "smooth_muscle"
	metadata$`Cell type formatted`[grep("neuron", metadata$`Biosample term name`, fixed=T) ] = "neural"
	metadata$`Cell type formatted`[grep("keratinocyte", metadata$`Biosample term name`, fixed=T) ] = "keratinocyte"
	metadata$`Cell type formatted`[grep("astrocyte", metadata$`Biosample term name`, fixed=T) ] = "astrocyte"
	metadata$`Cell type formatted`[grep("dendritic cell", metadata$`Biosample term name`, fixed=T) ] = "dendritic"
	metadata$`Cell type formatted`[grep("myocyte", metadata$`Biosample term name`, fixed=T) ] = "myocyte"
	metadata$`Cell type formatted`[grep("natural killer", metadata$`Biosample term name`, fixed=T) ] = "natural_killer"
	metadata$`Cell type formatted`[grep("CD4-positive helper T cell", metadata$`Biosample term name`, fixed=T) ] = "t_helper"
	metadata$`Cell type formatted`[grep("neural", metadata$`Biosample term name`, fixed=T) ] = "neural"
	metadata$`Cell type formatted`[grep("CD14-positive monocyte", metadata$`Biosample term name`, fixed=T) ] = "mono_derived"
	metadata$`Cell type formatted`[grep("hepatocyte", metadata$`Biosample term name`, fixed=T) ] = "hepatocyte"
	metadata$`Cell type formatted`[grep("hematopoietic multipotent progenitor cell", metadata$`Biosample term name`, fixed=T) ] = "stem_cell"
	metadata$`Cell type formatted`[grep("chondrocyte", metadata$`Biosample term name`)] = "chondrocyte"
	metadata$`Cell type formatted`[grep("CD8-positive, alpha-beta T cell", metadata$`Biosample term name`) ] = "t_CD8"
	metadata$`Cell type formatted`[grep("melanocyte", metadata$`Biosample term name`, fixed=T) ] = "melanocyte"
	metadata$`Cell type formatted`[grep("B cell", metadata$`Biosample term name`, fixed=T) ] = "b_cell"
	metadata$`Cell type formatted`[grep("osteoblast", metadata$`Biosample term name`, fixed=T) ] = "osteoblast"
	metadata$`Cell type formatted`[grep("naive B cell", metadata$`Biosample term name`, fixed=T) ] = "b_cell_naive"
	metadata$`Cell type formatted`[grep("T-cell", metadata$`Biosample term name`, fixed=T) ] = "t_cell"
	metadata = metadata[!is.na(metadata$`Cell type formatted`),]
	
	metadata %>%
	dplyr::mutate(sample = `File accession`) %>%
	mutate(`Cell type` = `Biosample term name`) %>%
		
	left_join(
		# Get counts from files
		foreach(f = metadata$`File accession`, .combine = bind_rows) %dopar% {
			read_delim(
				sprintf("/wehisan/home/allstaff/m/mangiola.s/PhD/deconvolution/ENCODE/%s.tsv", f), 
				"\t", 
				escape_double = FALSE, 
				trim_ws = TRUE
			) %>% 
				dplyr:::select(gene_id, expected_count) %>%
				mutate(sample = f)
		} %>%
			
		# Add symbol
		spread(sample, expected_count) %>%
		dplyr::rename(raw_geneID = gene_id) %>%
		separate(raw_geneID, c("ENSEMBL_ID", "dummy"), sep = "\\." ) %>%
		dplyr::select(-dummy) %>%
		mutate(
			symbol = 
				AnnotationDbi::mapIds(
					org.Hs.eg.db::org.Hs.eg.db,
					keys=ENSEMBL_ID,
					column="SYMBOL",
					keytype="ENSEMBL",
					multiVals="first"
				)
		) %>%
		gather(sample, `read count`, -ENSEMBL_ID, -symbol)
	) %>%
		
		dplyr::select(sample, `Cell type`, `Cell type formatted`, `read count`, symbol, ENSEMBL_ID) %>%
		mutate(`Data base` = "ENCODE")

}

integrate_RNAseq = function(){
	
	
	
	xx = 
		check_df_ENCODE() %>%
		bind_rows(check_df_BLUEPRINT()) %>%
		bind_rows(check_df_FANTOM5()) 
		
	xx %>% 
		group_by(`Cell type formatted`) %>%
		do({
			
			(.) %>%
				write_csv(sprintf("tibble_cellTypes_%s.csv", (.) %>% pull(`Cell type formatted`) %>% unique %>% { gsub("/| |,", "_", (.)) } ))
			
			tibble(`Cell type` = (.) %>% pull(`Cell type formatted`) %>% unique)
			
		})

}
