library(tidyverse)
library(magrittr)
library(foreach)
library(doParallel)
registerDoParallel()

get_FANTOM5 = function(){

  onto = ontologyIndex::get_ontology(
    "~/PhD/deconvolution/FANTOM5/ff-phase2-140729.obo.txt",
    propagate_relationships = "is_a",
    extract_tags = "minimal"
  )

  reads = read_delim(
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
    separate(short_description, c("dummy", "symbol"), sep="@", remove = F)

  # reads  %>%
  #   left_join(
  #
  #     (.) %>%
  #       distinct(onto_link) %>%
  #       rowwise() %>%
  #       do(
  #         ontologyIndex::get_term_property(
  #           ontology=onto,
  #           property="ancestors",
  #           term=sprintf("FF:%s", .$onto_link),
  #           as_names=TRUE
  #         ) %>%
  #           as.data.frame() %>%
  #           t() %>%
  #           as_tibble() %>%
  #           mutate(onto_link =  names(.)[length(names(.))] %>% gsub("FF:", "", .) ) %>%
  #           setNames(
  #             c(
  #               names(.)[-c(length(names(.)), length(names(.))-1)],
  #               c("FF:main_info", "onto_link")
  #             )
  #           ) %>%
  #           gather(onto_category, onto_value, -onto_link)
  #       )  %>%
  #
  #       # Filter onyl human samples
  #       right_join((.) %>% filter(onto_value == "human sample") %>% distinct(onto_link)  ) %>%
  #
  #       # Filter only main info
  #       filter(onto_category == "FF:main_info") %>%
  #
  #       # Establish cell types
  #       mutate(`Cell type` = onto_value)
  #   ) %>%
  #   distinct(onto_link ,  sample, `Cell type`) %>%
  #   write_csv("big_data/tibble_cellType_files/FANTOM5_annotation_cell_types.csv")

  reads %>%
    # Attach ontology
    left_join(
      read_csv("~/PostDoc/RNAseq-noise-model/big_data/tibble_cellType_files/FANTOM5_annotation_cell_types.csv")
    ) %>%
    dplyr::distinct(sample, symbol, `Cell type`, `Cell type formatted`, `read count`, entrezgene_id) %>%

    mutate(`Data base` = "FANTOM5")

}

get_BLUEPRINT = function(){
  # Parse BLUEPRINT data base


  # read_delim(
  #   "/wehisan/home/allstaff/m/mangiola.s/PhD/deconvolution/BLUEPRINT_db/blueprint_files.tsv",
  #   "\t",
  #   escape_double = FALSE,
  #   trim_ws = TRUE
  # ) %>%
  #   filter(`File type` == "Transcription quantification (Genes)") %>%
  #   filter(!is.na(`Cell type`)) %>%
  #
  #   # Get data from URL
  #   separate(URL, sep="/|\\.", sprintf("URL_%s", 1:30), remove = F) %>%
  #   mutate(`Cell type` = gsub("_", " ", URL_15)) %>%
  #   mutate(sample = URL_18) %>%
  #
  #   distinct(Group, `Sub-group`, `Cell type`, Tissue, sample) %>%
  #   write_csv("big_data/tibble_cellType_files/BLUEPRINT__annotation_cell_types.csv")

    # # Chenage cell type names
    # mutate(`Cell type formatted` = NA) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("plasma cell", `Cell type`, ignore.case=T), "plasma_cell", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("plasma cell", `Cell type`, ignore.case=T), "plasma_cell", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("band form neutrophil", `Cell type`, ignore.case=T), "neutrophil", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("mature neutrophil", `Cell type`, ignore.case=T), "neutrophil", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("segmented neutrophil of bone marrow", `Cell type`, ignore.case=T), "neutrophil", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("myeloid cell", `Cell type`, ignore.case=T), "myeloid", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("lymphocyte of B lineage", `Cell type`, ignore.case=T), "b_cell", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("CD14-positive, CD16-negative classical monocyte", `Cell type`, ignore.case=T), "monocyte", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("common lymphoid progenitor", `Cell type`, ignore.case=T), "lymphoid", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("hematopoietic multipotent progenitor cell", `Cell type`, ignore.case=T), "stem_cell", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("hematopoietic stem cell", `Cell type`, ignore.case=T), "stem_cell", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("CD38-negative naive B cell", `Cell type`, ignore.case=T), "b_cell_naive", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("cytotoxic CD56-dim natural killer cell", `Cell type`, ignore.case=T), "natural_killer", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("common myeloid progenitor", `Cell type`, ignore.case=T), "myeloid", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("inflammatory macrophage", `Cell type`, ignore.case=T), "macrophage_M1", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("macrophage", `Cell type`, ignore.case=T), "macrophage", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("endothelial", `Cell type`, ignore.case=T), "endothelial", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("alternatively activated macrophage", `Cell type`, ignore.case=T), "macrophage_M2", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("conventional dendritic cell", `Cell type`, ignore.case=T), "dendritic", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("germinal center B cell", `Cell type`, ignore.case=T), "b_cell_germinal_center", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("naive B cell", `Cell type`, ignore.case=T), "b_cell_naive", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("immature conventional dendritic cell", `Cell type`, ignore.case=T), "dendritic_immature", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("mature conventional dendritic cell", `Cell type`, ignore.case=T), "dendritic", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("osteoclast", `Cell type`, ignore.case=T), "osteoclast", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("class switched memory B cell", `Cell type`, ignore.case=T), "b_cell_memory", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("memory B cell", `Cell type`, ignore.case=T), "b_cell_memory", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("monocyte", `Cell type`, ignore.case=T), "monocyte", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("peripheral blood mononuclear cell", `Cell type`, ignore.case=T), "mono_derived", `Cell type formatted`)) %>%
    #
    # mutate(`Cell type formatted` = ifelse(`Cell type`==("CD8-positive alpha-beta T cell"), "t_CD8", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(`Cell type`==("CD4-positive alpha-beta T cell"), "t_CD4", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("central memory CD8-positive", `Cell type`, ignore.case=T), "t_CD8_memory_central", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("effector memory CD8-positive", `Cell type`, ignore.case=T), "t_CD8_memory_effector", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("central memory CD4-positive", `Cell type`, ignore.case=T), "t_memory_central", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("effector memory CD4-positive", `Cell type`, ignore.case=T), "t_memory_effector", `Cell type formatted`)) %>%
    # mutate(`Cell type formatted` = ifelse(grepl("regulatory T cell", `Cell type`, ignore.case=T), "t_reg", `Cell type formatted`)) %>%

    # Select info
    read_csv("big_data/tibble_cellType_files/BLUEPRINT__annotation_cell_types.csv") %>%

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
        separate(sample, sep="\\.", sprintf("sample_%s", 1:6), remove=F) %>%
        mutate(sample = sample_1) %>%
        select(gene_id, expected_count, sample)

    ) %>%

    rename(`read count` =  expected_count) %>%

    # Attach symbol names
    left_join(
      (.) %>%
        distinct(gene_id) %>%
        separate(gene_id, sep="\\.", c("ensembl_gene_id", "isoform"), remove=F)
      ) %>%
    left_join(
        # Query symbols
       (.) %>%
        distinct(ensembl_gene_id) %>%
        biomaRt::getBM(
          filters= "ensembl_gene_id",
          attributes= c("ensembl_gene_id","hgnc_symbol"),
          values=.$ensembl_gene_id,
          mart= biomaRt::useDataset(
            "hsapiens_gene_ensembl",
            biomaRt::useMart("ensembl")
          )
        ) %>%
        as_tibble() %>%
        rename(symbol = hgnc_symbol)
    ) %>%
    distinct %>%
    mutate(`Data base` = "BLUEPRINT")


}

get_bloodRNA = function(){

  load("~/PhD/deconvolution/test_ARMET/myTest_TME_first_run_pure_populations_TME/humanRNASeqCnts_CdG160630.rda")

  counts$counts %>%
    as_tibble(rownames="ensembl_gene_id") %>%
    gather(file, `read count`, -ensembl_gene_id) %>%
    mutate(`Cell type` = NA) %>%

    # Format cell types
    mutate(`Cell type` = ifelse(grepl("_mDC_", file), "dendritic_m", `Cell type`)) %>%
    mutate(`Cell type` = ifelse(grepl("_CD4_Tcells_", file), "t_CD4", `Cell type`)) %>%
    mutate(`Cell type` = ifelse(grepl("_Eosinophils_", file), "eosinophil", `Cell type`)) %>%
    mutate(`Cell type` = ifelse(grepl("_Memory_Bcells_", file), "b_memory", `Cell type`)) %>%
    mutate(`Cell type` = ifelse(grepl("_Monocytes_", file), "monocyte", `Cell type`)) %>%
    mutate(`Cell type` = ifelse(grepl("_Naive_Bcells_", file), "b_naive", `Cell type`)) %>%
    mutate(`Cell type` = ifelse(grepl("_Neutrophils_", file), "neutrophil", `Cell type`)) %>%
    mutate(`Cell type` = ifelse(grepl("_Nkcells_", file), "natural_killer", `Cell type`)) %>%
    mutate(`Cell type` = ifelse(grepl("_CD8_Tcells_", file), "t_CD8", `Cell type`)) %>%
    mutate(`Cell type` = ifelse(grepl("_Mem_Bcell_", file), "b_memory", `Cell type`)) %>%
    filter(`Cell type` %>% is.na %>% `!`) %>%

    # Merge same samples
    mutate(sample = gsub("_C1B73ACXX.+", "", file)) %>%
    group_by(sample, ensembl_gene_id, `Cell type`) %>%
    summarise(`read count` = `read count` %>% median(na.rm=T)) %>%
    ungroup() %>%


    # Cell type formatted
    mutate(`Cell type formatted` = `Cell type`) %>%

    # Attach symbol names
    left_join(
      (.) %>%
      distinct(ensembl_gene_id) %>%
      biomaRt::getBM(
        filters= "ensembl_gene_id",
        attributes= c("ensembl_gene_id","hgnc_symbol"),
        values=.$ensembl_gene_id,
        mart= biomaRt::useDataset(
          "hsapiens_gene_ensembl",
          biomaRt::useMart("ensembl")
        )
      ) %>%
        as_tibble()
    ) %>%
    rename(symbol = hgnc_symbol) %>%
    distinct() %>%
    mutate(`Data base` = "bloodRNA")

}

get_ENCODE = function(){


  # read_delim("/wehisan/home/allstaff/m/mangiola.s/PhD/deconvolution/ENCODE/metadata.tsv",  "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  #   filter(`Output type` == "gene quantifications") %>%
  #   dplyr::mutate(sample = `File accession`) %>%
  #   mutate(`Cell type` = `Biosample term name`) %>%
  #   distinct(sample, `Cell type`, `Biosample type`) %>%
  #   write_csv("big_data/tibble_cellType_files/ENCODE__annotation_cell_types.csv")
  #
  #
  # metadata$`Cell type formatted` = NA
  # metadata$`Cell type formatted`[grep("epithelial", metadata$`Biosample term name`, ignore.case=T)] = "epithelial"
  # metadata$`Cell type formatted`[grep("stem cell", metadata$`Biosample term name`, fixed=T) ] = "stem_cell"
  # metadata$`Cell type formatted`[grep("fibroblast", metadata$`Biosample term name`, fixed=T) ] = "fibroblast"
  # metadata$`Cell type formatted`[grep("endothelial", metadata$`Biosample term name`, fixed=T) ] = "endothelial"
  # metadata$`Cell type formatted`[grep("smooth muscle", metadata$`Biosample term name`, fixed=T) ] = "smooth_muscle"
  # metadata$`Cell type formatted`[grep("neuron", metadata$`Biosample term name`, fixed=T) ] = "neural"
  # metadata$`Cell type formatted`[grep("keratinocyte", metadata$`Biosample term name`, fixed=T) ] = "keratinocyte"
  # metadata$`Cell type formatted`[grep("astrocyte", metadata$`Biosample term name`, fixed=T) ] = "astrocyte"
  # metadata$`Cell type formatted`[grep("dendritic cell", metadata$`Biosample term name`, fixed=T) ] = "dendritic"
  # metadata$`Cell type formatted`[grep("myocyte", metadata$`Biosample term name`, fixed=T) ] = "myocyte"
  # metadata$`Cell type formatted`[grep("natural killer", metadata$`Biosample term name`, fixed=T) ] = "natural_killer"
  # metadata$`Cell type formatted`[grep("CD4-positive helper T cell", metadata$`Biosample term name`, fixed=T) ] = "t_helper"
  # metadata$`Cell type formatted`[grep("neural", metadata$`Biosample term name`, fixed=T) ] = "neural"
  # metadata$`Cell type formatted`[grep("CD14-positive monocyte", metadata$`Biosample term name`, fixed=T) ] = "mono_derived"
  # metadata$`Cell type formatted`[grep("hepatocyte", metadata$`Biosample term name`, fixed=T) ] = "hepatocyte"
  # metadata$`Cell type formatted`[grep("hematopoietic multipotent progenitor cell", metadata$`Biosample term name`, fixed=T) ] = "stem_cell"
  # metadata$`Cell type formatted`[grep("chondrocyte", metadata$`Biosample term name`)] = "chondrocyte"
  # metadata$`Cell type formatted`[grep("CD8-positive, alpha-beta T cell", metadata$`Biosample term name`) ] = "t_CD8"
  # metadata$`Cell type formatted`[grep("melanocyte", metadata$`Biosample term name`, fixed=T) ] = "melanocyte"
  # metadata$`Cell type formatted`[grep("B cell", metadata$`Biosample term name`, fixed=T) ] = "b_cell"
  # metadata$`Cell type formatted`[grep("osteoblast", metadata$`Biosample term name`, fixed=T) ] = "osteoblast"
  # metadata$`Cell type formatted`[grep("naive B cell", metadata$`Biosample term name`, fixed=T) ] = "b_cell_naive"
  # metadata$`Cell type formatted`[grep("T-cell", metadata$`Biosample term name`, fixed=T) ] = "t_cell"
  # metadata = metadata[!is.na(metadata$`Cell type formatted`),]

  read_csv("big_data/tibble_cellType_files/ENCODE__annotation_cell_types.csv") %>%

    left_join(
      # Get counts from files
      foreach(f = .$sample, .combine = bind_rows) %dopar% {
        read_delim(
          sprintf("/wehisan/home/allstaff/m/mangiola.s/PhD/deconvolution/ENCODE/%s.tsv", f),
          "\t",
          escape_double = FALSE,
          trim_ws = TRUE
        ) %>%
          dplyr:::select(gene_id, expected_count) %>%
          mutate(sample = f)
      }
    ) %>%
    {
      samples_with_symbol = c("ENCFF060YNO", "ENCFF677SZA", "ENCFF708ZUJ", "ENCFF255ULI", "ENCFF717WSQ", "ENCFF118GPH", "ENCFF083PYO" ,"ENCFF712VOY" ,"ENCFF094ADI",
                              "ENCFF841AKS", "ENCFF331CDB", "ENCFF263OIE", "ENCFF798GKH" ,"ENCFF491YKJ" ,"ENCFF867RFN", "ENCFF461BKM", "ENCFF929RZY", "ENCFF440CJU",
                              "ENCFF246ZOR", "ENCFF680YEW")
      bind_rows(

        (.) %>%
          filter(sample %in% samples_with_symbol) %>%
          rename(symbol = gene_id),

        (.) %>% filter(!sample %in% samples_with_symbol) %>%
          left_join(
            (.) %>%
              separate(gene_id, c("ENSEMBL_ID", "dummy"), sep = "\\." , remove = F) %>%
              distinct(ENSEMBL_ID, gene_id) %>%
              mutate(
                symbol =
                  AnnotationDbi::mapIds(
                    org.Hs.eg.db::org.Hs.eg.db,
                    keys=ENSEMBL_ID,
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first"
                  ) %>%
                  unlist
              )
          )
      )
    } %>%
    rename(`read count` = expected_count) %>%
    dplyr::select(sample, `Cell type`, `Cell type formatted`, `read count`, symbol, ENSEMBL_ID) %>%
    mutate(`Data base` = "ENCODE")

}


ENCODE = get_ENCODE()
save(ENCODE, file="big_data/tibble_cellType_files/ENCODE.RData")

BLUEPRINT = get_BLUEPRINT()
save(BLUEPRINT, file="big_data/tibble_cellType_files/BLUEPRINT.RData")

FANTOM5 = get_FANTOM5()
save(FANTOM5, file="big_data/tibble_cellType_files/FANTOM5.RData")

bloodRNA = get_bloodRNA()
save(bloodRNA, file="big_data/tibble_cellType_files/bloodRNA.RData")

#     GSE115898
# DC, T cell

#  GSE107011
# 29 immune cell types

ENCODE %>%
  bind_rows(BLUEPRINT) %>%
  bind_rows(FANTOM5) %>%
  bind_rows(bloodRNA) %>%
  filter(symbol %>% is.na %>% `!`) %>%
  filter(`Cell type formatted` %>% is.na %>% `!`) %>%
  mutate(`read count` = `read count` %>% as.integer) %>%
  group_by(`Cell type formatted`) %>%
  do({

    (.) %>%
      write_csv(
        sprintf(
          "big_data/tibble_cellType_files/tibble_cellTypes_%s.csv",
          (.) %>%
            pull(`Cell type formatted`) %>%
            unique %>%
            { gsub("/| |,", "_", (.)) }
        )
      )

    tibble(`Cell type` = (.) %>% pull(`Cell type formatted`) %>% unique)

  })


