library(tidyverse)
library(magrittr)
library(foreach)
library(doParallel)
registerDoParallel()
source("https://gist.githubusercontent.com/stemangiola/90a528038b8c52b21f9cfa6bb186d583/raw/cdf04a5988ab44c10d8a05fa46f1e175d526b2de/tidyTranscriptionTools.R")
library(ttBulk)

forget_FANTOM5 = function(){

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
    separate(gene_id, sep="\\.", c("ensembl_gene_id", "isoform"), remove=F) %>%
    add_symbol_from_ensembl("ensembl_gene_id") %>%
    rename(symbol = hgnc_symbol) %>%
    distinct %>%
    mutate(`Data base` = "BLUEPRINT")


  #AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, "ENSG00000069712", 'SYMBOL', 'ENSEMBL')
}

get_bloodRNA = function(){

  load("~/PhD/deconvolution/test_ARMET/myTest_TME_first_run_pure_populations_TME/humanRNASeqCnts_CdG160630.rda")

  counts$counts %>%
    as_tibble(rownames="ensembl_gene_id") %>%
    gather(file, `read count`, -ensembl_gene_id) %>%
    mutate(`Cell type` = NA) %>%

    # Format cell types
    mutate(`Cell type` = ifelse(grepl("_mDC_", file), "dendritic_myeloid", `Cell type`)) %>%

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
    add_symbol_from_ensembl("ensembl_gene_id") %>%
    rename(symbol = hgnc_symbol) %>%
    distinct %>%

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
          separate(gene_id, sep="\\.", c("ensembl_gene_id", "isoform"), remove=F) %>%
          add_symbol_from_ensembl("ensembl_gene_id") %>%
          rename(symbol = hgnc_symbol)
      )
    } %>%
    rename(`read count` = expected_count) %>%
    dplyr::select(sample, `Cell type`, `Cell type formatted`, `read count`, symbol, ensembl_gene_id) %>%
    distinct %>%

    mutate(`Data base` = "ENCODE")

}

get_N52_TME = function(){
  load("/wehisan/home/allstaff/m/mangiola.s/PhD/TME/N52_plus_12_run/fastq/run_all/code_repository/input_parallel_TABI.RData")
  ex %>% gather(sample, `read count`, -symbol) %>%
    left_join(
      annot %>%
        left_join(
          tibble(
            `Cell type` = c("EpCAM", "CD90+ CD31-", "CD45+CD3+", "CD45+CD16+"),
            `Cell type formatted` = c("epithelial", "fibroblast", "t_cell", "myeloid"),
            cell_type_formatted = c("E", "F", "T", "M")
          )
        ) %>%
        distinct(file, `Cell type`, `Cell type formatted`) %>%
        rename(sample = file)
    ) %>%
    mutate(`Data base` = "N52 prostate")
}

get_immune_singapoor = function(){

  # Emailed the author  GSE107011
  load("/wehisan/home/allstaff/m/mangiola.s/PhD/deconvolution/immune_singapoor/Genes_TPM_counts.RData")

  kallisto_gene %$% counts %>%
    as_tibble(rownames = "gene_id") %>%
    gather(sample, `read count`, -gene_id) %>%
    mutate(`read count` = `read count` %>% as.integer) %>%

    # Add cell type
    left_join(
      kallisto_gene %$% samples %>% as_tibble %>%
        distinct(sample, group10) %>%
        rename(`Cell type` = group10) %>%
        mutate(`Cell type` = gsub('[\"]', '', `Cell type`) %>% trimws  ) %>%
        left_join(
          read_csv("/wehisan/home/allstaff/m/mangiola.s/PhD/deconvolution/immune_singapoor/cell_type_conversion.csv")
        )
    ) %>%

    # Add symbol
    separate(gene_id, sep="\\.", c("ensembl_gene_id", "isoform"), remove=F) %>%
    ttBulk::annotate_symbol(ensembl_gene_id) %>%
    rename(symbol = transcript) %>%
    distinct() %>%

    # Add data base label
    mutate(`Data base` = "Immune Singapoor")

}

get_dendritic_STIMULATED_NOT_GOOD = function(){
  read_csv("/wehisan/home/allstaff/m/mangiola.s/PhD/deconvolution/mDC_database/GSE89442_Mathan_et_al_RNAseq_data_raw_values.csv") %>%
    gather(sample, `read count`, -Symbol,-GeneID) %>%
    rename(symbol = Symbol) %>%
    filter(grepl("unstimulated", sample)) %>%
    mutate(`Cell type` = sample) %>%
    mutate(sample = paste("GSE89442", sample)) %>%
    mutate(`Cell type formatted` = ifelse(grepl("mDC", sample), "dendritic_myeloid", "dendritic_plasmacytoid")) %>%
    mutate(`Data base` = "dendritic GSE89442")
}

get_influenza_immune = function(){

  # PRJNA271578
  # Download fastq pair enda
  # cat SRR_Acc_List.txt | parallel --eta -j 10  '~/third_party_sofware/sratoolkit.2.9.6-1-centos_linux64/bin/fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip {}'
  # module load cutadapt; module load trimgalore; ls fastq/ | awk '{split($0,a,"_"); print a[1]}' | uniq | parallel --eta -j1 trim_galore --paired {}_pass_1.fastq.gz {}_pass_2.fastq.gz --fastqc > trim_galore.log
  # module load STAR; for i in $(ls *_val_*fq.gz | rev | cut -c 15- | rev | uniq); do mkdir -p alignment_hg38/$i; cd alignment_hg38/$i; STAR --genomeDir ~/third_party_sofware/hg38/karyotypic/ --readFilesIn ../../$i'_1_val_1.fq.gz' ../../$i'_2_val_2.fq.gz' --readFilesCommand zcat --genomeLoad LoadAndKeep --runThreadN 42 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 40000000000 --outReadsUnmapped Fastx; cd ../../; done

  read_csv("~/PhD/deconvolution/immune_influenza_RNAseq_PBMC/counts_paired_end.csv") %>%
    separate(sample, c("dummy", "sample"), sep=c("hg38.")) %>%
    separate(sample, c("sample", "dummy"), sep=c("_pass")) %>%
    select(-dummy) %>%
    inner_join(
      read_csv("~/PhD/deconvolution/immune_influenza_RNAseq_PBMC/annot.csv") %>%
        rename(sample = Run, `Cell type` = `source name`) %>%
        filter(`Cell type` != "PBMC") %>%
        filter(time == "0 d") %>%
        select(`Cell type`, sample)
    ) %>%
    left_join(
      tibble(
        `Cell type` = c("T cells", "NK cells",  "Neutrophils", "Monocytes", "myeloid DC", "B cells"),
        `Cell type formatted` = c("t_cell", "natural_killer",  "neutrophil", "monocyte", "dendritic_myeloid", "b_cell")
      )
    ) %>%
    mutate(`Data base` = "influenza immune")

}

get_immune_skin = function(){

  # PRJNA476386
  # Download fastq pair enda
  # cat SRR_Acc_List.txt | parallel --eta -j 10  '~/third_party_sofware/sratoolkit.2.9.6-1-centos_linux64/bin/fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip {}'
  # module load cutadapt; module load trimgalore; module load parallel;  ls fastq/ | awk '{split($0,a,"_"); print a[1]}' | uniq | parallel --eta -j1 trim_galore --paired {}_pass_1.fastq.gz {}_pass_2.fastq.gz --fastqc > trim_galore.log
  # for i in $(ls *_val_*fq.gz | rev | cut -c 15- | rev | uniq); do qsub -l nodes=1:ppn=24,mem=41gb,walltime=120:00:00 STAR_HPC.sh -F "$i"; done


  # DOING

  read_csv("~/PhD/deconvolution/immune_skin/counts_paired_end.csv") %>%
    separate(sample, c("sample", "dummy"), sep=c("\\.pass")) %>%
    mutate(sample = gsub("\\.", "", sample)) %>%
    inner_join(
      read_csv("~/PhD/deconvolution/immune_skin/annot.csv") %>%
        rename(sample = Run, `Cell type` = `cell type`) %>%
        select(`Cell type`, sample)
    ) %>%
    left_join(
      tibble(
        `Cell type` = c("CD8+ T cell", "Keratinocyte",  "Dendritic cell", "CD4+ T effector"),
        `Cell type formatted` = c("t_CD8", "keratinocyte",  "dendritic_myeloid", "t_CD4_effector")
      )
    ) %>%
    mutate(`Data base` = "skin immune")

}

# from here

get_macro = function(){

  # PRJNA339309
  # Download fastq pair enda
  # cat SRR_Acc_List.txt | parallel --eta -j 10  '~/third_party_sofware/sratoolkit.2.9.6-1-centos_linux64/bin/fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip {}'
  # module load parallel; module load trimgalore; module load cutadapt; ls ./ | awk '{split($0,a,"_"); print a[1]}' | uniq | parallel --eta -j1 trim_galore --paired {}_pass_1.fastq.gz {}_pass_2.fastq.gz --fastqc > trim_galore.log
  # for i in $(ls *_val_*fq.gz | rev | cut -c 15- | rev | uniq); do qsub -l nodes=1:ppn=24,mem=41gb,walltime=120:00:00 STAR_HPC.sh -F "$i"; done



  create_tt_from_bam_sam_bulk



}

get_myeloid_differentiation = function(){

  # Emailed the author GSM2084371
  0

}

get_monocytes_brain = function(){

  # Emailed authors
  # Data set GSE88888 from  doi: 10.1038/s41598-018-28986-7

  GSE88888 <-
    read_delim("PhD/deconvolution/immune_monocyte_2/raw.txt",  "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    separate(`Annotation/Divergence`, c("transcript"), extra = "drop") %>%
    gather(sample, count, -(1:8)) %>%
    mutate(sample = gsub("/gpfs/data01/glasslab/home/jtao/analysis/pd_analysis//tag_directories/", "", sample)) %>%
    filter(grepl("RNA_Control", sample))

  GSE88888 %>%
    ttBulk(sample, transcript, count) %>%
    aggregate_duplicates() %>%
    scale_abundance() %>%
    reduce_dimensions(method="MDS") %>%
    select(contains("Dim"), sample) %>%
  distinct %>%
  ggplot(aes(x = `Dim 1`, y = `Dim 2`)) + geom_point()


}


ENCODE = get_ENCODE()
save(ENCODE, file="big_data/tibble_cellType_files/ENCODE.RData")

BLUEPRINT = get_BLUEPRINT()
save(BLUEPRINT, file="big_data/tibble_cellType_files/BLUEPRINT.RData")

# FANTOM5 = get_FANTOM5()
# save(FANTOM5, file="big_data/tibble_cellType_files/FANTOM5.RData")

bloodRNA = get_bloodRNA()
save(bloodRNA, file="big_data/tibble_cellType_files/bloodRNA.RData")

immune_singapoor = get_immune_singapoor()
save(immune_singapoor, file="big_data/tibble_cellType_files/immune_singapoor.RData")

N52_TME = get_N52_TME()
save(N52_TME, file="big_data/tibble_cellType_files/N52_TME.RData")

# dendritic = get_dendritic()
# save(dendritic, file="big_data/tibble_cellType_files/dendritic.RData")

influenza_immune = get_influenza_immune()
save(influenza_immune, file="big_data/tibble_cellType_files/influenza_immune.RData")

immune_skin = get_immune_skin()
save(immune_skin, file="big_data/tibble_cellType_files/immune_skin.RData")

load("big_data/tibble_cellType_files/ENCODE.RData")
load("big_data/tibble_cellType_files/BLUEPRINT.RData")
load("big_data/tibble_cellType_files/bloodRNA.RData")
load("big_data/tibble_cellType_files/immune_singapoor.RData")
load("big_data/tibble_cellType_files/N52_TME.RData")
load("big_data/tibble_cellType_files/influenza_immune.RData")
load("big_data/tibble_cellType_files/immune_skin.RData")


ENCODE %>%
  bind_rows(BLUEPRINT) %>%
  bind_rows(bloodRNA) %>%
  bind_rows(immune_singapoor) %>%
  bind_rows(N52_TME) %>%
  bind_rows(influenza_immune) %>%
  bind_rows(immune_skin) %>%

  #print statistics
  {
    (.) %>% distinct(`Cell type formatted`, `Data base`, sample) %>% count(`Cell type formatted`, `Data base`)%>% drop_na %>% arrange(n %>% desc) %>% spread(`Data base`, n) %>% print(n=99)
    (.)
  }

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


# (
#   ARMET::ARMET_ref %>%
#     filter(level==3) %>%
#     distinct(`Data base`, sample, symbol, `Cell type category`, `read count normalised bayes`)  %>%
#      inner_join(rr$signatures[[3]] %>% distinct(symbol)) %>%
#     reduce_dimensions(sample, symbol, `read count normalised bayes`, method = "tSNE", .dims=2) %>%
#     select(contains("tSNE"), `Data base`, `Cell type category`) %>%
#     distinct %>%
#     ggplot(aes(x = `tSNE 1`, y = `tSNE 2`, color=`Cell type category`, shape=`Data base`)) + geom_point(size=2) +
#     theme_bw() +
#     theme(
#       panel.border = element_blank(),
#       axis.line = element_line(),
#       panel.grid.major = element_line(size = 0.2),
#       panel.grid.minor = element_line(size = 0.1),
#       text = element_text(size=12),
#       legend.position="bottom",
#       aspect.ratio=1,
#       #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#       strip.background = element_blank(),
#       axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
#       axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
#     )
# ) %>% plotly::ggplotly()



(
  ARMET::ARMET_ref %>%
    distinct(`Data base`, sample, symbol, `Cell type formatted`, `read count`)  %>%
    ttBulk(sample, symbol, `read count`) %>%
    aggregate_duplicates() %>%
    scale_abundance() %>%
    #inner_join(rr$signatures[[3]] %>% distinct(symbol)) %>%
    reduce_dimensions(sample, symbol, `read count normalised`, method = "tSNE", .dims=2) %>%
    select(contains("tSNE"), `Data base`, `Cell type formatted`) %>%
    distinct %>%
    ggplot(aes(x = `tSNE 1`, y = `tSNE 2`, color=`Cell type formatted`, shape=`Data base`)) + geom_point(size=2) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line = element_line(),
      panel.grid.major = element_line(size = 0.2),
      panel.grid.minor = element_line(size = 0.1),
      text = element_text(size=12),
      legend.position="bottom",
      aspect.ratio=1,
      #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.background = element_blank(),
      axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
      axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
    )
) %>% plotly::ggplotly()

counts_annot =
  counts %>%
  filter(`read count` %>% is.na %>% `!`) %>%
  filter(symbol %>% is.na %>% `!`) %>%
  left_join(
    ARMET::tree %>%
      ARMET::ToDataFrameTypeColFull(fill = F) %>%
      rename(`Cell type formatted` = level_5)
  ) %>%
  ttBulk(sample, symbol, `read count`) %>%
  aggregate_duplicates() %>%
  scale_abundance() %>%
  reduce_dimensions(method = "MDS" ) %>%
  reduce_dimensions(method = "tSNE")
