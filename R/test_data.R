load_test_data <- function(N, G, counts_source) {
  # Evaluate whether to recompute input data frame
  temp_file_counts = here::here("local_temp_data",paste0("subset_", counts_source,".rds"))
  if(!file.exists(temp_file_counts)) {
    recompute_counts_tidy = TRUE
  } else {
    counts_data <- readRDS(temp_file_counts)
    if(counts_data$N != N || counts_data$G != G) {
      recompute_counts_tidy = TRUE
    } else {
      counts_tidy = counts_data$counts_tidy
      recompute_counts_tidy = FALSE
    }
  }

  if(recompute_counts_tidy) {
    # Import gene transcription data set
    counts_tidy =
      read_csv(here("big_data", "tibble_TCGA_files", paste0(counts_source, ".csv"))) %>%
      filter(sample %in% ( (.) %>% distinct(sample) %>% head(n=N) %>% pull(sample)) ) %>%
      filter(ens_iso %in% ( (.) %>% distinct(ens_iso) %>% head(n=G) %>% pull(ens_iso)) ) %>%
      mutate_if(is.character, as.factor) %>%

      # Add numerical indexes
      left_join(
        (.) %>%
          distinct(ens_iso) %>%
          arrange(ens_iso) %>%
          mutate(ens_iso_idx = 1:n())
      ) %>%
      left_join(
        (.) %>%
          distinct(sample) %>%
          arrange(sample) %>%
          mutate(sample_idx = 1:n())
      )

    counts_data = list(N = N, G = G, counts_tidy = counts_tidy)
    saveRDS(counts_data, file = temp_file_counts)
  }

  counts_tidy
}


data_for_stan_from_tidy <- function(N, G, counts_tidy) {
  list(
    N = N,
    G = G,
    counts =
      counts_tidy %>%
      select(ens_iso,`read count`,sample) %>%
      spread(ens_iso, `read count`) %>%
      select(-sample) %>%
      as.matrix(),

    # Info on data set
    my_prior = c(0,5),
    holdout = array(0,N), #Don't hold anything out
    generate_quantities = 1,
    generate_log_lik = 1,
    is_prior_asymetric = 0,
    exposure =
      counts_tidy %>%
      group_by(sample) %>%
      summarise(`Total counts` = sum(`read count`)) %>%
      pull(`Total counts`)
  )
}
