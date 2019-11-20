prepare_kfold <- function(K, model_defs, base_data, seed = NULL) {
  model_names <- unique(model_defs$model_name)

  models =
    foreach(model_name = model_names) %do% {
      library(rstan) #On Windows, %dopar%-ed code does not share the main session
      library(here)
      rstan_options(auto_write = TRUE)
      stan_model(here("stan",sprintf("%s.stan", model_name)))
    } %>%
    setNames(model_names)

  model_defs <- model_defs %>% mutate(id = 1:n())
  model_defs_kfold = model_defs %>% crossing(fold = 1:K)

  if(!is.null(seed)) {
    set.seed(seed)
  }

  folds = kfold_split_random(K, base_data$N)
  data_list = list()
  for(i in 1:nrow(model_defs_kfold)) {
    data_list[[i]] = base_data
    data_list[[i]]$holdout = (folds == model_defs_kfold$fold[i])
  }

  if(!is.null(model_defs$adapt_delta)) {
    control_list = map(model_defs_kfold$adapt_delta, function(x) { list(adapt_delta =x)})
  } else {
    control_list = NULL
  }

  list(
    base_data = base_data,
    data_list = data_list,
    control_list = list,
    model_defs = model_defs,
    model_defs_kfold = model_defs_kfold,
    models_list = models[as.character(model_defs_kfold$model_name)],
    models_unique = models
  )

}

run_kfold <- function(kfold_def, name_for_cache) {
  cache_dir = here("local_temp_data",name_for_cache)
  if(!dir.exists(cache_dir)) {
    dir.create(cache_dir)
  }
  fits = sampling_multi(kfold_def$models_list,
                        kfold_def$data_list,
                        control_per_item = kfold_def$control_list,
                        cache_dir = cache_dir,
                        cores = getOption("mc.cores"))

  kfold_log_lik <- kfold_def$model_defs$id %>% map(function(id) {
    indices = kfold_def$model_defs_kfold$id == id
    fits_for_model <- fits[indices]
    holdout_for_model <- lapply(kfold_def$data_list[indices], '[[', "holdout")
    extract_log_lik_K(fits_for_model, holdout_for_model, "log_lik")
  })

  res_names <- kfold_def$model_defs %>% unite(full_name, -id) %>% pull(full_name)
  kfold_res <- kfold_log_lik %>%
    map(kfold) %>%
    set_names(res_names)

  list(
    fits = fits,
    log_lik = kfold_log_lik,
    kfold = kfold_res
  )
}

extract_holdout_ranks <- function(list_of_stanfits, list_of_holdout, base_data) {
  K <- length(list_of_stanfits)
  N_samples <- 100
  all_genes_holdout_gen <- array(NA_real_, c(N_samples, base_data$N, base_data$G))
  for(i in 1:K) {
    counts_gen <- rstan::extract(list_of_stanfits[[i]], pars = "counts_gen_geneWise")$counts_gen_geneWise
    holdout <- list_of_holdout[[i]]
    # I should better get samples at fixed steps, but ignoring for now
    samples_to_use <- sample(1:(dim(counts_gen)[1]), N_samples)
    all_genes_holdout_gen[, holdout, ] <- counts_gen[samples_to_use, holdout, ]
  }
  if(any(is.na(counts_gen))) {
    stop("Inconsistent holdout data")
  }

  holdout_ranks_less <- sweep(all_genes_holdout_gen, MARGIN = c(2,3), STATS = base_data$counts, FUN = "<") %>%
    apply(MARGIN = c(2,3), FUN = sum)

  holdout_ranks_equal <- sweep(all_genes_holdout_gen, MARGIN = c(2,3), STATS = base_data$counts, FUN = "==") %>%
    apply(MARGIN = c(2,3), FUN = sum)

  #If there are equal values, sample the rank randomly over the equal values
  #Using a trich with rounded uniform random numbers to get that in a vectorized way
  holdout_ranks = holdout_ranks_less +
    round((holdout_ranks_equal + 1) * array(runif(length(holdout_ranks_equal)), dim(holdout_ranks_equal)) - 0.5)

  dimnames(holdout_ranks) <- dimnames(base_data$counts)
  holdout_ranks
}

extract_holdout_ranks_all <- function(kfold_def, kfold_res) {
  kfold_def$model_defs$id %>% map(function(id) {
    indices = kfold_def$model_defs_kfold$id == id
    fits_for_model <- kfold_res$fits[indices]
    holdout_for_model <- lapply(kfold_def$data_list[indices], '[[', "holdout")
    extract_holdout_ranks(fits_for_model, holdout_for_model, kfold_def$base_data) %>%
      as.tibble() %>%
      rownames_to_column("sample") %>%
      gather("gene","rank", - sample) %>%
      mutate(model =  kfold_def$model_defs$model_name[kfold_def$model_defs$id == id])
  }) %>% do.call(rbind, args = .)
}

plot_holdout_ranks <- function(ranks, binwidth = 1, facet = ~ model) {
  if(100 %% binwidth != 0) {
    stop("binwidth has to divide 100")
  }

  n_ranks <- aggregate(update.formula(facet, rank ~ .) , ranks, length)
  if(length(unique(n_ranks$rank)) != 1) {
    stop("Unequal number of observations per group")
  }

  n_ranks <- n_ranks$rank[1]

  CI = qbinom(c(0.005,0.5,0.995), size=n_ranks,prob  =  binwidth / 100)
  lower = CI[1]
  mean = CI[2]
  upper = CI[3]


  ranks %>% ggplot(aes(x = rank)) +
    geom_segment(aes(x=0,y=mean,xend=100,yend=mean),colour="grey25") +
    geom_polygon(data=data.frame(x=c(-10,0,-10,110,100,110,-10),y=c(lower,mean,upper,upper,mean,lower,lower)),aes(x=x,y=y),fill="grey45",color="grey25",alpha=0.5) +
    geom_histogram(breaks =  seq(1, 101, by = binwidth), closed = "left" ,fill="#A25050",colour="black") +
    facet_wrap(facet, scales = "free_y") +
    ggtitle("Posterior ranks of heldout observations")

}

#Functions extract_log_lik and kfold taken from
#https://github.com/stan-dev/stancon_talks/blob/master/2017/Contributed-Talks/07_nicenboim/kfold.Rmd

extract_log_lik_K <- function(list_of_stanfits, list_of_holdout, ...){
  K <- length(list_of_stanfits)
  list_of_log_liks <- plyr::llply(1:K, function(k){
    extract_log_lik(list_of_stanfits[[k]],...)
  })
  # `log_lik_heldout` will include the loglike of all the held out data of all the folds.
  # We define `log_lik_heldout` as a (samples x N_obs) matrix
  # (similar to each log_lik matrix)
  log_lik_heldout <- list_of_log_liks[[1]] * NA
  for(k in 1:K){
    log_lik <- list_of_log_liks[[k]]
    samples <- dim(log_lik)[1]
    N_obs <- dim(log_lik)[2]
    # This is a matrix with the same size as log_lik_heldout
    # with 1 if the data was held out in the fold k
    heldout <- matrix(rep(list_of_holdout[[k]], each = samples), nrow = samples)
    # Sanity check that the previous log_lik is not being overwritten:
    if(any(!is.na(log_lik_heldout[heldout==1]))){
      warning("Heldout log_lik has been overwritten!!!!")
    }
    # We save here the log_lik of the fold k in the matrix:
    log_lik_heldout[heldout==1] <- log_lik[heldout==1]
  }
  return(log_lik_heldout)
}

kfold <- function(log_lik_heldout)  {
  library(matrixStats)
  logColMeansExp <- function(x) {
    # should be more stable than log(colMeans(exp(x)))
    S <- nrow(x)
    colLogSumExps(x) - log(S)
  }
  # See equation (20) of @VehtariEtAl2016
  pointwise <-  matrix(logColMeansExp(log_lik_heldout), ncol= 1)
  colnames(pointwise) <- "elpd"
  # See equation (21) of @VehtariEtAl2016
  elpd_kfold <- sum(pointwise)
  se_elpd_kfold <-  sqrt(ncol(log_lik_heldout) * var(pointwise))

  estimates <- matrix(NA_real_, nrow = 1, ncol = 2)
  rownames(estimates) <- "elpd_loo"
  colnames(estimates) <- c("Estimate","SE")
  estimates["elpd_loo", "Estimate"] <- elpd_kfold
  estimates["elpd_loo", "SE"] <- se_elpd_kfold

  out <- list(
    pointwise = pointwise,
    estimates = estimates
  )
  structure(out, class = "loo")
}

