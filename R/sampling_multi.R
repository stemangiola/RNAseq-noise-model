#' @param models either a single stanmodel object which is used for
#'   all data or a list of models the same length as data
#' @param data a list of datasets to pass to each fit
#' @param map_fun a function to process each single chain fit. It has to take 3 params:
#'   a fit (single chain), data index and chain index. By default it is a noop.
#'   This function is run in parallel.
#'   Note that this function must be runnable in new RStudio sessions. You may
#'   use the `R_session_init` or `map_fun_dependencies` parameters to ensure it is loaded.
#' @param combine_fun a function called to combine the results of applying `map_fun`
#'   to individual chain fits. It should accept a list of values of the `map_fun` calls.
#'   Default value is [sflist2stanfit()] which creates a combined multichain fit.
#'   As of now, this function is run serially after all fits have bin computed.
#' @param ids_to_compute if given the computation is done only for certain items in
#'   the data list (useful for restarting long computations after a failure)
#' @param init_per_item optional list of inits of the same length as data
#' @param control_per_item optional list of control arguments of the same length as data
#' @param map_fun_dependencies a list of package names that need te be loaded for
#'   `map_fun` to run. IMPORTANT: when developing packages, you need to install
#'   the latest version (the packages are loaded from the default library)
#' @param ... Other params to be passed the [sampling()].
#'
#' @return A list of length `length(data)` containing the result of applying
#'   `combine_fun` to each set of `chains` result of applying `map_fun`
#'   to each fit.
sampling_multi <- function(models, data, map_fun = sampling_multi_noop, combine_fun = sflist2stanfit, chains = 4, cores = parallel::detectCores(),
                           init = NULL, control = NULL, init_per_item = NULL, control_per_item = NULL,
                           map_fun_dependencies = c(),
                           R_session_init_expr = NULL,
                           ids_to_compute = 1:length(data),  ...) {

  #TODO: argument validation

  cl <- parallel::makeCluster(cores, useXDR = FALSE)
  on.exit(parallel::stopCluster(cl))

  #Prepare argument lists for all calls
  .dotlists_per_item <- list()
  for(data_id in 1:length(data)) {
    if(is.list(models)) {
      model = models[[data_id]]
    } else {
      model = models
    }

    for(chain_id in 1:chains) {

      .dotlist <- c(list(model, chains = 1L, cores = 0L), list(...))

      if(is.list(init_per_item)) {
        .dotlist$init <- init_per_item[[data_id]]
        #TODO check that init is not set
      } else {
        .dotlist$init <- init
      }

      #Handle per-chain inits
      if(is.list(.dotlist$init)) .dotlist$init <- .dotlist$init[chain_id]

      if(is.list(control_per_item)) {
        .dotlist$control <- control_per_item[[data_id]]
        #TODO check that control is not set
      } else {
        .dotlist$control <- control
      }

      .dotlist$chain_id <- chain_id

      .dotlist$data <- data[[data_id]]


      .dotlists_per_item[[(data_id - 1) * chains + chain_id]] <- .dotlist
    }
  }


  fit_fun <- function(i) {
    out <- do.call(rstan::sampling, args = .dotlists_per_item[[i]])

    #TODO should catch error from map_fun
    map_fun(out, data_id, chain_id)
  }



  dependencies <- c("rstan", "Rcpp", map_fun_dependencies)
  .paths <- unique(c(.libPaths(), sapply(dependencies, FUN = function(d) {
    dirname(system.file(package = d))
  })))
  .paths <- .paths[.paths != ""]
  parallel::clusterExport(cl, varlist = c(".paths","dependencies"), envir = environment())
  parallel::clusterEvalQ(cl, expr = .libPaths(.paths))
  parallel::clusterEvalQ(cl, expr =
                           for(dep in dependencies) {
                             suppressPackageStartupMessages(require(dep, quietly = TRUE, character.only = TRUE))
                           }
                         )

  parallel::clusterExport(cl, varlist = ".dotlists_per_item", envir = environment())

  parallel::clusterExport(cl, varlist =
                            c("map_fun", "R_session_init_expr"),
                          envir = environment())

  parallel::clusterEvalQ(cl, expr = R_session_init_expr)

  ids <- rep(ids_to_compute, each = chains)
  items <- ((rep(ids_to_compute, each = chains) - 1) * chains  + rep(1:chains, times = length(ids_to_compute)))

  results_flat <-  parallel::parSapplyLB(cl, X = items, FUN = fit_fun)

  results <- list()
  for(i in 1:length(data)) {
    results[[i]] <- combine_fun(results_flat[ids == i])
  }
  results
}

sampling_multi_noop <- function(fit, data_id, chain_id) {
  fit
}

#' @return Returns a function to be used as `map_fun` in [sampling_multi()]
#'   that stores all fits in RDS files
sampling_multi_store_file_generator <- function(base_dir, base_name) {
  function(fit, data_id, chain_id) {
    filename = paste0(base_dir,"/",base_name, "_", data_id, "_", chain_id,".rds")
    saveRDS(fit, filename)
    filename
  }
}
