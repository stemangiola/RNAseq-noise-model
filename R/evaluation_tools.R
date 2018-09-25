missing_arg <- function() quote(expr=)

evaluate_single_param_indices <- function(samples, param_name, indices, true_value) {
  if(is.null(indices)) {
    param_samples = samples[[param_name]];
  }
  else {
    #A magic form to run samples[[param_name]][,indices[1], ... , indices[N]] based on the length of indices
    param_samples = do.call(`[`,append(list(samples[[param_name]],missing_arg()),indices))
  }

  if(length(indices) > 0) {
    indices_str = do.call(paste, append(indices, list(sep = ",")))
    fullName =   paste0(param_name,"[", indices_str, "]")
  } else {
    fullName = param_name
  }

  mad_val = mad(param_samples, center = true_value)
  rmse_val = sqrt(mean((param_samples - true_value) ^ 2))
  return(data.frame(
    param_name = fullName,
    true_value = true_value,
    median = median(param_samples),
    IQR = IQR(param_samples),
    quantile = ecdf(param_samples)(true_value)
    # mad = mad_val,
    # relative_mad = mad_val / true_value,
    # relative_rmse = rmse_val / true_value
  ))
}

evaluate_single_param <- function(samples, param_name, param_values)
{
  result = list();
  dimensions <- dim(samples[[param_name]])[-1] #The first dimension is the number of samples
  num_dimensions <- length(dimensions)
  next_element = 1
  if(num_dimensions == 0) {
    result[[next_element]] = evaluate_single_param_indices(samples, param_name, NULL, param_values)
    next_element = next_element + 1
  } else if (num_dimensions == 1) {
    for(i in 1:dimensions[1]) {
      result[[next_element]] = evaluate_single_param_indices(samples, param_name, list(i), param_values[i])
      next_element = next_element + 1
    }
  }
  else if(num_dimensions == 2) {
    for(i in 1:dimensions[1]) {
      for(j in 1:dimensions[2]) {
        result[[next_element]] = evaluate_single_param_indices(samples, param_name, list(i,j), param_values[i,j])
        next_element = next_element + 1
      }
    }
  } else {
    stop("3+ dimensional parameters not supported yet");
  }
  return(do.call(rbind.data.frame, result))
}

evaluate_all_params <- function(samples, true_params) {
  result = list();
  next_element = 1;
  for(param_name in names(true_params)) {
    if(!param_name %in% names(samples)) {
      next;
    }
    param_values = get(param_name, true_params);
    result[[next_element]] = evaluate_single_param(samples, param_name, param_values)
    next_element = next_element + 1
  }
  return(do.call(rbind.data.frame, result));
}

evaluation_summary <- function(fit, true_params, printParamsResults = TRUE) {
  samples = rstan::extract(fit)
  eval_result = evaluate_all_params(samples, true_params);

  if(printParamsResults) {
    #Add convergence diagnostics
    diags = rstan::summary(fit)$summary[eval_result$param_name ,]

    eval_result = eval_result %>% mutate(n_eff = diags[,"n_eff"], Rhat = diags[,"Rhat"])

    print(eval_result);
  }
  quantiles = eval_result$quantile;
  within25 = mean(quantiles >= 0.375 & quantiles <= 0.625);
  within50 = mean(quantiles >= 0.25 & quantiles <= 0.75);
  within95 = mean(quantiles >= 0.025 & quantiles <= 0.975);
  cat("\nWithin 25% interval:", within25,"\nWithin 50% interval:", within50, "\nWithin 95% interval:",within95,"\n")
}

averageSamplingTime <- function(fits)
{
  timeList = lapply(fits, get_elapsed_time)
  allTimes = Reduce(rbind,timeList, array(0,c(0,2)))
  warmupTimes = allTimes[,"warmup"]
  sampleTimes = allTimes[,"sample"]
  return(list(total = mean(warmupTimes + sampleTimes), sample = mean(sampleTimes)))
}


launch_shinystan_nonblocking <- function(fit) {
  library(future)
  library(shinystan)
  plan(multisession)
  future(launch_shinystan(fit))
}
