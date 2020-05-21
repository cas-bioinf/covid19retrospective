missing_arg <- function() quote(expr=)

n_missing_args <- function(n) {
  res <- list()
  for(i in 1:n) {
    res[[i]] <- missing_arg()
  }
  res

}

order_within_samples <- function(x, samples, n_samples = 99) {
  thinned_samples <- sample(samples, n_samples)
  # For discrete parameters, I need to uniformly samples across the positions that are equal
  sum(thinned_samples < x) + rdunif(1, a = 0, b = sum(thinned_samples == x)) + 1
}

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
    mean = mean(param_samples),
    IQR = IQR(param_samples),
    quantile = ecdf(param_samples)(true_value),
    order_within = order_within_samples(true_value, param_samples)
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
    diags = rstan::summary(fit)$summary[eval_result$param_name , , drop = FALSE]

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


# Code from https://betanalpha.github.io/assets/case_studies/rstan_workflow.html

#' Analyse divergences in a stanfit object
#'
#' @param sf stanfit object.
#' @param nupars either the string 'all', or an integer reflecting how many pars
#' (from first to nupars) to use.
#'
#' @return A list of four matrices. $locationsort and $sdsort contian the bivariate interactions of
#' unconstrained parameters, sorted by either the relative location of any divergences, or the relative standard deviation.
#' $locationmeans and $sdmeans collapse across the bivariate interactions to return the means for each parameter.
#' @export
#'
#' @examples
stan_checkdivergences <- function(sf,nupars = 'all'){

	samplerps <- get_sampler_params(sf)
	skeleton <- get_inits(sf)[[1L]]
	if('all' %in% nupars) nupars <- get_num_upars(sf)
	e <- rstan::extract(sf)
	ea <- as.array(sf) #[,,1:nupars]

	ucsnames <- dimnames(ea)$parameters[1:nupars]

	#unconstrain pars into iter * chain * par array
	ucs <- array(NA,dim=c(dim(ea)[c(1,2)],nupars))
	for(i in 1:dim(ea)[1]){
		for(j in 1:dim(ea)[2]){
			ucs[i,j,] <-  unconstrain_pars(sf,  relist(ea[i,j,],skeleton))[1:nupars]
		}}

	#forget chains
	ucs2 <- matrix(ucs,nrow=prod(dim(ucs)[c(1,2)]), ncol=dim(ucs)[3])

	#get divergences
	dvg <- unlist(lapply(samplerps, function(x) x[(nrow(x)-dim(ea)[1]+1):nrow(x),'divergent__']))

	#seperate samples into diverging and not
	ucsbad <- ucs2[dvg==1,]
	ucsgood <- ucs2[dvg==0,]

	#calculate bivariate summary stats
	ints <- matrix(NA,nrow=6, ncol= (ncol(ucs2)^2-ncol(ucs2))/2)
	count <- 0
	cnames <- c()
	for(i in 1:(ncol(ucs2)-1)){
		for(j in (i+1):ncol(ucs2)){
			count = count +1
			ibadgc <- (ucsbad[,i]-mean(ucsgood[,i])) * (ucsbad[,j]-mean(ucsgood[,j])) #divergent interactions centered around good interactions
			igoodgc <- (ucsgood[,i]-mean(ucsgood[,i])) * (ucsgood[,j]-mean(ucsgood[,j]))#centered regular interactions centered around good interactions
			ibadbc <- (ucsbad[,i]-mean(ucsbad[,i])) * (ucsbad[,j]-mean(ucsbad[,j])) #divergent interactions centered around bad interactions
			ints[1,count] <- mean(ibadgc)
			ints[2,count] <- mean(igoodgc)

			ints[3,count] <- sd(ibadbc)
			ints[4,count] <- sd(igoodgc)
			cnames <- rbind(cnames,c(ucsnames[i], ucsnames[j]))
		}
	}

	colnames(cnames) <- c('par1','par2')
	rownames(ints) <- c('mean_dvg','mean_good','sd_dvg_badc','sd_good_goodc','dvg_rel_location','dvg_rel_sd')

	ints[5,] <-  (ints[1,] - ints[2,])/ ints[4,] #location of mean of divergences relative to mean / sd of good samples
	ints[6,] <- ints[3,] / ints[4,] #ratio of sd of divergences to sd of good samples

	#get ordering of par interactions by stat
	badmean <- order(abs(ints[5,]),decreasing = TRUE)
	badsd <- order(ints[6,])

	ints <- cbind(cnames,t(ints))
	sdsort=ints[badsd,c('par1','par2','dvg_rel_location','dvg_rel_sd')]
	locationsort = ints[badmean,c('par1','par2','dvg_rel_location','dvg_rel_sd')]

	# sdscores <- cbind(ucsnames,unlist(lapply(ucsnames,function(x) sum(c(which(sdsort[,'par1'] %in% x),which(sdsort[,'par2'] %in% x))))))
	# locationscores <- cbind(ucsnames,unlist(lapply(ucsnames,function(x) sum(c(which(sdsort[,'par1'] %in% x),which(locationsort[,'par2'] %in% x))))))

	sdmeans <- cbind(ucsnames,unlist(lapply(ucsnames,function(x) {
		mean(as.numeric(c(ints[which(ints[,'par1'] %in% x),'dvg_rel_sd'],
											ints[which(ints[,'par2'] %in% x),'dvg_rel_sd'])))
	})))

	locationmeans <- cbind(ucsnames,unlist(lapply(ucsnames,function(x) {
		mean(as.numeric(c(ints[which(ints[,'par1'] %in% x),'dvg_rel_location'],
											ints[which(ints[,'par2'] %in% x),'dvg_rel_location'])))
	})))


	sdmeans <- sdmeans[order(sdmeans[,2]),]
	locationmeans <- locationmeans[order(abs(as.numeric(locationmeans[,2])),decreasing = TRUE),]

	rownames(sdmeans) <- sdmeans[,1]
	rownames(locationmeans) <- locationmeans[,1]
	sdmeans <- sdmeans[,2,drop=FALSE]
	locationmeans <- locationmeans[,2,drop=FALSE]
	sdmeans[,1] <- as.numeric(sdmeans[,1])
	locationmeans[,1] <- as.numeric(locationmeans[,1])

	#output
	out<-list(sdmeans=sdmeans,locationsort=locationsort,sdsort=sdsort,locationmeans=locationmeans)

	return(out)
}

# Check transitions that ended with a divergence
check_div <- function(fit) {
	sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
	divergent <- do.call(rbind, sampler_params)[,'divergent__']
	n = sum(divergent)
	N = length(divergent)

	print(sprintf('%s of %s iterations ended with a divergence (%s%%)',
								n, N, 100 * n / N))
	if (n > 0)
		print('  Try running with larger adapt_delta to remove the divergences')
}

# Check transitions that ended prematurely due to maximum tree depth limit
check_treedepth <- function(fit, max_depth = 10) {
	sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
	treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
	n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
	N = length(treedepths)

	print(sprintf('%s of %s iterations saturated the maximum tree depth of %s (%s%%)',
								n, N, max_depth, 100 * n / N))
	if (n > 0)
		print('  Run again with max_depth set to a larger value to avoid saturation')
}

# Checks the energy Bayesian fraction of missing information (E-BFMI)
check_energy <- function(fit) {
	sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
	no_warning <- TRUE
	for (n in 1:length(sampler_params)) {
		energies = sampler_params[n][[1]][,'energy__']
		numer = sum(diff(energies)**2) / length(energies)
		denom = var(energies)
		if (numer / denom < 0.2) {
			print(sprintf('Chain %s: E-BFMI = %s', n, numer / denom))
			no_warning <- FALSE
		}
	}
	if (no_warning)
		print('E-BFMI indicated no pathological behavior')
	else
		print('  E-BFMI below 0.2 indicates you may need to reparameterize your model')
}

# Checks the effective sample size per iteration
check_n_eff <- function(fit) {
  fit_summary <- rstan::summary(fit, probs = c(0.5))$summary
	N <- dim(fit_summary)[[1]]

	iter <- dim(rstan::extract(fit)[[1]])[[1]]

	no_warning <- TRUE
	for (n in 1:N) {
		ratio <- fit_summary[,5][n] / iter
		if (ratio < 0.001) {
			print(sprintf('n_eff / iter for parameter %s is %s!',
										rownames(fit_summary)[n], ratio))
			no_warning <- FALSE
		}
	}
	if (no_warning)
		print('n_eff / iter looks reasonable for all parameters')
	else
		print('  n_eff / iter below 0.001 indicates that the effective sample size has likely been overestimated')
}

# Checks the potential scale reduction factors
check_rhat <- function(fit) {
	fit_summary <- rstan::summary(fit, probs = c(0.5))$summary
	N <- dim(fit_summary)[[1]]

	no_warning <- TRUE
	for (n in 1:N) {
		rhat <- fit_summary[,6][n]
		if (rhat > 1.1 || is.infinite(rhat) || is.nan(rhat)) {
			print(sprintf('Rhat for parameter %s is %s!',
										rownames(fit_summary)[n], rhat))
			no_warning <- FALSE
		}
	}
	if (no_warning)
		print('Rhat looks reasonable for all parameters')
	else
		print('  Rhat above 1.1 indicates that the chains very likely have not mixed')
}

check_all_diagnostics <- function(fit) {
	check_n_eff(fit)
	check_rhat(fit)
	check_div(fit)
	check_treedepth(fit)
	check_energy(fit)
}

# Returns parameter arrays separated into divergent and non-divergent transitions
partition_div <- function(fit) {
	nom_params <- rstan::extract(fit, permuted=FALSE)
	n_chains <- dim(nom_params)[2]
	params <- as.data.frame(do.call(rbind, lapply(1:n_chains, function(n) nom_params[,n,])))

	sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
	divergent <- do.call(rbind, sampler_params)[,'divergent__']
	params$divergent <- divergent

	div_params <- params[params$divergent == 1,]
	nondiv_params <- params[params$divergent == 0,]

	return(list(div_params, nondiv_params))
}

# Added later



