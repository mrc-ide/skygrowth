#' derive timeseries of coalescent and ltt along appropriate time axis 
#' @param haxis vector retrospective time axis increasing order 
.tre2df <- function( tre, haxis){
	res <- length(haxis)
	n <- length( tre$tip )
	D <- node.depth.edgelength( tre )
	rh <- max( D[1:n] )
	sts <- D[1:n]
	
	shs <- max(sts) - sts
	
	inhs <- max(sts) - D[ (n+1):(n + tre$Nnode) ]
	
	ltt.h <- function(h) sum( shs < h ) - sum( inhs < h )
	nco <- sapply( 2:res, function(i) sum( ( inhs >  haxis[i-1] ) & ( inhs <= haxis[i] ) ) )
	ltt <- sapply( 2:res, function(i) sqrt( ltt.h( haxis[i-1] ) * ltt.h( haxis[i] )) )
		
	# derive terms for coalescent likelihood 
	#alpha = {A choose 2} ; gamma = 1/ Ne
	#alpha gamma exp( -alpha gamma dh )
	# log(alpha) + log(gamma)^nc - alpha gamma dh 
	events <- cbind( c( haxis[-1], inhs, shs ), c( rep(0, length(haxis)-1), rep(1, length(inhs)), rep(2, length(shs)) ) )
	events <- events[ order( events[,1]) , ]
	lterms <- matrix(0.,  nrow = length(haxis) -1, ncol = 2 )
	dh <- abs( diff( haxis)[1] )
	.h2i <- function(hh) min( 1 + floor( hh / dh ), length( haxis) -1 )
	lasth <- 0
	difh <- NA
	for (k in 1:nrow(events)){
		et <- events[k,2] # event type
		hh <- events[k,1] 
		i <- .h2i(  hh ) 
		if (i > 0 & i <= length(haxis) )
		{
			difh <- hh - lasth
			if (et == 1 ){ #co
				A <- max(1, ltt.h( hh ) )
				alpha <- ( A * (A - 1)/2 )
				lterms[i,1] <- lterms[i,1] + log(alpha) # + gamma
				lterms[i,2] <- lterms[i,2] + difh * alpha # * gamma 
				lasth <- hh 
			}
			if (et == 2){ #sample
				A <- max(1, ltt.h( hh ) )
				alpha <- ( A * (A - 1)/2 )
				lterms[i,1] <- lterms[i,1] + 0 # no change 
				lterms[i,2] <- lterms[i,2] + difh * alpha # * gamma 
				lasth <- hh 
			}
		}	
		
	}
	# NOTE lterms on forward axis
	lterms <- lterms[ rev( 1:nrow(lterms)), ]
	
	# NOTE not counting most recent sample 
	data.frame( heights = haxis[-1],  nco = nco , ltt = ltt, lterms = lterms)
}


.process.tau_logprior <- function(tau_logprior, tau0){
	if (is.character(tau_logprior)){
		if ( tau_logprior == 'exponential' ){
			tau_logprior <- function(x) dexp(x, 1/tau0, log=T)
		} else if( tau_logprior == 'gamma' ){
			#tau_logprior <- function(x) dgamma( 1/x, shape=sqrt(1/tau0), scale=sqrt(1/tau0), log=T)
			# inflate variance :
			tau_logprior <- function(x) dgamma( x, shape=(tau0)^(1/10), scale=(tau0)^(9/10), log=T)
		} else if( tau_logprior == 'oneOnX'){
			tau_logprior <- function(x) -log(x)
		}
		else{
			stop('tau_logprior must be exponential, gamma, oneOnX, NULL, or user defined log density' )
		}
	}
	tau_logprior 
}

.meanlogexp <- function(x)
{
	m <- max( x)
	log( mean( exp( x - m ) ) ) + m
}

#' Maximum a posteriori estimate of effective size through time with a non-parametric growth model
#'
#' @param tre A dated phylogeny in ape::phylo format (see documentation for ape) or a list of phylogenies or multi.phylo object
#' @param tau0 Initial guess of the precision parameter
#' @param tau_logprior Prior for precision parameter (character string (gamma or exponential) or function)
#' @param res Number of time intervals (integer)
#' @param quiet Provide verbose output? 
#' @param maxiter Maximum number of iterations
#' @param abstol Criterion for convergence of likelihood
#' @param control List of options passed to optim
#' @param maxHeight The earliest time (furthest time in to the past) to estimate Ne. If not provided, will use the median of root heights for input trees 
#' @param ... Not implemented
#' @return A fitted model including effective size through time
#' @export
#' @examples
#' \doNotRun{
#' require(skygrowth)
#' require(ape)
#' load( system.file( package='skygrowth', 'NY_flu.rda' , mustWork=TRUE) ) 
#' # NOTE branch lengths in weeks  / 13 years in all
#' fit <- skygrowth.map( NY_flu 
#'  , res = 24*13  # Ne changes every 2 weeks
#'  , tau0 = .1    # Smoothing parameter. If prior is not specified, 
#'                 # this will also set the scale of the prior
#' )
#' plot( fit ) + scale_y_log10()
#' }
skygrowth.map <- function(tre
  , tau0 = 10
  , tau_logprior = 'exponential'
  , res = 50 
  , quiet = FALSE
  , maxiter = 20
  , abstol = 1e-2
  , control = NULL
  , maxHeight = NULL 
  , ... #not implemented
){
	if (class(tre)=='phylo'){
		tres <- list( tre )
	}else if( class(tre)=='multi.phylo' | class(tre)=='list'){
		tres <- tre
	} else{
		stop('*tre* must be a ape::phylo or multi.phylo or list of ape::phylo')
	}
	
	n <- Ntip(tres[[1]])
	if (is.null(maxHeight)){
		maxHeight <- median(sapply( tres, function(tr){ #NOTE median not max
			D <- node.depth.edgelength( tr )
			max( D[1:n] )
		})) 
	}
	haxis <- seq( 0, maxHeight, length.out = res+1 )
	dh <- abs( haxis[2] - haxis[1] )
	## ne0 
	ne0 <- median( sapply( tres, function(tre ){
		coint <- coalescent.intervals( tre )
		with( coint , {
			abs(interval.length) * ( lineages * (lineages-1) / 2) 
		}) -> ne
		ne[ ne == 0 ] <- NA
		median( ne, na.rm=T)
	}))
	
	tredats <- lapply( tres, function(tre) .tre2df( tre, haxis  ))
	lterms_list  <- lapply( tredats, function(tredat) cbind( tredat$lterms.1, tredat$lterms.2) )
	
	tau_logprior <- .process.tau_logprior( tau_logprior , tau0)
	
	#  update to use lterms  / maybe call on c code
	.of1.2 <- function(  logne , xtau = tau0)
	{
		fwdne <- exp(logne)
		grs <- ( diff ( fwdne ) / ( fwdne[-res] ) / dh )
		ll <- 0
		lls <- sapply( 1:length(tredats), function(k)
		{
			ll <-  sum(  lterms_list[[k]][,1] + log( 1/fwdne ) * rev( tredats[[k]]$nco ) )
			ll <- ll - sum(  lterms_list[[k]][,2] / fwdne )
			ll <- ll + sum( dnorm( diff( grs ), 0, sqrt( dh / xtau ) , log = TRUE) )  + tau_logprior( xtau )
		})
		ll <- .meanlogexp( lls )
		if (is.na(ll) | is.infinite(ll))  ll <- -1e12
		ll
	}
	
	.of2.2 <- function( logtau, xne =ne)
	{
		.of1.2( log( xne ) , xtau = exp( logtau ) )
	}
	
	ne <- rlnorm( res , log( ne0 ), .2 ) # add some jitter
	
	optim( par = log(ne), fn = .of1.2
	  , method = 'BFGS'
	  , control = list( trace = 1, fnscale = -1, parscale = rep(median( abs(log(ne))), length(ne)  ) )
	) -> fit
	
	trace <- matrix( NA, nrow = maxiter, ncol = 2 + res )
	colnames(trace ) <- c('loglik', 'tau', paste('ne', 1:res, sep=''))
	tau <- tau0
	lastll <- -Inf
	{
		for (iter in 1:maxiter){
			ne <- exp(fit$par )
			optim( par =log(tau), fn = .of2.2
				  , method = 'BFGS'
			  , control = list( trace = 1, fnscale  = -1)
			  , xne = ne
			) -> fit_tau
			tau <- exp( fit_tau$par )
			optim( par = log(ne), fn = .of1.2
			  , method = 'BFGS'
			  , control = list( trace = 1, fnscale  = -1, parscale = abs(rep(median( log(ne)), length(ne)  ) ) )
			  , xtau = tau
			) -> fit
			trace[ iter,1] <-fit$value
			trace[ iter,2] <- tau
			trace[ iter,3:ncol(trace)] <- ( exp(fit$par) )
			if ( fit$value - lastll < abstol) break;
			lastll <- fit$value
			
			cat( 'iter\n')
			print( iter)
			print( paste( tau, fit$value) )
		}
	}
	trace <- trace[1:iter,] 
	cat( 'Computing hessian...\n')
	optim( par = log(ne), fn = .of1.2
	  , method = 'BFGS'
	  , control = list( trace = 1, fnscale  = -1, parscale = abs(rep(median( log(ne)), length(ne)  ) ) )
	  , xtau = tau
	  , hessian=TRUE 
	) -> f
	fi <- tryCatch( solve( -f$hessian), error = function(e) {
		warning('Hessian could not be computed. Will not compute CIs.')
		NA
	})
	fsigma <- if (!any(is.na(fi))) {
		sqrt( diag( fi ))
	} else{
		NA
	}
	
	ne <- (exp(fit$par) ) #note forward time now
	nelb <- exp( (fit$par) - fsigma*1.96 )
	neub <- exp( (fit$par) + fsigma*1.96 )
	ne_ci <- cbind( nelb, ne, neub )
	
	growthrate <-  c( diff ( exp( fit$par) ) / ( exp(fit$par)[-res] ) / dh , NA) 
	
	rv <- list( 
		ne =  ne
	  , ne_ci = ne_ci  
	  , growthrate =  growthrate
	  , tau = tau
	  , time = haxis[-1]
	  , tredat = tredats
	  , gamma = NA
	  , control = control
	  , tre = tre	
	  , trace = trace 
	  , loglik = fit$value
	  , sigma = (fsigma )
	)
	
	class(rv) <- 'skygrowth.map'
	rv
}

#' Maximum a posteriori estimate of effective size through time using covariate data
#'
#' @param tre A dated phylogeny in ape::phylo format (see documentation for ape) or a list of phylogenies or multi.phylo object
#' @param formula An R formula with empty left-hand-side; the right-hand-side specifies relationship of covariates with growth rate of Ne
#' @param data A data.frame, must include 'time' column
#' @param maxSampleTime The scalar time that the most recent sample was collected
#' @param tau0 Initial guess of the precision parameter
#' @param tau_logprior Prior for precision parameter (character string (gamma or exponential) or function)
#' @param beta_logpriors Optional list of functions providing log density for coefficients (must correspond to data)
#' @param res Number of time points (integer)
#' @param quiet Provide verbose output? 
#' @param maxiter Maximum number of iterations
#' @param abstol Criterion for convergence of likelihood
#' @param control List of options passed to optim
#' @param maxHeight The earliest time (furthest time in to the past) to estimate Ne
#' @param ... Not implemented
#' @return A fitted model including effective size through time
#' @export
skygrowth.map.covar =skygrowth.map.covars <- function(tre
  , formula # should not have left hand side 
  , data # data.frame must include 'time' 
  , maxSampleTime # required to relate to covars 
  , tau0 = 10
  , tau_logprior = 'exponential'
  , beta_logpriors = list()
  , res = 50 
  , quiet = F
  , maxiter = 20
  , abstol = 1e-2
  , control = NULL
  , maxHeight = NULL
  , ... #not implemented
){
	stopifnot( class(formula)=='formula' )
	if (!('time' %in% colnames(data))) stop('covariate data must include *time* of observation' )
	# fit w/o covars first 
	mapfit <- skygrowth.map(tre
	  , tau0 = tau0
	  , tau_logprior = tau_logprior
	  , res = res
	  , quiet = quiet
	  , maxiter = min( 20, maxiter )
	  , abstol = 1e-4
	  , control = control
	  , maxHeight = maxHeight 
	)
	tredats <- mapfit$tredat
	lterms_list  <- lapply( tredats, function(tredat) cbind( tredat$lterms.1, tredat$lterms.2) )
	dh <- abs(diff(mapfit$time)[1] )
	
	
	#ne <- tredat$ne0
	ne <- ( mapfit$ne )
	tau0 <- mapfit$tau 
		
	X0 <- as.data.frame( model.matrix(  formula , data ) )
	betanames <- colnames( X0 )[-1]
	X0 <- cbind( time = data$time 
	 , X0 )
	
	mapfit$time <- maxSampleTime  + mapfit$time - max(mapfit$time)
	#rtredat <- tredat[ rev(1:nrow(tredat)), ]
	#covar.df <- data.frame( time = mapfit$time, gr0 = c( mapfit$growthrate, NA) )
	covar.df <- data.frame( time =mapfit$time, gr0 =  mapfit$growthrate )
	for ( bn in betanames ){
		itime <- setdiff( order( X0$time ), which(is.na( X0[[bn]] )) )
		#covar.df[[bn]] <- approx( X0$time[itime] , X0[[bn]][itime], xout = mapfit$time, rule = 2)$y
		# allow NA
		covar.df[[bn]] <- approx( X0$time[itime] , X0[[bn]][itime], xout = mapfit$time)$y
	}
	covar.df <- covar.df[ order( covar.df$time) , ] 
	lmfit <- lm ( formula(
		   paste0( paste( 'gr0', paste( betanames , collapse = '+'), sep='~' ) , '-1' )
		 )
		 , data = covar.df 
		)
	beta0 <- coef(lmfit)
	
	if (length(beta_logpriors)==0){
		warning('No prior distributions provided for regression coefficients. Will use normal distribution with unit variance. Covariates should be centred and rescaled before being passed to this function.')
		beta_logpriors <- setNames( lapply( betanames, function(x) function(xx) dnorm( xx, 0, 1, log=T) ) , betanames)
	}
	beta.logprior <- function(beta){
		sum( sapply( names(beta), function(x)  beta_logpriors[[x]]( beta[x] ) ) )
	}
	
	tau_logprior <- .process.tau_logprior( tau_logprior, tau0 )
	
	#  update to use lterms  / maybe call on c code
	.prior.gr0 <- function ( xtau, xbeta, zxb, xne){
		# fix endpoints 
		dzxb <- diff( zxb )[-(res-1)]
		dzxb[is.na(dzxb)] <- 0
		grs <- ( diff ( xne ) / ( xne[-res] ) / dh )
		ll <- sum( dnorm( diff( grs ), dzxb, sqrt( dh / xtau ) , log = TRUE) )  + tau_logprior( xtau )
		ll
	}
	
	.of3.1 <- function(  logne ,zxb,  xbeta,  xtau = tau0)
	{
		fwdne = xne <- exp(logne)
		
		logprior <- .prior.gr0 ( xtau, xbeta, zxb, xne )
		ll <- 0
		lls <- sapply( 1:length(tredats), function(k)
		{
			ll <-  sum(  lterms_list[[k]][,1] + log( 1/fwdne ) * rev( tredats[[k]]$nco ) )
			ll <- ll - sum(  lterms_list[[k]][,2] / fwdne )
			ll <- ll + logprior
		})
		ll <- .meanlogexp( lls )
		
		ll
	}
	
	.of3.logtau <- function( logtau, logne ,zxb,  xbeta)
	{
		.of3.1( logne , zxb, xbeta, xtau = exp( logtau ) )
	}
	
	.of3.beta <- function( xbeta,  xtau, logne)
	{
		zxb <- beta2zxb ( xbeta )
		.of3.1(logne = logne , zxb = zxb, xbeta = xbeta, xtau = xtau )
	}

	##
	beta2zxb <- function( beta ){
		if ( length( betanames ) > 1 ){
			zedCrossBeta <- as.vector( as.matrix(covar.df[, betanames]) %*% beta )
		} else{
			zedCrossBeta <- covar.df[, betanames] * beta 
		}
		zedCrossBeta
	}
	
	
	## 
	beta <- beta0
	zxb <- beta2zxb ( beta )
	tau <- tau0
	optim( par = log(ne), fn = .of3.1
	  , method = 'BFGS'
	  , control = list( trace = 1, fnscale  = -1, parscale = rep(median( abs(log(ne))), length(ne)  ) )
	  , zxb = zxb, xbeta = beta, xtau = tau
	) -> fit
	
	np <- length(betanames )
	trace <- matrix( NA, nrow = maxiter, ncol = 2 + np + res )
	colnames(trace ) <- c('loglik', 'tau', betanames, paste('ne', 1:res, sep=''))
	
	lastll <- -Inf
	{
		for (iter in 1:maxiter){
			logne <- (fit$par )
			
			optim( par =log(tau), fn = .of3.logtau
				  , method = 'BFGS'
			  , control = list( trace = 1, fnscale  = -1)
			  , logne = logne ,zxb = zxb,  xbeta = beta
			) -> fit_tau
			tau <- exp( fit_tau$par )
			
			#.of3.beta <- function( xbeta,  xtau, logne)
			optim( par =beta, fn = .of3.beta
				  , method = 'BFGS'
			  , control = list( trace = 1, fnscale  = -1)
			  , xtau = tau , logne = logne
			) -> fit_beta
			beta <- fit_beta$par
			zxb <- beta2zxb( beta )
			
			optim( par = logne, fn = .of3.1
			  , method = 'BFGS'
			  , control = list( trace = 1, fnscale  = -1, parscale = rep( abs(median( logne)), length(logne)  ) )
			  , xtau = tau, zxb = zxb , xbeta = beta
			) -> fit
			
			trace[ iter,1] <-fit$value
			trace[ iter,2] <- tau
			trace[ iter, 3:(3+np-1)] <- beta
			trace[ iter,(3+np):ncol(trace)] <- ( exp(fit$par) )
			if ( fit$value - lastll < abstol) break;
			lastll <- fit$value
			
			cat( 'iter\n')
			print( iter)
			print( paste( c(tau, beta, fit$value) ))
		}
	}
	trace <- trace[1:iter,] 
	cat( 'Computing hessian...\n')
	optim( par = logne, fn = .of3.1
	  , method = 'BFGS'
	  , control = list( trace = 1, fnscale  = -1, parscale = rep( abs(median( logne)), length(logne)  ) )
	  , xtau = tau, zxb = zxb , xbeta = beta
	  , hessian = TRUE
	) -> f
	fi <- solve( -f$hessian)
	fsigma <- sqrt( diag( fi ))
	
	ne <- (exp(fit$par) ) #note forward time now
	nelb <- exp( (fit$par) - fsigma*1.96 )
	neub <- exp( (fit$par) + fsigma*1.96 )
	ne_ci <- cbind( nelb, ne, neub )
	growthrate <-  c( diff ( exp( fit$par) ) / ( exp(fit$par)[-res] ) / dh , NA)
	
	covar.df$gr <- growthrate 
	covar.df$ne <- ne
	rv <- list( 
		ne =  ne
	  , ne_ci = ne_ci  
	  , growthrate =  growthrate
	  , tau = tau
	  , time = mapfit$time #rev(-tredat$heights )
	  , tredat = tredats
	  , gamma = NA
	  , control = control
	  , tre = tre	
	  , trace = trace 
	  , loglik = fit$value
	  , sigma = (fsigma )
	  , covar.df = covar.df
	  , beta = beta 
	)
	
	class(rv) <- c('skygrowth.map.covar', 'skygrowth.map' )
	rv
}



########################################################################



.default_mcmc_control <- list(ne_samp_size = 100
    , prop_log_tau_sd = NULL
    , ne_steps_per_tau_step = 100 
    , logne_proposal_sd_factor = 1/20
    , thin = 100
    , burnin_percent = 20)

##

#' A gibbs-metropolis algorithm for sampling Ne(t) with a non-parametric growth model
#'
#' @param tre A dated phylogeny in ape::phylo format (see documentation for ape)
#' @param mhsteps Number of mcmc steps
#' @param tau0 Initial guess of the precision parameter
#' @param tau_logprior Prior for precision parameter (character string (gamma or exponential) or function)
#' @param res Number of time points (integer)
#' @param quiet Provide verbose output? 
# @param abstol Criterion for convergence of likelihood
#' @param control List of options passed to optim
#' @param ... Not implemented
#' @return A fitted model including effective size through time
#' @export
skygrowth.mcmc <- function(tre
  , mhsteps = 1e5
  , tau0 = 10
  , tau_logprior = 'exponential'
  , res = 50
  , quiet = F
  , control = NULL
  , ... #not implemented
){
	if( is.null( control)) {
			control <- .default_mcmc_control
	} else{
		x <- .default_mcmc_control
		x[names(control)] <- control
		control <- x
	}
with( control, {
	tredat <- .tre2df( tre, res )
	lterms <- cbind( tredat$lterms.1, tredat$lterms.2) ;
	dh <- abs(diff(tredat$heights)[1] )
	
	tau_logprior <- .process.tau_logprior( tau_logprior , tau0)
	
	mapfit <- skygrowth.map(tre
	  , tau0 = tau0
	  , tau_logprior = tau_logprior
	  , res = res
	  , quiet = quiet
	  , maxiter = 1000
	  , abstol = 1e-4
	  , control = NULL
	)
	if (is.null( prop_log_tau_sd )) prop_log_tau_sd <- .2 + abs( log(mapfit$tau) ) / 5
	#ne <- tredat$ne0
	ne <- ( mapfit$ne )
	tau0 <- mapfit$tau 
#~ ne <- tredat$ne0
	logne_proposal_sd <- median(ne) * logne_proposal_sd_factor
	
	# output to save
	i_burnin <- floor( (burnin_percent / 100) * mhsteps )
	datadim <- floor( (mhsteps-i_burnin) / thin  ) 
	NE <- matrix( NA, nrow = datadim, ncol = res )
	NE[1,] <- ne
	TAU <- rep(NA, datadim )
	TAU[1] <- tau0
	ACCEPT <- rep( 0, datadim )
	GROWTHRATE <- matrix( NA, nrow = datadim, ncol = res-1 )
	
	
	.of1.2 <- function(  xtau , ne)
	{
		fwdlogne <- log( ne )
		
		ll <-  sum(  lterms[,1] + log( 1/ne ) * rev( tredat$nco ) )
		ll <- ll - sum( lterms[,2] / ne )
		
		#grs <- diff ( fwdlogne ) / ( fwdlogne[-res] ) / dh #
		grs <- diff ( ne ) / ( ne[-res] ) / dh #
		ll <- ll + sum( dnorm( diff( grs ), 0, sqrt( dh / xtau ) , log = TRUE) )  + tau_logprior( xtau )
		ll
	}

	
	beta.logprior <- function ( beta ) 0
	
	tau <- tau0
	n_accept <- 0
	for (istep in 1:mhsteps) {
		#gibbs ne
		ne <- (as.vector( 
		  mh_sample_ne1_2( ( ne) , rev( tredat$nco ),  tau, mapfit$sigma, dh, lterms) 
		 ))
		
		#mh move tau
		if (!is.null( tau_logprior )){
		#if (FALSE){
			if (istep > 1 & istep %% ne_steps_per_tau_step == 0){
				proptau <- exp( log(tau) + rnorm( 1, 0, prop_log_tau_sd) )
				lltau <- .of1.2( tau,  ne )
				llproptau <- .of1.2( proptau,  ne )
				if ( runif(1) < exp(llproptau - lltau) ) {
					n_accept <- n_accept + 1
					tau <- proptau
				}
			}
		}
		
		if ( ((istep-i_burnin)>0) & ((istep-i_burnin) %% thin == 0)){
			i <- (istep-i_burnin) / thin
			ACCEPT[i] <- n_accept * ne_steps_per_tau_step/ istep
			TAU[i] <- tau
			NE[i, ]<- ( ne ) 
			#~  		logne <- log(ne ) 
			#~ 			GROWTHRATE[i, ] <-  diff ( (logne) ) / ( (logne[-res]) ) / dh 
			GROWTHRATE[i, ] <-  diff ( (ne) ) / ( (ne[-res]) ) / dh 
		}
		
		if ( istep %% thin == 0){
			if (!quiet){
				print(base::date())
				print(istep)
				print( tau )
				print( n_accept * ne_steps_per_tau_step/ istep )
			}
		}
	}
	
	GROWTHRATE <- cbind( GROWTHRATE, NA ) 
	## compute quantiles
	ne_ci <- t( sapply(1:ncol( NE ), function(i) quantile( NE[, i], prob = c( .025 , .5, .975 ), na.rm=T) ) )
	growthrate_ci <- t( sapply(1:ncol( GROWTHRATE ), function(i) quantile( GROWTHRATE[, i], prob = c( .025 , .5, .975 ), na.rm=T) ) )
	rv <- list( ne = NE
	  , growthrate = GROWTHRATE
	  , ne_ci = ne_ci
	  , growthrate_ci = growthrate_ci 
	  , tau = TAU
	  , time = rev(-tredat$heights ) #+ abs( diff( tredat$heights))[1]
	  , accept = ACCEPT
	  , tredat = tredat
	  , gamma = NA
	  , control = control
	  , tre = tre
	  , tau_logprior = tau_logprior 
	  , mapfit = mapfit 
	)
	
	class(rv) <- 'skygrowth.mcmc'
	rv
})}


#' Continue MCMC chain of previously fitted model
#'
#' @param fit A skygrowth.mcmc fit
#' @param mhsteps Additional MCMC steps to compute
#' @param quiet Verbose output? 
#' @param ... Passed on
#' @return A skygrowth.mcmc fit 
#' @export
continue.skygrowth.mcmc <- function( fit 
 , mhsteps
 , quiet = F
 , ...){
	stopifnot(inherits(fit, "skygrowth.mcmc"))
	skygrowth.mcmc( fit$tre
	  , mhsteps = mhsteps
	  , tau0 = tail(fit$tau,1)
	  , tau_logprior = fit$tau_logprior 
	  , res = fit$res
	  , quiet = quiet
	  , control = fit$control
	  , ... 
	)
}

########################################################################


##


#' A gibbs-metropolis algorithm for sampling Ne(t) with a 1st order moving average model and using covariate data
#'
#' @param tre A dated phylogeny in ape::phylo format (see documentation for ape)
#' @param formula An R formula with empty left-hand-side; the right-hand-side specifies relationship of covariates with growth rate of Ne
#' @param data A data.frame, must include 'time' column
#' @param maxSampleTime The scalar time that the most recent sample was collected
#' @param iter iter
#' @param iter0 iter0
# @param mhsteps Number of mcmc steps
#' @param tau0 Initial guess of the precision parameter
#' @param tau_logprior Prior for precision parameter (character string (gamma or exponential) or function)
#' @param res Number of time points (integer)
#' @param beta_logpriors Optional list of functions providing log density for coefficients (must correspond to data)
#' @param prop_beta_sd Standard deviation of beta proposal kernel
#' @param quiet Provide verbose output? 
# @param abstol Criterion for convergence of likelihood
#' @param control List of options passed to optim
#' @param ... Not implemented
#' @return A fitted model including effective size through time
#' @export
skygrowth.mcmc.covar = skygrowth.mcmc.covars <- function(tre
  , formula # should not have left hand side 
  , data # data.frame must include 'time' 
  , maxSampleTime # required to relate to covars 
  , iter = 1e5 
  , iter0 = 10
  , tau0 = 10
  , tau_logprior = 'exponential' #gamma
  , res = 50
  , beta_logpriors = list()
  , prop_beta_sd = NULL
  , quiet = F
  , control = NULL
  , ... #not implemented
){
	if (!('time' %in% colnames(data))) stop('covariate data must include *time* of observation' )
	# initial fit for guess of ne growth 

	skygrowth.map.covar ( tre, 
	 , formula = formula
	 , data = data
	 , maxSampleTime = maxSampleTime
	 , tau0 = tau0
	 , tau_logprior = tau_logprior
	 , res = res
	 , quiet = quiet
	 , maxiter = iter0
	 , abstol = 1e-4
	 , control = NULL
	) -> mapfit
	
	#ne <- tredat$ne0
	ne = ne0 <- ( mapfit$ne )
	gr0 <- mapfit$growthrate
	tau0 <- mapfit$tau
	
	tredat <- mapfit$tredat
	lterms <- cbind( tredat$lterms.1, tredat$lterms.2) ;
	X0 <- as.data.frame( model.matrix(  formula , data ) )
	betanames <- colnames( X0 )[-1]
	X0 <- cbind( time = data$time 
	 , X0 )
	
	tredat$times <- maxSampleTime - tredat$heights
	tredat$gr0  <- rev( gr0 ) #
	tredat$ne0 <- rev( ne0 )
	for ( bn in betanames ){
		itime <- setdiff( order( X0$time ), which(is.na( X0[[bn]] )) )
		#tredat[[bn]] <- approx( X0$time[itime] , X0[[bn]][itime], xout = tredat$times, rule = 2)$y
		#allow NA
		tredat[[bn]] <- approx( X0$time[itime] , X0[[bn]][itime], xout = tredat$times)$y
	}

	if( is.null( control)) {
			control <- .default_mcmc_control
	} else{
		x <- .default_mcmc_control
		x[names(control)] <- control
		control <- x
	}
	
	beta0 <- mapfit$beta
	
	#proposals
	if (is.null( prop_beta_sd )){
		prop_beta_sd <- setNames( abs(mapfit$beta/10), betanames )
		print( prop_beta_sd ) 
	}
	#proposal for ne
	mapfit_sigma <- pmax( .25*min(abs(log(mapfit$ne))), mapfit$sigma )
	mapfit_sigma <- pmin( .75*max(abs(log(mapfit$ne))), mapfit$sigma )
	if (length(beta_logpriors)==0){
		warning('No prior distributions provided for regression coefficients. Will use normal distribution with unit variance. Covariates should be centred and rescaled before being passed to this function.')
		beta_logpriors <- setNames( lapply( betanames, function(x) function(xx) dnorm( xx, 0, 1, log=T) ) , betanames)
	}
	beta.logprior <- function(beta){
		sum( sapply( names(beta), function(x)  beta_logpriors[[x]]( beta[x] ) ) )
	}
	
with( control, {
	lterms <- cbind( tredat$lterms.1, tredat$lterms.2) ;
	dh <- abs(diff(tredat$heights)[1] )
	
	tau_logprior <- .process.tau_logprior( tau_logprior , tau0)
	
	if (is.null( prop_log_tau_sd )) prop_log_tau_sd <- .2+ abs( log(tau0) ) / 5
	
	tau0 <- mapfit$tau 
	
	# output to save
	i_burnin <- floor( (burnin_percent / 100) * iter )
	datadim <- floor( (iter-i_burnin) / thin  ) 
	NE <- matrix( NA, nrow = datadim, ncol = res )
	NE[1,] <- ne
	BETA <- matrix( NA, nrow = datadim, ncol = length(betanames ))
	TAU <- rep(NA, datadim )
	TAU[1] <- tau0
	ACCEPT <- rep( 0, datadim )
	LOGPO <- rep( NA, datadim )
	GROWTHRATE <- matrix( NA, nrow = datadim, ncol = res-1 )
	
	# 
	.prior.gr0 <- function ( xtau, xbeta, zxb, xne){
		# fix endpoints 
		dzxb <- diff( zxb )
		grs <- ( diff ( xne ) / ( xne[-res] ) / dh )
		dzxb <- dzxb[1:(length( grs)-1)]
		dzxb[is.na(dzxb)] <- 0
		ll <- sum( dnorm( diff( grs ), dzxb, sqrt( dh / xtau ) , log = TRUE) )  + tau_logprior( xtau )
		ll
	}
	.of3.1 <- function( xtau, xbeta, zxb, fwdne)
	{
		ll <-  sum(  lterms[,1] + log( 1/fwdne ) * rev( tredat$nco ) )
		ll <- ll - sum( lterms[,2] / fwdne )
		ll <- ll + .prior.gr0 ( xtau, xbeta, zxb, fwdne )
		ll
	}
	
	beta2zxb <- function( beta ){
		# note tredat in rev order 
		if ( length( betanames ) > 1 ){
			zedCrossBeta <- as.vector( as.matrix(tredat[, betanames]) %*% beta )
		} else{
			zedCrossBeta <- tredat[, betanames] * beta 
		}
		rev( zedCrossBeta )
	}
	
	tau <- tau0
	beta <- beta0
	n_accept <- 0
	ll <- NA
	for (istep in 1:iter) {
		zxb <- beta2zxb( beta )
		#gibbs ne
		ne <- (as.vector( 
		  mh_sample_ne2_2( ne , rev( tredat$nco ) , tau ,  mapfit_sigma , dh, zxb , lterms ) 
		))
		
		#mh move tau
		if (!is.null( tau_logprior )){
			if (istep > 1 & istep %% ne_steps_per_tau_step == 0){
				proptau <- exp( log(tau) + rnorm( 1, 0, prop_log_tau_sd) )
				lltau <- .of3.1(tau, beta, zxb, ne )
				llproptau <- .of3.1(proptau, beta, zxb, ne )
				if ( runif(1) < exp(llproptau - lltau) ) {
					n_accept <- n_accept + 1
					tau <- proptau
					lltau <- llproptau
				}
			}
		}		
		
		# mh move beta's 
		if (istep > 1 & istep %% ne_steps_per_tau_step == 0) {
			ll <- .of3.1(tau, beta,  zxb, ne )
			for (k in 1:length(betanames)){
				x <- rnorm( 1, beta[k], prop_beta_sd[k])
				propbeta <- beta
				propbeta[k] <- unname( x )
				.zxb <- beta2zxb( propbeta )
				.ll <- .of3.1( tau, propbeta , .zxb, ne )
				
				if ( runif(1) < exp(.ll - ll) ) {
					n_accept <- n_accept + 1
					beta <- propbeta
					ll <- .ll
				}
			}
		}
		
		if ( ((istep-i_burnin)>0) & ((istep-i_burnin) %% thin == 0) ){
			i <- (istep-i_burnin) / thin
			ACCEPT[i] <- n_accept * ne_steps_per_tau_step/ istep / 2
			LOGPO[i] <- ll
			TAU[i] <- tau
			NE[i, ]<- ( ne ) 
			BETA[i,] <- beta
			#~  			logne <- log(ne ) 
			#~ 			GROWTHRATE[i, ] <-  diff ( (logne) ) / ( (logne[-res]) ) / dh 
			GROWTHRATE[i, ] <-  diff ( (ne) ) / ( (ne[-res]) ) / dh 
		}
		
		if ( istep %% thin == 0){
			if (!quiet){
				print( base::date() )
				print( istep )
				print( tau )
				print( beta )
				print( n_accept * ne_steps_per_tau_step/ istep / 2 )
			}
		}
	}
	
	GROWTHRATE <- cbind( GROWTHRATE, NA ) 
	## compute quantiles
	ne_ci <- t( sapply(1:ncol( NE ), function(i) quantile( NE[, i], prob = c( .025 , .5, .975 ), na.rm=T) ) )
	growthrate_ci <- t( sapply(1:ncol( GROWTHRATE ), function(i) quantile( GROWTHRATE[, i], prob = c( .025 , .5, .975 ), na.rm=T) ) )
	rv <- list( ne = NE
	  , growthrate = GROWTHRATE
	  , ne_ci = ne_ci
	  , growthrate_ci = growthrate_ci 
	  , tau = TAU
	  , time = rev(-tredat$heights ) 
	  , accept = ACCEPT
	  , beta = BETA 
	  , tredat = tredat
	  , gamma = NA
	  , control = control
	  , tre = tre
	  , tau_logprior = tau_logprior 
	  , mapfit = mapfit 
	  , logpo = LOGPO
	)
	
	class(rv) <- c('skygrowth.mcmc.covar', 'skygrowth.mcmc')
	rv
})}





########################################################################
##

#' Compute the reproduction number through time provided model fit and generation time
#'
#' @param fit A model fit 
#' @param gamma Per-capita death or recovery rate; equivalent to 1/generation time
#' @return Fitted object with $R attribute
#' @export
computeR <- function(fit, gamma){
  if (inherits(fit, "skygrowth.map")) return(computeR.skygrowth.map(fit,gamma))
  if (inherits(fit, "skygrowth.mcmc")) return(computeR.skygrowth.mcmc(fit,gamma))
}

computeR.skygrowth.map <- function(fit, gamma )
{
	stopifnot(inherits(fit, "skygrowth.map"))
	D <- 1/gamma
	# r = (R0 - 1) / D
	fit$gamma = gamma 
	fit$R <- fit$growth * D + 1
	fit
}

computeR.skygrowth.mcmc <- function(fit, gamma )
{
	stopifnot(inherits(fit, "skygrowth.mcmc"))
	D <- 1/gamma
	# r = (R0 - 1) / D
	fit$gamma = gamma 
	fit$R <- fit$growthrate * D + 1
	fit$R_ci <- t( sapply(1:ncol( fit$R ), function(i) quantile( fit$R[, i], prob = c( .025 , .5, .975 ), na.rm=T) ) )
	fit
}

## plots

#' Plot effective size through time 
#'
#' @param fit A fitted object (eg skygrowth.map or skygrowth.mcmc)
#' @param logy  If TRUE, the plot is returned with logarithmic y-axis
#' @param ggplot  If TRUE, returns a ggplot2 figure
#' @param ... Additional parameters are passed to ggplot or the base plotting function
#' @return A ggplot2 plot
#' @export
neplot <- function(fit, ggplot=TRUE, logy=TRUE, ... ){
  if (inherits(fit, "skygrowth.map")) return(neplot.skygrowth.map(fit, ggplot, logy, ... ))
  if (inherits(fit, "skygrowth.mcmc")) return(neplot.skygrowth.mcmc(fit, ggplot, logy, ...))
}

#' Plot growth rate of effective size through time 
#'
#' @param fit A fitted object (eg skygrowth.map)
#' @param logy  If TRUE, the plot is returned with logarithmic y-axis
#' @param ggplot  If TRUE, returns a ggplot2 figure
#' @param ... Additional parameters are passed to ggplot or the base plotting function
#' @return A ggplot2 plot
#' @export
growth.plot <- function(fit , ggplot=TRUE, logy=FALSE, ...){
  if (inherits(fit, "skygrowth.map")) return(growth.plot.skygrowth.map(fit,ggplot,logy,...))
  if (inherits(fit, "skygrowth.mcmc")) return(growth.plot.skygrowth.mcmc(fit,ggplot,logy,...))
}

#' Plot reproduction number through time 
#'
#' @param fit A fitted object (eg skygrowth.map)
#' @param gamma Value of gamma
#' @param ggplot If TRUE, returns a ggplot2 figure
#' @param ... Additional parameters are passed to ggplot or the base plotting function
#' @return A ggplot2 plot
#' @export
R.plot <- function(fit, gamma = NA, ggplot=TRUE,...){
  if (inherits(fit, "skygrowth.map")) return(R.plot.skygrowth.map(fit,gamma,ggplot))
  if (inherits(fit, "skygrowth.mcmc")) return(R.plot.skygrowth.mcmc(fit,gamma,ggplot))
}

neplot.skygrowth.map <- function( fit, ggplot=TRUE, logy=TRUE, ... )
{
	stopifnot(inherits(fit, "skygrowth.map"))
	ne <- fit$ne
	if ( 'ggplot2' %in% installed.packages()  & ggplot)
	{
		pldf <- data.frame( t = fit$time,  nemed = ne, nelb = fit$ne_ci[,1], neub = fit$ne_ci[,3] )
		pl <- ggplot2::ggplot( pldf, ggplot2::aes( x = t, y = nemed) , ...) + ggplot2::geom_line()+ ggplot2::ylab('Effective population size') + ggplot2::xlab('Time before most recent sample')
		pl <- pl + ggplot2::geom_ribbon( ggplot2::aes( ymin = nelb, ymax = neub), fill = 'blue', alpha = .2)
		if (logy) pl <- pl + ggplot2::scale_y_log10() 
		return(pl)
	} else{
		if (logy)
			plot( fit$time, ne, ylim=range(fit$ne_c[,1:3],na.rm=T),lwd =2, col = 'black', type = 'l', log='y', xlab='Time', ylab='Effective population size', ...)
		else
			plot( fit$time, ne, ylim=range(fit$ne_c[,1:3],na.rm=T),lwd =2, col = 'black', type = 'l', xlab='Time', ylab='Effective population size', ...)
		lines( fit$time, fit$ne_c[,1] , lty=3)
		lines( fit$time, fit$ne_c[,3] , lty=3)
		invisible(fit)
	}
}

growth.plot.skygrowth.map <- function( fit , ggplot=TRUE, logy=FALSE, ...)
{
	stopifnot(inherits(fit, "skygrowth.map"))
	if ( 'ggplot2' %in% installed.packages()  & ggplot)
	{
		pldf <- data.frame( t = fit$time, gr = fit$growth)
		pl <- ggplot2::ggplot( pldf, ggplot2::aes( x = t, y = gr), ... ) + ggplot2::geom_line() + ggplot2::ylab('Growth rate') + ggplot2::xlab('Time before most recent sample')
		if (logy) pl <- pl + ggplot2::scale_y_log10() 
		return(pl)
	} else{
		if (logy)
			plot( fit$time, fit$growth, lwd =2, col = 'black', type = 'l', log='y', xlab='Time', ylab='Growth rate',...)
		else
			plot( fit$time, fit$growth, lwd =2, col = 'black', type = 'l', xlab='Time', ylab='Growth rate', ...)
		
		invisible(fit)
	}
}

R.plot.skygrowth.map <- function(fit, gamma = NA , ggplot=TRUE, ...)
{
	stopifnot(inherits(fit, "skygrowth.map"))
	if ( is.na(fit$gamma) & is.na(gamma)) stop('Removal rate (gamma) must be supplied')
	if (is.na(gamma)) gamma <- fit$gamma
	fit <- computeR.skygrowth.map( fit, gamma )
	if ( 'ggplot2' %in% installed.packages()  & ggplot)
	{
		i <- 1:(length(fit$time)-1)
		pldf <- data.frame( t = fit$time[1:length(fit$R)],R = fit$R)
		ggplot2::ggplot( pldf, ggplot2::aes( x = t, y = R) , ...) + ggplot2::geom_line() + ggplot2::ylab('Reproduction number') + ggplot2::xlab('Time before most recent sample')
	} else{
		plot( fit$time, fit$R, lwd =2, col = 'black', type = 'l',xlab='Time', ylab='Reproduction number', ...)
		invisible(fit)	
	}
	
}

#' @export
plot.skygrowth.map <- function( x,  ... ){
	neplot( x, ...) 
}

neplot.skygrowth.mcmc <- function( fit, ggplot=TRUE, logy = TRUE , ... )
{
	stopifnot(inherits(fit, "skygrowth.mcmc"))
	ne <- fit$ne_ci
	if ( 'ggplot2' %in% installed.packages()  & ggplot)
	{
		pldf <- data.frame( t = fit$time, nelb = ne[,1], nemed = ne[,2], neub = ne[,3] )
		pl <- ggplot2::ggplot( pldf, ggplot2::aes( x = t, y = nemed), ... ) + ggplot2::geom_line() + ggplot2::geom_ribbon( ggplot2::aes( ymin = nelb, ymax = neub), fill = 'blue', alpha = .2) + ggplot2::ylab('Effective population size') + ggplot2::xlab('Time before most recent sample')
		if (logy) pl <- pl + ggplot2::scale_y_log10()
		return(pl)
	} else{
		if (logy)
			plot( fit$time, ne[,2], ylim=range(ne[,1:3],na.rm=T),lwd =2, col = 'black', type = 'l', log='y',xlab='Time', ylab='Effective population size', ...)
		else
			plot( fit$time, ne[,2], ylim=range(ne[,1:3],na.rm=T),lwd =2, col = 'black', type = 'l',xlab='Time', ylab='Effective population size', ...)
		lines( fit$time, ne[,1] , lty=3)
		lines( fit$time, ne[,3] , lty=3)
		invisible(fit)
	}
}

growth.plot.skygrowth.mcmc <- function( fit ,  ggplot=TRUE, logy = FALSE , ...)
{
	stopifnot(inherits(fit, "skygrowth.mcmc"))
	x <- fit$growthrate_ci
	if ( 'ggplot2' %in% installed.packages()  & ggplot)
	{
		pldf <- data.frame( t = fit$time, lb = x[,1], med = x[,2], ub = x[,3] )
		pl <- ggplot2::ggplot( pldf, ggplot2::aes( x = t, y = med), ... ) + ggplot2::geom_line() + ggplot2::geom_ribbon( ggplot2::aes( ymin = lb, ymax = ub), fill = 'blue', alpha = .2) + ggplot2::ylab('Growth rate') + ggplot2::xlab('Time before most recent sample')
		if (logy) pl <- pl + ggplot2::scale_y_log10()
		return(pl) 
	} else{
		if (logy)
			plot( fit$time, x[,2], ylim=range(x[,1:3],na.rm=T),lwd =2, col = 'black', type = 'l', log='y', xlab='Time', ylab='Growth rate', ...)
		else
			plot( fit$time, x[,2], ylim=range(x[,1:3],na.rm=T),lwd =2, col = 'black', type = 'l',xlab='Time', ylab='Growth rate', ...)
		#
		lines( fit$time, x[,1] , lty=3)
		lines( fit$time, x[,3] , lty=3)
		invisible(fit)
	}
}

R.plot.skygrowth.mcmc <- function(fit, gamma = NA, ggplot=TRUE )
{
	stopifnot(inherits(fit, "skygrowth.mcmc"))
	fit <- computeR.skygrowth.mcmc( fit, gamma )
	x <- fit$R_ci
	if ( 'ggplot2' %in% installed.packages()  & ggplot)
	{
		if ( is.na(fit$gamma) & is.na(gamma)) stop('Removal rate (gamma) must be supplied')
		if (is.na(gamma)) gamma <- fit$gamma
		
		
		pldf <- data.frame( t = fit$time, lb = x[,1], med = x[,2], ub = x[,3] )
		ggplot2::ggplot( pldf, ggplot2::aes( x = t, y = med) ) + ggplot2::geom_line() + ggplot2::geom_ribbon( ggplot2::aes( ymin = lb, ymax = ub), fill = 'blue', alpha = .2) + ggplot2::ylab('Reproduction number') + ggplot2::xlab('Time before most recent sample')
	} else{
		plot( fit$time, x[,2], ylim=range(x[,1:3],na.rm=T),lwd =2, col = 'black', type = 'l',xlab='Time', ylab='Reproduction number')
		
		lines( fit$time, x[,1] , lty=3)
		lines( fit$time, x[,3] , lty=3)
		invisible(fit)
	}
}

#' @export
plot.skygrowth.mcmc <- function( x, ... ){
	neplot( x, ...) 
}


## print and summary methods

#' @export
print.skygrowth.mcmc <- function(x, ...)
{
	stopifnot(inherits(x, "skygrowth.mcmc"))
	cat( 'Bayesian non parametric moving average phylodynamic model\n')
	cat(paste( 'Effective population size bins:', ncol(x$ne), '\n'))
	cat(paste( 'Iterations:', ncol(x$ne), '\n'))
	cat('Effective population size at last iteration:\n')
	i <- 1:length(x$time) #seq(1, length(x$tim), le = 5)
	print( data.frame( Time=x$time[ i] , Ne=x$ne[i] ) )
	invisible( x )
}

#' @export
summary.skygrowth.mcmc <- function(object, ...){
	print.skygrowth.mcmc( object )
}

#' @export
print.skygrowth.map <- function(x, ...)
{
	stopifnot(inherits(x, "skygrowth.map"))
	cat( 'Maximum a posteriori non parametric moving average phylodynamic model\n')
	cat(paste( 'Effective population size bins:', length(x$ne), '\n'))
	cat('Effective population size at last iteration:\n')
	i <- 1:length(x$time) #seq(1, length(x$tim), le = 5)
	print( data.frame( Time=x$time[ i] , Ne=x$ne[i] ) )
	cat('Drift pararameter (tau):\n')
	print( x$tau )
	invisible( x )
}

#' @export
summary.skygrowth.map <- function(object, ...){
	print( object )
}


