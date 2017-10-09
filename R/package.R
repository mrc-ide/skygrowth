#' @name skygrowth
#' @title Phylodynamic inference with nonparametric growth rate models
#' 
#' @description This package provides methods for inferring effective population size for time-scaled 
#' phylogenies using Bayesian MCMC and maximum a posteriori. 
#' 
#' The main functions of the package are:
#' \itemize{
#' 
#' \item phylo.mapma and phylo.bnpma provide maximum and a posteriori methods respectively for 
#' inferring effective population size through time given a dated genealogy. 
#' \item phylo.mapma.covars and phylo.bnpma.covars allow you to include a data.frame with 
#' time-varying covariates which may improve the inference. Covariates should be included if 
#' they are likely to be associated with the growth rate of effective population size. 
#' }
#' 
#' @author Erik M Volz and Xavier Didelot
#' Maintainer Erik M Volz \email{erik.volz@gmail.com}
#' 
#' @references Erik M Volz and Xavier Didelot, Manuscript in prepration. 
#' @seealso https://github.com/mrc-ide/skygrowth

#' @importFrom Rcpp evalCpp
#' @useDynLib skygrowth
#' @import stats
#' @import ape
#' @importFrom graphics plot lines
#' @importFrom utils installed.packages tail
NULL