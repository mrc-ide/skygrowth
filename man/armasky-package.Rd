\name{armasky-package}
\alias{armasky-package}
\alias{armasky}
\docType{package}
\title{
  Phylodynamic inference with nonparametric moving average models
}
\description{
  This package provides methods for inferring effective population size for time-scaled phylogenies using Bayesian MCMC and maximum a priori. The transformed effective population size is modeled as a stochastic process with autoregressive and moving average components. 
}
\details{
  phylo.mapma and phylo.bnpma provide maximum and a posteriori methods respectively for inferring effective population size through time given a dated genealogy. 
  phylo.mapma.covars and phylo.bnpma.covars allow you to include a data.frame with time-varying covariates which may improve the inference. Covariates should be included if they are likely to be associated with the growth rate of effective population size. 
}
\author{
Erik M Volz and Xavier Didelot

Maintainer: Erik M Volz <erik.volz@gmail.com>
}
\references{
  Erik M Volz and Xavier Didelot, Manuscript in prepration. 
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{ package }
\seealso{
  https://github.com/mrc-ide/armasky
}
\examples{
  
}
