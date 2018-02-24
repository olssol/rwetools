#' R tools for synthesizing real-world evidence
#'
#' @docType   package
#' @name      rwetools-package
#' @aliases   rwetools
#' @useDynLib rwetools, .registration = TRUE
#'
#' @importFrom rstan     sampling extract stanc rstan_options traceplot stan_rhat
#' @importFrom grDevices colors
#' @importFrom graphics  axis box legend lines par plot points text
#' @importFrom loo       extract_log_lik loo
#' @importFrom parallel  detectCores
#' @importFrom MCMCpack  rdirichlet
#' @importFrom mvtnorm   rmvnorm
#'
#' @import stats
#' @import Rcpp
#' @import methods
#' @import ggplot2
#'
#' @description
#'
#' This package contains the functions for synthesizing real-world evidence in
#' single-arm medical device studies.
NULL



#' Parameters for simulating data
#'
#' @name simupara
#'
#' @param nPat     number of patients
#' @param muCov    mean vector of covariates
#' @param sdCov    standard deviation vector of covariates
#' @param corCov   correlation of covariates
#' @param regCoeff regression coefficients
#' @param type     distributions of the random error
#' @param ysig     standard error of the random error
#' @param skew.n   parameter of negative bionomial distribution
#' @param skew.p   parameter of negative bionomial distribution
#' @param b0       intercept in regession model
#' @param bin.mu      mean of the binary outcomes used to compute b0
#' @param formula.z   formula of the treatement assignment model. No intercept term.
#' @param formula.y   formula of the outcome model. No intercept term.
#'
#'
NULL
