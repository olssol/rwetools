# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Get Moment Constraints for Logistic Regression
NULL

#' Get Moment Constraints for balance in covariates
NULL

#' Get Moment Constraints
#'
#' return N row k column matrix
#' Row:   subject
#' Column: Moment
#'
NULL

#' Get PS
#'
#' @export
c_get_ps <- function(beta, mat_x, tol = 1e-6) {
    .Call(`_rwetools_c_get_ps`, beta, mat_x, tol)
}

#' Test Rcpp function
#'
#'
#' @param test test parameter
#'
#' @export
crtTest <- function(test) {
    .Call(`_rwetools_crtTest`, test)
}

#' Match by nearest neighbor
#'
#' Match subjects in group candidate with subject in the target group
#'
#' @param target  vector of data from the target group
#' @param candidates vector of data from the candidate group
#' @param ratio  1:ratio match
#'
cMatch <- function(target, candidate, ratio) {
    .Call(`_rwetools_cMatch`, target, candidate, ratio)
}

#' Get ALL Moment Constraints
#'
#' return N row k column matrix
#' Row:   subject
#' Column: Moment
#'
#'
c_ps_gmm_g <- function(beta, mat_grp_x, att = FALSE) {
    .Call(`_rwetools_c_ps_gmm_g`, beta, mat_grp_x, att)
}

#' Get Gbar and Sigma
#'
#' return a List
#'
#'
c_ps_gmm_gbar <- function(beta, mat_x, ps, mgrp, mz, mr, n3) {
    .Call(`_rwetools_c_ps_gmm_gbar`, beta, mat_x, ps, mgrp, mz, mr, n3)
}

#' Get Gbar and Sigma
#'
#' return a List
#'
#' @export
c_ps_gmm_sigma <- function(beta, mat_x, ps, mgrp, mz, mr, n3) {
    .Call(`_rwetools_c_ps_gmm_sigma`, beta, mat_x, ps, mgrp, mz, mr, n3)
}

#' Get Derivative of Moment Constraints
#'
#' return k row p column matrix
#' Row:    1..k moment functions
#' Column: 1..p coefficients
#'
#'
c_ps_gmm_dg <- function(beta, mat_grp_x, att = FALSE) {
    .Call(`_rwetools_c_ps_gmm_dg`, beta, mat_grp_x, att)
}

