##-----------------------------------------------------------------------------
##     FUNCTIONS RELATED TO USING MARGINAL STATISTICS TO
##     1) SAMPLE PATIENTS
##     2) CALCULATE LIKELIHOOD
##     3)
##-----------------------------------------------------------------------------

#' Simulate covariates based on summary statistics
#'
#'
#' @return data frame with n patients, each column represents a covariate and
#'     each row represents a patient
#'
#' @export
#'
rwe_margin_simu <- function(lst_stats, n = 500) {
    rst <- lapply(lst_stats, private_margin_simu)
    rst <- data.frame(rst)
}

#' Calculate log-likelihood based on summary statistics
#'
#'
#' @param y matrix of data, each column represents a covariate and each row
#'     represents a patient
#' @return data frame for all y
#'
#' @export
#'
rwe_margin_ll <- function(lst_stats, y) {
    rst <- NULL
    for (i in seq_len(length(lst_stats))) {
        cur_stats <- lst_stats[[i]]
        cur_y     <- y[[names(lst_stats)[i]]]
        cur_ll    <- private_margin_ll(cur_stats, cur_y)
        rst       <- cbind(rst, cur_ll)
    }

    apply(rst, 1, sum)
}


#' Extract summary statistics
#'
#' Extract summary statistics from a data frame based on an existing summary
#' statistics list
#'
#' @return list of summary statistics
#'
#' @export
#'
rwe_extract_stats <- function(lst_stats, y) {
    rst <- list()
    for (i in seq_len(length(lst_stats))) {
        cur_v   <- names(lst_stats)[i]
        cur_y   <- y[[cur_v]]
        cur_s   <- lst_stats[[i]]

        cur_rst <- switch(
            cur_s$type,
            discrete   = tkExtractStats(cur_y,
                                        xlev = cur_s$values),
            continuous = tkExtractStats(cur_y, type = "continuous"),
            quants     = tkExtractStats(cur_y, type = "quants",
                                        quants = cur_s$quants[, 1]))
        rst[[cur_v]] <- cur_rst
    }
    rst
}

#' Calculate differences in summary statistics
#'
#' @return vector of differences in summary statistics
#'
#' @export
#'
rwe_margin_stat_diff <- function(lst_stats, y) {
    lst_stats_y <- rwe_extract_stats(lst_stats, y)

    rst <- NULL
    for (i in seq_len(length(lst_stats))) {
        cur_diff <- private_stat_diff(lst_stats[[i]],
                                      lst_stats_y[[i]])
        rst      <- c(rst, cur_diff)
    }

    rst
}



## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
##             PRIVATE FUNCTIONS
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

#' Simulate a covariate based on summary statistics
#'
#'
#' @return column vector with n patients
#'
private_margin_simu <- function(lst_stats, n = 500) {
    ## sample based on quantiles
    fsmp_quants <- function(n, range, quants) {
        qmat <- private_quantile_mat(range, quants)
        rstq <- sample(seq_len(nrow(qmat)),
                       size = n, replace = TRUE, prob = qmat[, 1])
        rst  <- sapply(rstq, function(x) {
            runif(1, min = qmat[x, 2], max = qmat[x, 3])
        }, simplify = TRUE)

        rst
    }

    rst <- switch(lst_stats$type,
                  discrete = {
                      smps <- sample(seq_len(length(lst_stats$values)),
                                     size    = n,
                                     replace = TRUE,
                                     prob    = lst_stats$probs)
                      factor(lst_stats$values[smps])
                  },
                  continuous = rnorm(n,
                                     lst_stats$mean,
                                     lst_stats$sd),
                  quants = fsmp_quants(n,
                                       lst_stats$range,
                                       lst_stats$quants),
                  NULL)
    rst
}

#' Calculate log likelihood
#'
#' Calculate likelihood of a covariate based on its summary statistics
#'
#' @param y vector of data
#' @param lst_stats summary statistics
#'
#' @return column vector with n patients
#'
private_margin_ll <- function(lst_stats, y) {
    ## ll based proportions
    f_dis <- function(y, values, probs) {
        lp <- log(probs)
        sapply(y, function(x) {
            lp[which(x == values)]
        }, simplify = TRUE)
    }

    rst <- switch(lst_stats$type,
                  discrete   = f_dis(y, lst_stats$values, lst_stats$probs),
                  continuous = dnorm(y, lst_stats$mean, lst_stats$sd,
                                     log = TRUE),
                  NULL)
    rst
}

#' Quantile matrix
#'
#' Get matrix data of quantiles
#'
#' @return matrix with 3 columns. 1: probabilities, 2: lower bound,
#' 3: upper bound
#'
private_quantile_mat <- function(range, quants) {
    qmat  <- NULL
    lastq <- 0
    lastx <- range[1]
    for (i in seq_len(nrow(quants))) {
        curq <- quants[i, 1]
        curx <- quants[i, 2]

        cmat  <- c(curq - lastq, lastx, curx)
        qmat  <- rbind(qmat, cmat)

        lastq <- curq
        lastx <- curx
    }

    cmat <- c(1 - lastq, lastx, range[2])
    qmat <- rbind(qmat, cmat)
}

#' Difference in summary statistics
#'
#'
private_stat_diff <- function(stat0, stat1) {

    stopifnot(stat0$type == stat1$type)

    rst <- switch(stat0$type,
                  discrete = {
                      stopifnot(identical(stat0$values,
                                          stat1$values))
                      stat0$probs - stat1$probs
                  },
                  continuous = {
                      c(stat0$mean - stat1$mean,
                        stat0$sd   - stat1$sd)
                  },
                  quants = {
                      stopifnot(identical(stat0$quants[, 1],
                                          stat1$quants[, 1]))
                      stat0$quants[, 2] - stat1$quant[, 2]
                  },
                  ... = NULL)
}
