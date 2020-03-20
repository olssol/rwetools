##-----------------------------------------------------------------------------
##     FUNCTIONS RELATED TO USING MARGINAL STATISTICS TO
##     1) SAMPLE PATIENTS
##     2) CALCULATE LIKELIHOOD
##     3)
##-----------------------------------------------------------------------------

#' Sample patient based on summary statistics
#'
#' @param target_stats target summary statistics
#' @param dta_ext external data
#' @param method sampling methods
#' @param n     number of patients to be drawn in each try
#' @param weights weights based on likelihood
#'
#' @export
#'
rwe_margin_sample <- function(target_stats, dta_ext,
                              method = c("genetic", "sa",
                                         "naive", "weighted", "ps"),
                              n_min = 300, max_steps = 10000, epsilon = NULL,
                              seed = NULL, ...) {

    method <- match.arg(method)

    if (!is.null(seed))
        set.seed(seed)

    rst <- tkCallFun(c("private_margin_sample_", method),
                     target_stats, dta_ext, n_min, max_steps, epsilon,
                     ...)

    ## return
    best_selected <- rst$best_selected
    best_val      <- rst$best_val

    list(selected = dta_ext[1 == best_selected, ],
         ext_ind  = best_selected,
         distance = best_val$distance,
         stats    = best_val$stats,
         target   = lst_stats)
}


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

    ll      <- apply(rst, 1, sum)
    weights <- exp(ll)
    weights <- weights / sum(weights)

    cbind(ll = ll, weights = weights)
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
            discrete   = tkExtractStats(cur_y, xlev = cur_s$values),
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
rwe_margin_stat_diff <- function(lst_stats, y, type = c("max", "sum"),
                                 epsilon = NULL, ...) {

    type        <- match.arg(type)
    lst_stats_y <- rwe_extract_stats(lst_stats, y)

    ## all differences
    rst <- NULL
    for (i in seq_len(length(lst_stats))) {
        cur_diff <- private_stat_diff(lst_stats[[i]],
                                      lst_stats_y[[i]])
        rst      <- c(rst, cur_diff)
    }

    ## distance summary
    dist <- switch(type,
                   max = max(abs(rst)),
                   sum = sum(abs(rst)),
                   9999)

    ## acceptable or not
    yn <- FALSE
    if (!is.null(epsilon))
        yn <- dist < epsilon

    list(stats      = lst_stats_y,
         diff       = rst,
         distance   = dist,
         acceptable = yn)
}



## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
##             PRIVATE FUNCTIONS
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

## Simulate a covariate based on summary statistics
##
##
## @return column vector with n patients
##
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

## Calculate log likelihood
##
## Calculate likelihood of a covariate based on its summary statistics
##
## @param y vector of data
## @param lst_stats summary statistics
##
## @return column vector with n patients
##
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

## Quantile matrix
##
## Get matrix data of quantiles
##
## @return matrix with 3 columns. 1: probabilities, 2: lower bound,
## 3: upper bound
##
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

## Difference in summary statistics
private_stat_diff <- function(stat0, stat1) {

    stopifnot(stat0$type == stat1$type)

    rst <- switch(stat0$type,
                  discrete = {
                      stopifnot(identical(stat0$values,
                                          stat1$values))
                      diff <- stat0$probs - stat1$probs
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

## A tool function for getting distance
f_dist <- function(target_stats, dta_ext, ind_selected, ...) {
    cur_sel <- which(1 == ind_selected)

    if (0 == length(cur_sel))
        return(NULL)

    rst_dist <- rwe_margin_stat_diff(target_stats,
                                     y = dta_ext[cur_sel, ],
                                     ...)
    rst_dist
}

## Naive Sampling from RWD
private_margin_sample_naive <- function(target_stats, dta_ext, n_min = 300,
                                        max_steps = 10000, ...,
                                        weighted = FALSE, ps = FALSE) {

    f_select <- function() {
        cur_n   <- sample(n_min:n_max, 1)

        if (!ps) {
            cur_inx <- sample(1:n_max, cur_n, replace = FALSE,
                              prob = weights)
        } else {
            simu_target <- rwe_margin_simu(n = cur_n, target_stats)
            cur_inx     <- rwe_match_ps(simu_target,
                                        dta_ext,
                                        n_match = 1,
                                        v.covs  = names(target_stats))
        }

        cur_selected           <- rep(0, n_max)
        cur_selected[cur_inx]  <- 1
        cur_val                <- f_dist(target_stats, dta_ext, cur_selected,
                                         ...)

        list(selected = cur_selected,
             val      = cur_val)
    }

    ## total patients
    n_max <- nrow(dta_ext)

    ## set weights
    if (weighted) {
        weights <- rwe_margin_ll(target_stats, dta_ext)[, "weights"]
    } else {
        weights <- rep(1/n_max, n_max)
    }

    ## initiate
    cur_rst       <- f_select()
    best_selected <- cur_rst$selected
    best_val      <- cur_rst$val

    k <- 1
    while (k < max_steps & !best_val$acceptable) {
        cur_rst <- f_select()

        if (best_val$distance > cur_rst$val$distance) {
            best_selected <- cur_rst$selected
            best_val      <- cur_rst$val

            cat("-----", k, sum(best_selected), best_val$distance, "\n")
        }

        k <- k + 1
    }

    list(best_selected = best_selected,
         best_val      = best_val)
}

## Weighted Sampling from RWD
private_margin_sample_weighted <- function(...) {
    private_margin_sample_naive(..., weighted = TRUE, ps = FALSE)
}


## PS Sampling from RWD
private_margin_sample_ps <- function(...) {
    private_margin_sample_naive(..., weighted = FALSE, ps = TRUE)
}

## Simulated Annealing for Sampling from RWD
private_margin_sample_sa <- function(target_stats, dta_ext, n_min = 300,
                                     max_steps = 10000, ..., alpha = 0.1,
                                     weighted = TRUE, swap_size = 5,
                                     p_extend = 0.5) {

    ## select samples from available or selected pool
    f_swap <- function(ind_selected) {
        ## extend or not
        to_extend <- rbinom(1, 1, p_extend)

        if (0 == to_extend | all(1 == ind_selected)) {
            ## remove current selected
            inx <- which(1 == ind_selected)
            rst <- sample(inx, size = min(swap_size, length(inx) - n_min),
                          prob = 1 - weights[inx])
            ind_selected[rst] <- 0
        } else {
            ## add to selected patients
            inx <- which(0 == ind_selected)
            rst <- sample(inx, size = min(swap_size, length(inx)),
                          prob = weights[inx])
            ind_selected[rst] <- 1
        }
        ind_selected

    }

    ## total patients
    n_max <- nrow(dta_ext)

    ## set weights
    if (weighted) {
        weights <- rwe_margin_ll(lst_stats, dta_ext)[, "weights"]
    } else {
        weights <- rep(1/nrow(dta_ext), nrow(dta_ext))
    }

    ## indicator of selected patients
    cur_n                 <- sample(n_min:n_max, 1)
    cur_inx               <- sample(1:n_max, cur_n, replace = FALSE,
                                    prob = weights)
    cur_selected          <- rep(0, n_max)
    cur_selected[cur_inx] <- 1
    cur_val               <- f_dist(target_stats, dta_ext, cur_selected, ...)

    ## initiate best values
    best_selected <- cur_selected
    best_val      <- cur_val

    k <- 1
    while (k < max_steps & !best_val$acceptable) {
        temp_T        <- alpha * k / max_steps
        next_selected <- f_swap(cur_selected)
        next_val      <- f_dist(target_stats, dta_ext, next_selected, ...)

        diff_val      <- next_val$distance - cur_val$distance
        p_update      <- exp( - diff_val / temp_T)

        if (p_update >= runif(1)) {
            cur_selected <- next_selected
            cur_val      <- next_val
        }

        if (best_val$distance > cur_val$distance) {
            best_selected <- cur_selected
            best_val      <- cur_val

            cat("-----", k, sum(best_selected), best_val$distance, "\n")
        }

        k <- k + 1
    }

    list(best_selected = best_selected,
         best_val      = best_val)
}

## Genetic algorithm for sampling
private_margin_sample_genetic <- function(target_stats, dta_ext, n_min = 300,
                                          max_steps = 10000, ..., monitor = FALSE) {

    fitness <- function(ind_selected, ...) {
        rst <- f_dist(target_stats, dta_ext, ind_selected, ...)
        return(-rst$distance)
    }

    rst <- ga(type = "binary",
              fitness = fitness, nBits = nrow(dta_ext),
              keepBest = TRUE, maxiter = max_steps,
              monitor = monitor, ...)

    best_selected <- as.numeric(rst@solution)
    best_val      <- f_dist(target_stats, dta_ext, best_selected, ...)

    list(best_selected = best_selected,
         best_val      = best_val)
}

## Random forest for Sampling from RWD
private_margin_sample_rf <- function(target_stats, dta_ext, n_min = 300,
                                     max_steps = 10000, ...) {

    ## current distance
    f_dis <- function(ind_selected) {
        cur_sel  <- which(1 == ind_selected)
        rst_dist <- rwe_margin_stat_diff(lst_stats,
                                         y = dta_ext[cur_sel, ],
                                         type = type)
        rst_dist
    }

    ## initial random seed
    if (!is.null(seed))
        set.seed(seed)

    ## set weights
    if (weighted) {
        weights <- rwe_margin_ll(lst_stats, dta_ext)[, "weights"]
    } else {
        weights <- rep(1/nrow(dta_ext), nrow(dta_ext))
    }

    ## indicator of selected patients
    cur_selected <- rep(1, nrow(dta_ext))
    cur_val      <- f_dis(cur_selected)
    best_selected <- cur_selected
    best_val      <- cur_val

    ind_tree <- rep(1, nrow(dta_ext))
    k        <- 1
    while (sum(cur_selected) > n_min
           & best_val$distance > epsilon
           & k < max_steps) {
               cur_size <- sample(1:swap_size, 1)
               inx_cand <- which(1 == cur_selected)
               cur_inx  <- sample(inx_cand,
                                  size = cur_size,
                                  prob = 1 - weights[inx_cand])
               next_selected          <- cur_selected
               next_selected[cur_inx] <- 0
               next_val               <- f_dis(next_selected)
               diff_val               <- next_val$distance - cur_val$distance

               if (diff_val < 0) {
                   cur_selected <- next_selected
                   cur_val      <- next_val

                   if (best_val$distance > cur_val$distance) {
                       best_selected <- cur_selected
                       best_val      <- cur_val

                       cat("iter", k, "distance", best_val$distance, "\n")
                   }
               }

               ## keep node
               ind_tree[cur_inx] <- 0
               ## update steps
               k                 <- k + 1
           }
}



