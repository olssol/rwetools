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
         diff     = best_val$diff,
         stats    = best_val$stats,
         target   = target_stats)
}


#' Simulate covariates based on summary statistics
#'
#'
#' @return data frame with n patients, each column represents a covariate and
#'     each row represents a patient
#'
#' @export
#'
rwe_margin_simu <- function(target_stats, n = 500) {
    rst <- lapply(target_stats, private_margin_simu)
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
rwe_margin_ll <- function(target_stats, y) {
    rst <- NULL
    for (i in seq_len(length(target_stats))) {
        cur_stats <- target_stats[[i]]
        cur_y     <- y[[names(target_stats)[i]]]
        cur_ll    <- private_margin_ll(cur_stats, cur_y)
        rst       <- cbind(rst, cur_ll)
    }

    ll      <- apply(rst, 1, sum)
    weights <- exp(ll)
    weights <- weights / sum(weights)

    cbind(ll = ll, weights = weights)
}

#' Calculate PS score based on simulated target data
#'
#'
#'
#' @export
#'
rwe_margin_ps <- function(target_stats, dta_ext, n_cur = NULL,
                          reps = 10, ps_thresh = 0.01) {

    n_ext <- nrow(dta_ext)
    if (is.null(n_cur))
        n_cur <- n_ext

    rst_ps <- NULL
    for (i in 1:reps) {
        dta_cur         <- rwe_margin_simu(n = n_cur, target_stats)
        v_covs          <- colnames(dta_cur)
        dta_ext         <- dta_ext[, v_covs]
        dta_cur$grp_tmp <- 1
        dta_ext$grp_tmp <- 0
        dta             <- rbind(dta_ext, dta_cur)

        ## get ps
        dta_ps <- rwePS(data = dta, v.grp = "grp_tmp", v.covs = v_covs)

        ## return dta_ext PS
        cur_ps <- dta_ps$data[["_ps_"]][1:n_ext]

        cur_ps[which(cur_ps < ps_thresh)] <- ps_thresh
        rst_ps <- cbind(rst_ps, cur_ps)
    }

    ## return
    apply(rst_ps, 1, mean)
}

#' Calculate entropy score
#'
#' Get entropy balancing score
#'
#' @export
#'
rwe_margin_entropy <- function(target_stats, dta_ext,
                               tol = 1e-8, max_it = 10000,
                               print_level = 0) {

    line.searcher <- function (Base.weight, Co.x, Tr.total,
                               coefs, Newton, ss) {
        weights.temp <- c(exp(Co.x %*% (coefs - (ss * Newton))))
        weights.temp <- weights.temp * Base.weight
        Co.x.agg     <- c(weights.temp %*% Co.x)
        maxdiff      <- max(abs(Co.x.agg - Tr.total))
        return(maxdiff)
    }

    x_m         <- rwe_extract_stats_covx(target_stats, dta_ext)
    covx        <- cbind(1, x_m$covx)
    n_ext       <- nrow(covx)
    tr_total    <- c(n_ext, x_m$constraints * n_ext)
    base_weight <- rep(1, n_ext)
    coefs       <- rep(0, ncol(covx))

    ## optimize
    converged   <- FALSE
    for (iter in 1:max_it) {
        weights_temp <- c(exp(covx %*% coefs))
        weights_ebal <- weights_temp * base_weight
        covx_agg     <- c(weights_ebal %*% covx)
        gradient     <- covx_agg - tr_total

        if (1 == iter)
            init_stats <- covx_agg


        if (max(abs(gradient)) < tol) {
            converged <- TRUE
            break
        }

        if (print_level >= 2) {
            cat("Iteration", iter, "maximum deviation is =",
                format(max(abs(gradient)), digits = 4), "\n")
        }

        hessian  <- t(covx) %*% (weights_ebal * covx)
        Coefs    <- coefs
        newton   <- solve(hessian, gradient)
        coefs    <- coefs - newton

        loss.new <- line.searcher(Base.weight = base_weight,
                                  Co.x = covx,
                                  Tr.total = tr_total,
                                  coefs = coefs,
                                  Newton = newton, ss = 1)

        loss.old <- line.searcher(Base.weight = base_weight,
                                  Co.x = covx,
                                  Tr.total = tr_total,
                                  coefs = Coefs,
                                  Newton = newton, ss = 0)

        if (is.nan(loss.new) |
            loss.old <= loss.new) {
            ss.out <- optimize(line.searcher, lower = 1e-05,
                               upper = 1, maximum = FALSE,
                               Base.weight = base_weight,
                               Co.x = covx, Tr.total = tr_total, coefs = Coefs,
                               Newton = newton)

            coefs <- Coefs - ss.out$minimum * solve(hessian, gradient)
        }
    }

    ## return
    weights <- weights_ebal / n_ext
    stats   <- rbind(as.numeric(init_stats)[-1] / n_ext,
                     x_m$constraints,
                     as.numeric(weights %*% covx)[-1])
    rownames(stats) <- c("raw", "target", "adjusted")

    list(maxdiff      = max(abs(gradient)),
         coefs        = coefs,
         weights      = weights,
         stats        = stats,
         converged    = converged)
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
rwe_extract_stats <- function(target_stats, y, ...) {
    rst <- list()
    for (i in seq_len(length(target_stats))) {
        cur_v   <- names(target_stats)[i]
        cur_y   <- y[[cur_v]]
        cur_s   <- target_stats[[i]]

        cur_rst <- switch(
            cur_s$type,
            discrete   = tk_extract_stats(cur_y, xlev = cur_s$values, ...),
            continuous = tk_extract_stats(cur_y, type = "continuous", ...),
            quants     = tk_extract_stats(cur_y, type = "quants", ...,
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
rwe_margin_stat_diff <- function(target_stats, y, type = c("max", "sum"),
                                 epsilon = NULL, ...) {

    type           <- match.arg(type)
    target_stats_y <- rwe_extract_stats(target_stats, y, ...)

    ## all differences
    rst <- NULL
    for (i in seq_len(length(target_stats))) {
        cur_diff <- private_stat_diff(target_stats[[i]],
                                      target_stats_y[[i]])
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

    list(stats      = target_stats_y,
         diff       = rst,
         distance   = dist,
         acceptable = yn)
}

#' Get estimating equation
#'
#' Get moment equations based on existing statistics
#'
#' @return list of summary statistics
#'
#' @export
#'
rwe_extract_stats_covx <- function(target_stats, y,
                                   cont_stat = c("mean", "ex2")) {

    ## continuous covariates summary statistics
    cont_stat <- match.arg(cont_stat, several.ok = TRUE)

    f_disc <- function(cur_y, stats) {

        xlev  <- stats$values
        probs <- stats$probs

        c_rst <- NULL
        c_m   <- NULL
        for (i in seq_len(length(xlev) - 1)) {
            x     <- xlev[i]
            c_rst <- cbind(c_rst,
                           as.numeric(x == cur_y))
            c_m   <- c(c_m, probs[i])
        }

        list(c_rst, c_m)
    }

    f_cont <- function(cur_y, stats, cont_stat) {
        c_rst <- NULL
        c_m   <- NULL
        for (i in cont_stat) {
            cur_v <- stats[[i]]
            if (is.null(cur_v))
                next;

            if ("mean" == i) {
                y <- cur_y
            } else if ("ex2" == i) {
                y <- y^2
            }

            c_rst <- cbind(c_rst, y)
            c_m   <- c(c_m, cur_v)
        }

        list(c_rst, c_m)
    }

    f_quan <- function(cur_y, stats) {
        quants <- stats$quants

        c_m   <- quants[, "quants"]
        c_rst <- NULL
        for (i in seq_len(nrow(quants))) {
            qt    <- quants[i, "x_quants"]
            c_rst <- cbind(c_rst,
                           as.numeric(cur_y <= qt))
        }

        list(c_rst, c_m)
    }

    rst_x <- NULL
    rst_m <- NULL
    for (i in seq_len(length(target_stats))) {
        cur_v   <- names(target_stats)[i]
        cur_s   <- target_stats[[i]]
        cur_y   <- y[[cur_v]]

        if (!is.null(cur_s$scale))
            cur_y <- cur_y * cur_s$scale

        cur_rst <- switch(
            cur_s$type,
            discrete   = f_disc(cur_y, cur_s),
            continuous = f_cont(cur_y, cur_s, cont_stat),
            quants     = f_quan(cur_y, cur_s)
        )

        rst_x <- cbind(rst_x, cur_rst[[1]])
        rst_m <- c(rst_m, cur_rst[[2]])
    }

    ## return
    list(covx        = rst_x,
         constraints = rst_m)
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
private_margin_simu <- function(target_stats, n = 500) {
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

    rst <- switch(target_stats$type,
                  discrete = {
                      smps <- sample(seq_len(length(target_stats$values)),
                                     size    = n,
                                     replace = TRUE,
                                     prob    = target_stats$probs)
                      factor(target_stats$values[smps])
                  },
                  continuous = rnorm(n,
                                     target_stats$mean,
                                     target_stats$sd),
                  quants = fsmp_quants(n,
                                       target_stats$range,
                                       target_stats$quants),
                  NULL)
    rst
}

## Calculate log likelihood
##
## Calculate likelihood of a covariate based on its summary statistics
##
## @param y vector of data
## @param target_stats summary statistics
##
## @return column vector with n patients
##
private_margin_ll <- function(target_stats, y) {
    ## ll based proportions
    f_dis <- function(y, values, probs) {
        lp <- log(probs)
        sapply(y, function(x) {
            lp[which(x == values)]
        }, simplify = TRUE)
    }

    rst <- switch(target_stats$type,
                  discrete   = f_dis(y, target_stats$values, target_stats$probs),
                  continuous = dnorm(y, target_stats$mean, target_stats$sd,
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
                      diff[-1]
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
        weights <- rep(1 / n_max, n_max)
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

            ## cat("-----", k, sum(best_selected), best_val$distance, "\n")
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
        weights <- rwe_margin_ll(target_stats, dta_ext)[, "weights"]
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

            ## cat("-----", k, sum(best_selected), best_val$distance, "\n")
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
        rst_dist <- rwe_margin_stat_diff(target_stats,
                                         y = dta_ext[cur_sel, ],
                                         type = type)
        rst_dist
    }

    ## initial random seed
    if (!is.null(seed))
        set.seed(seed)

    ## set weights
    if (weighted) {
        weights <- rwe_margin_ll(target_stats, dta_ext)[, "weights"]
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
