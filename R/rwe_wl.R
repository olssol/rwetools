#' PS-Integrated Weighted Likelihood Estimation
#'
#' For all stratum. Variance estimated by bootstrap or Jack-Knife. Works for
#' one- or two-arm studies when there is one RWD
#'
#' @param data class DWITHPS data frame
#' @param ... parameters for \code{rweWL}
#' @param bs.n number of bootstraps
#' @param lambdas power parameter without standardization by ns0
#' @param m.var method to get variance: jackknife or bootstrap
#' @param seed random seed
#'
#' @export
#'
rwePsWl <- function(data, lambdas, v.outcome = "Y", m.var = c("jk", "bs"),
                    bs.n = 1000, seed = NULL, ...) {

    stopifnot(v.outcome %in% colnames(data));

    m.var <- match.arg(m.var);

    if (!is.null(seed))
        set.seed(seed);

    ## prepare data
    data   <- data[!is.na(data[["_strata_"]]), ]
    S      <- max(data[["_strata_"]]);

    ## find mwle
    rst.theta <- NULL;
    for (i in 1:S) {
        cur.d1 <- data[data[["_strata_"]] == i & data[["_grp_"]] == 1, v.outcome];
        cur.d0 <- data[data[["_strata_"]] == i & data[["_grp_"]] == 0, v.outcome];

        ns1 <- length(cur.d1);
        ns0 <- length(cur.d0);

        if (0 == ns1) {
            stop(paste("Stratum ", i, " contains no subjects from group 1", sep = ""));
        }

        cur.lambda <- lambdas[i];
        cur.theta  <- rweWL(cur.data = cur.d1, ext.data = cur.d0, lambda = cur.lambda, ...);

        ##bootstrap or jackknife
        var.theta  <- NULL;
        if ("bs" == m.var) {
            for (j in 1:bs.n) {
                cur.d1.bs <- sample(cur.d1, replace = TRUE);
                cur.bs    <- rweWL(cur.data = cur.d1.bs, ext.data = cur.d0,
                                   lambda = cur.lambda, ...);
                var.theta <- c(var.theta, cur.bs);
            }

            var.mle <- var(var.theta);
        } else if ("jk" == m.var) {
            for (j in 1:ns1) {
                cur.bs <- rweWL(cur.data = cur.d1[-j], ext.data = cur.d0,
                                lambda = cur.lambda, ...);
                var.theta <- c(var.theta, cur.bs);
            }

            if (ns0 > 0) {
                for (j in 1:ns0) {
                    cur.bs    <- rweWL(cur.data = cur.d1, ext.data = cur.d0[-j],
                                       lambda = cur.lambda, ...);
                    var.theta <- c(var.theta, cur.bs);
                }
            }

            var.mle <- (ns1+ns0-1)/(ns1+ns0)*sum((var.theta - cur.theta)^2);
        } else {
            var.mle <- 0
        }

        rst.theta <- rbind(rst.theta, c(ns1, cur.theta, var.mle, ns0));
    }

    ##mwle
    ws       <- rst.theta[,1] / sum(rst.theta[,1]);
    rst.mwle <- sum(ws * rst.theta[,2]);
    ##rst.bs   <- apply(rst.theta[,-(1:2), drop = FALSE], 2, function(x) sum(ws*x));
    rst.var  <- sum(ws^2 * rst.theta[,3]);

    list(mwle        = rst.mwle,
         var         = rst.var,
         mwle.strata = rst.theta[, 2],
         var.strata  = rst.theta[, 3],
         ns1         = rst.theta[, 1],
         ns0         = rst.theta[, 4])
}



#' PS-Integrated Composite Likelihood Estimation
#'
#' Return all Jack-Knife values for all stratum. Works for
#' one--arm studies when there is one RWD
#'
#' @param data class DWITHPS data frame
#' @param ... parameters for \code{rweWL}
#' @param lambdas power parameter without standardization by ns0
#' @param seed random seed
#'
#' @export
#'
rwePsJkWl <- function(data, lambdas, v.outcome = "Y", ...) {

    stopifnot(v.outcome %in% colnames(data));

    ## prepare data
    data[["_id_"]] <- 1:nrow(data)

    data   <- data[!is.na(data[["_strata_"]]), ]
    S      <- max(data[["_strata_"]]);

    ## find mwle
    rst_theta <- NULL;
    for (i in 1:S) {
        cur.d1 <- data[data[["_strata_"]] == i &
                       data[["_grp_"]]    == 1, ]
        cur.d0 <- data[data[["_strata_"]] == i &
                       data[["_grp_"]]    == 0, ]

        ns1 <- nrow(cur.d1)
        ns0 <- nrow(cur.d0)

        if (0 == ns1) {
            stop(paste("Stratum ", i, " contains no subjects from group 1",
                       sep = ""))
        }

        ## overall estimate
        cur.lambda <- lambdas[i];
        cur.theta  <- rweWL(cur.data = cur.d1[, v.outcome],
                            ext.data = cur.d0[, v.outcome],
                            lambda = cur.lambda, ...)

        rst_theta  <- rbind(rst_theta,
                            c(i, NA, NA, ns1, ns0, cur.theta))

        for (j in 1:ns1) {
            cur.theta <- rweWL(cur.data = cur.d1[-j, v.outcome],
                               ext.data = cur.d0[, v.outcome],
                               lambda = cur.lambda, ...)

            rst_theta  <- rbind(rst_theta,
                                c(i, 1, cur.d1[j, "_id_"], ns1 - 1, ns0, cur.theta))
        }

        if (ns0 > 0) {
            for (j in 1:ns0) {
                cur.theta <- rweWL(cur.data = cur.d1[, v.outcome],
                                   ext.data = cur.d0[-j, v.outcome],
                                   lambda = cur.lambda, ...)
                rst_theta  <- rbind(rst_theta,
                                    c(i, 0, cur.d0[j, "_id_"], ns1, ns0-1, cur.theta))
            }
        }

    }
    colnames(rst_theta) <- c("Strata", "Group", "Inx", "N1", "N0", "Theta")
    data.frame(rst_theta)
}


#' PS-Integrated Weighted Likelihood Estimation
#'
#' For all stratum. Variance estimated by Jack-Knife. Works for
#' multiple external RWDs.
#'
#' @param data class D_GPS data frame
#' @param ... parameters for \code{rweWL_G}
#' @param lambdas matrix of number of patients to be borrowed
#'
#' @export
#'
rweGpsWl <- function(data, lambdas, v.outcome = "Y", ...) {

    stopifnot(v.outcome %in% colnames(data));
    ## prepare data
    data   <- data[!is.na(data[["_strata_"]]), ]
    S      <- max(data[["_strata_"]])
    nd     <- max(data[["_grp_"]])

    ## find mwle
    rst_theta <- NULL;
    for (i in 1:S) {
        cur_d <- list()
        ns_d  <- NULL
        for (j in 1:nd) {
            cur_d[[j]] <- data[data[["_strata_"]] == i & data[["_grp_"]] == j, v.outcome]
            cur_n      <- length(cur_d[[j]])
            if (0 == cur_n) {
                stop(paste("Stratum ", i, " contains no subjects from group ", "j", sep = ""))
            }
            ns_d <- c(ns_d, cur_n)
        }
        cur_lambda <- c(ns_d[1], lambdas[i, ])

        ## estimate
        est_theta  <- rweWL_G(cur_d, lambda = cur_lambda, ...)

        ## jackknife
        jk_theta   <- NULL;
        for (j in 1:nd) {
            for (k in 1:ns_d[j]) {
                cur_djk      <- cur_d
                cur_djk[[j]] <- cur_djk[[j]][-k]
                cur_theta    <- rweWL_G(cur_djk, lambda = cur_lambda, ...)
                jk_theta     <- c(jk_theta, cur_theta);
            }
        }
        var_theta <- (length(jk_theta)-1) / (length(jk_theta)) * sum((jk_theta - est_theta)^2)

        rst_theta <- rbind(rst_theta,
                           c(est_theta, var_theta, ns_d))
    }

    ##mwle
    ws       <- rst_theta[, 3] / sum(rst_theta[, 3])
    rst_mwle <- sum(ws * rst_theta[, 1])
    rst_var  <- sum(ws^2 * rst_theta[, 2])

    list(mwle        = rst_mwle,
         var         = rst_var,
         mwle.strata = rst_theta[, 1],
         var.strata  = rst_theta[, 2],
         ns1         = rst_theta[, 3],
         nsj         = rst_theta[, -(1:3)])
}



#' Weighted Likelihood Estimation for one external RWD
#'
#' @param cur.data data from current study
#' @param ext.data data from external study
#' @param type     type of outcomes
#' @param lambda   power parameter without standardization of ns0
#' @param equal.sd boolean. whether sd is the same between the current and external study
#'
#' @export
#'
rweWL <- function(cur.data, ext.data, lambda, type = c("continuous", "binary"), equal.sd = TRUE) {

    f.ll <- function(pars) {
        theta  <- pars[1];
        sig2.1 <- pars[2];
        sig2.0 <- pars[3];

        ll <- - n1 * log(sig2.1) / 2;
        ll <- ll - n1     * mean((cur.data - theta)^2)/2/sig2.1;
        ll <- ll - lambda * log(sig2.0) / 2;
        ll <- ll - lambda * mean((ext.data - theta)^2)/2/sig2.0;

        ll
    }

    f.gradient <- function(pars) {
        theta  <- pars[1];
        sig2.1 <- pars[2];
        sig2.0 <- pars[3];

        g <- numeric(length(pars));
        ## d logl / d theta
        g[1] <- n1/sig2.1*(mean(cur.data) - theta) + lambda/sig2.0*(mean(ext.data) - theta);
        ## d logl / d sig2.1
        g[2] <- - n1/2/sig2.1     + n1     * mean((cur.data - theta)^2)/2/sig2.1/sig2.1;
        ## d logl / d sig2.0
        g[3] <- - lambda/2/sig2.0 + lambda * mean((cur.data - theta)^2)/2/sig2.0/sig2.0;

        return(g)
    }


    type <- match.arg(type);
    n1   <- length(cur.data);

    ## ignore external data
    if (0 == lambda) {
        ext.data <- cur.data; ## placeholder
        equal.sd <- TRUE;
    }

    init.theta <- (n1/(n1 + lambda)) * mean(cur.data) + (lambda/(n1 + lambda)) * mean(ext.data);

    if (("continuous" == type & equal.sd) | "binary" == type) {
        rst <- init.theta;
    } else {
        init.sig2.1 <- mean((cur.data - init.theta)^2);
        init.sig2.0 <- mean((ext.data - init.theta)^2);
        rst         <- optim(c(init.theta, init.sig2.1, init.sig2.0),
                             method = "L-BFGS-B",
                             fn     = f.ll,
                             lower  = c(-Inf, 1e-6, 1e-6), upper = rep(Inf,3),
                             control=list(fnscale=-1))$par[1];
    }

    rst;
}


#' Weighted Likelihood Estimation for general number of external RWDs
#'
#' @param dta list of data from each source including current study
#' @param type type of outcomes
#' @param lambda number of patients to be included or borrowed. For current
#'     study, it should be equal to the number of patients
#'
#' @export
#'
rweWL_G <- function(lst_dta, lambda, type = c("continuous", "binary")) {
    type     <- match.arg(type);
    all_mean <- NULL
    for (i in seq_len(length(lst_dta))) {
        all_mean <- c(all_mean, mean(lst_dta[[i]]))
    }

    sum(all_mean * lambda / sum(lambda))
}


## #' Weighted Likelihood Estimation
## #'
## #' @param cur.data data from current study
## #' @param ext.data data from external study
## #' @param type     type of outcomes
## #' @param ext.ps   ps of the external study
## #' @param equal.sd boolean. whether sd is the same between the current and external study
## #'
## #' @export
## #'
## rweWL2 <- function(cur.data, ext.data, ext.ps, type = c("continuous", "binary"), equal.sd = FALSE) {

##     type       <- match.arg(type);
##     n1         <- length(cur.data);
##     init.theta <- mean(c(cur.data, ext.data));
##     init.sig   <- sd(c(cur.data, ext.data));

##     f.cont <- function(pars) {
##         theta  <- pars[1];
##         sig.1  <- exp(pars[2]);
##         sig.0  <- exp(pars[3]);

##         ll.1 <- dnorm(cur.data, theta, sd = sig.1, log = TRUE);
##         ll.0 <- dnorm(ext.data, theta, sd = sig.0, log = TRUE);

##         sum(ll.1) + sum(ll.0 * ext.ps);
##     }

##     if ("continuous" == type) {
##         rst <- optim(c(init.theta, log(init.sig), log(init.sig)),
##                      method = "Nelder-Mead",
##                      fn     = f.cont,
##                      control=list(fnscale=-1))$par[1];
##     }  else {
##         A   <- sum(cur.data) + sum(ext.ps * ext.data);
##         B   <- sum(1 - cur.data) + sum( ext.ps * (1-ext.data));
##         rst <- A/(A+B);
##     }

##     rst;
## }

