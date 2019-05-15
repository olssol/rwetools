#' PS-Integrated Weighted Likelihood Estimation for all stratum by bootstrap
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
rwePsWL <- function(data, lambdas, v.outcome = "Y", m.var = c("jk", "bs"),
                    bs.n = 1000, seed = NULL, ...) {

    stopifnot(v.outcome %in% colnames(data));

    m.var <- match.arg(m.var);

    if (!is.null(seed))
        set.seed(seed);

    ## prepare data
    data   <- data[!is.na(data[["_strata_"]]),];
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
        }

        rst.theta <- rbind(rst.theta, c(ns1, cur.theta, var.mle, ns0));
    }

    ##mwle
    ws       <- rst.theta[,1]/sum(rst.theta[,1]);
    rst.mwle <- sum(ws * rst.theta[,2]);
    ##rst.bs   <- apply(rst.theta[,-(1:2), drop = FALSE], 2, function(x) sum(ws*x));
    rst.var  <- sum(ws^2 * rst.theta[,3]);

    list(mwle        = rst.mwle,
         var         = rst.var,
         mwle.strata = rst.theta[,2],
         var.strata  = rst.theta[,3],
         ns1         = rst.theta[,1],
         ns0         = rst.theta[,4]);
}


#' Weighted Likelihood Estimation
#'
#' @param cur.data data from current study
#' @param ext.data data from external study
#' @param type     type of outcomes
#' @param lambda   power parameter
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



#' Weighted Likelihood Estimation
#'
#' @param cur.data data from current study
#' @param ext.data data from external study
#' @param type     type of outcomes
#' @param ext.ps   ps of the external study
#' @param equal.sd boolean. whether sd is the same between the current and external study
#'
#' @export
#'
rweWL2 <- function(cur.data, ext.data, ext.ps, type = c("continuous", "binary"), equal.sd = FALSE) {

    type       <- match.arg(type);
    n1         <- length(cur.data);
    init.theta <- mean(c(cur.data, ext.data));
    init.sig   <- sd(c(cur.data, ext.data));

    f.cont <- function(pars) {
        theta  <- pars[1];
        sig.1  <- exp(pars[2]);
        sig.0  <- exp(pars[3]);

        ll.1 <- dnorm(cur.data, theta, sd = sig.1, log = TRUE);
        ll.0 <- dnorm(ext.data, theta, sd = sig.0, log = TRUE);

        sum(ll.1) + sum(ll.0 * ext.ps);
    }

    if ("continuous" == type) {
        rst <- optim(c(init.theta, log(init.sig), log(init.sig)),
                     method = "Nelder-Mead",
                     fn     = f.cont,
                     control=list(fnscale=-1))$par[1];
    }  else {
        A   <- sum(cur.data) + sum(ext.ps * ext.data);
        B   <- sum(1 - cur.data) + sum( ext.ps * (1-ext.data));
        rst <- A/(A+B);
    }

    rst;
}

