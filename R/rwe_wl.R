#' PS-Integrated Weighted Likelihood Estimation for all stratum
#'
#' @param data class DWITHPS data frame
#' @param A target number of subjects to be borrowed
#' @param RS parameters for dirichelet prior
#' @param ... parameters for \code{rweWL}
#' @param bs.n number of bootstraps
#' @param bs.conf bootstrap confidence interval level
#' @param seed random seed
#'
#' @export
#'
rwePsWL <- function(data, RS = NULL, A = 0, v.outcome = "Y", bs.n = 1000, bs.conf = c(0.025, 0.975),
                    seed = NULL, ...) {

    stopifnot(v.outcome %in% colnames(data));

    if (!is.null(seed))
        set.seed(seed);

    ## prepare data
    data   <- data[!is.na(data[["_strata_"]]),];
    S      <- max(data[["_strata_"]]);

    ## distance rs
    if (is.null(RS))
        RS <- rep(1, S);

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

        cur.lambda <- min(ns0, A * RS[i]/sum(RS));
        cur.theta  <- rweWL(cur.data = cur.d1, ext.data = cur.d0, lambda = cur.lambda, ...);

        bs.theta <- NULL;
        for (j in 1:bs.n) {
            cur.d1.bs <- sample(cur.d1, replace = TRUE);
            cur.bs    <- rweWL(cur.data = cur.d1.bs, ext.data = cur.d0, lambda = cur.lambda, ...);
            bs.theta  <- c(bs.theta, cur.bs);
        }

        rst.theta  <- rbind(rst.theta, c(ns1, cur.theta, bs.theta));
    }

    ##mwle
    ws       <- rst.theta[,1]/sum(rst.theta[,1]);
    rst.mwle <- sum(ws * rst.theta[,2]);
    rst.bs   <- apply(rst.theta[,-(1:2), drop = FALSE], 2, function(x) sum(ws*x));

    list(mwle   = rst.mwle,
         bs.ci  = quantile(rst.bs, probs = bs.conf),
         bs.var = var(rst.bs),
         mwle.strata = rst.theta[, 2]);
}


#' Weighted Likelihood Estimation
#'
#' @param cur.data data from current study
#' @param ext.data data from external study
#' @param type     type of outcomes
#' @param lambda    power parameter
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

    f.gradient <- function(par) {
        theta  <- pars[1];
        sig2.1 <- pars[2];
        sig2.0 <- pars[3];

        g <- numeric(length(pars));
        ## d logl / d theta
        g[1] <- n1/sig2.1*(mean(cur.data) - theta) + lambda/sig2.0*(mean(ext.data) - theta);
        ## d logl / d sig2.1
        g[2] <- - n1/2/sig2.1     + n1     * mean((cur.data - theta)^2)/2/sig2.1/sig2.1;
        ## d logl / d sig2.0
        g[2] <- - lambda/2/sig2.0 + lambda * mean((cur.data - theta)^2)/2/sig2.0/sig2.0;

        return(g)
    }


    type <- match.arg(type);
    n1   <- length(cur.data);

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
