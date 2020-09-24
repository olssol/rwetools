#' Get unbalance in baseline Covariates
#'
#' @param diff If TRUE, get the difference in covariates between groups.
#'     Otherwise, get covariates for each group separately
#'
#' @param var.group Column name in pts that corresponds to treat group
#'
#' @inheritParams simupara
#'
#' @return If diff is TRUE, return a dataset with columns V and Diff. Otherwise,
#'     return a dataset with columns V, Z and Value. In both cases, column V
#'     represents the name of covariates.
#'
#' @export
#'
rweUnbalance <- function(nPat, ..., pts = NULL, covs = NULL, diff = TRUE,
                         var.group   = "Z",
                         cov.pattern = "^V[0-9]+$") {
    if (is.null(pts)) {
        pts <- rweSimuTwoArm(nPat, ...);
    }

    ##unbalance
    inx_0 <- which(0 == pts[, var.group]);
    if (is.null(covs)) {
        c.xy  <- colnames(pts);
        c.xy  <- c.xy[grep(cov.pattern, c.xy)];
    } else {
        c.xy <- covs;
    }

    unb  <- NULL;
    for (i in 1:length(c.xy)) {
        x0     <- sample(pts[inx_0,  c.xy[i]], size = nPat, replace = TRUE)
        x1     <- sample(pts[-inx_0, c.xy[i]], size = nPat, replace = TRUE)

        if (diff) {
            x.diff <- x1 - x0;
            unb    <- rbind(unb, data.frame(V    = c.xy[i],
                                            Diff = x1 - x0));
        } else {
            unb    <- rbind(unb, data.frame(V = c.xy[i],
                                            Z = 1,
                                            Value = x1));
            unb    <- rbind(unb, data.frame(V = c.xy[i],
                                            Z = 0,
                                            Value = x0));
        }
    }

    ## make group column factor if it exists
    if (!diff) {
        unb$Z <- as.factor(unb$Z)
    }

    unb
}

#' Get unbalance after PS adjustment
#'
#' @param pts data frame of datasets
#' @param var_group Column name in pts that corresponds to treat group
#' @param var_group Column name in pts that corresponds to propensity score
#'
#' @return
#'
#' @export
#'
rwePSUnbalance <- function(pts, covs = NULL, diff = TRUE,
                         var.group   = "Z",
                         cov.pattern = "^V[0-9]+$") {
    if (is.null(pts)) {
        pts <- rweSimuTwoArm(nPat, ...);
    }

    ##unbalance
    inx_0 <- which(0 == pts[, var.group]);
    if (is.null(covs)) {
        c.xy  <- colnames(pts);
        c.xy  <- c.xy[grep(cov.pattern, c.xy)];
    } else {
        c.xy <- covs;
    }

    unb  <- NULL;
    for (i in 1:length(c.xy)) {
        x0     <- sample(pts[inx_0,  c.xy[i]], size = nPat, replace = TRUE)
        x1     <- sample(pts[-inx_0, c.xy[i]], size = nPat, replace = TRUE)

        if (diff) {
            x.diff <- x1 - x0;
            unb    <- rbind(unb, data.frame(V    = c.xy[i],
                                            Diff = x1 - x0));
        } else {
            unb    <- rbind(unb, data.frame(V = c.xy[i],
                                            Z = 1,
                                            Value = x1));
            unb    <- rbind(unb, data.frame(V = c.xy[i],
                                            Z = 0,
                                            Value = x0));
        }
    }

    ## make group column factor if it exists
    if (!diff) {
        unb$Z <- as.factor(unb$Z)
    }

    unb
}


#' Get balance by different metrics
#'
#' @param covs covariate matrix
#' @param grp  group indicator
#' @param metric different metrics for measuring balance
#'
#' @export
#'
rwe_bal_metric <- function(covs, grp,
                           metric = c("std", "abd", "ovl", "ksd",
                                      "ley", "mhb")) {
    ## check parameters
    metric <- match.arg(metric, several.ok = TRUE)
    stopifnot(nrow(covs) == length(grp))

    ## get covariate names
    cnames <- colnames(covs)
    if (is.null(cnames))
        cnames <- paste("V", seq_len(ncol(covs)), sep = "")

    cov0 <- covs[0 == grp, ]
    cov1 <- covs[1 == grp, ]

    ## check balance
    rst <- NULL
    for (m in metric) {
        if (m == "mhb") {
            cur_rst <- get_metric(cov0, cov1, metric = m)
            cur_rst <- data.frame(metric = m,
                                  cov    = "",
                                  value  = cur_rst)
            rst     <- rbind(rst, cur_rst)
        } else {
            for (j in seq_len(ncol(cov0))) {
                cur_rst <- get_metric(cov0[, j], cov1[, j], metric = m)
                cur_rst <- data.frame(metric = m,
                                      cov    = cnames[j],
                                      value  = cur_rst)
                rst     <- rbind(rst, cur_rst)
            }
        }
    }

    data.frame(rst)
}


#' Get PS-adjusted balance by different metrics
#'
#' @inheritParams rwe_bal_metric
#'
#'
#' @export
#'
rwe_bal_metric_ps <- function(covs, grp, ps,
                              adjust = c("matching", "stratification"),
                              n_strata = 5,
                              m_ratio = 3,
                              ...) {
    ## check parameters
    adjust <- match.arg(adjust)
    stopifnot(nrow(covs) == length(grp))
    stopifnot(nrow(covs) == length(ps))

    inx_1 <- which(1 == grp)
    cov1  <- covs[inx_1, ]
    ps1   <- ps[inx_1]
    cov0  <- covs[-inx_1, ]
    ps0   <- ps[-inx_1]

    if ("matching" == adjust) {
        r_match <- get_match(target = ps1, candidate = ps0, ratio = m_ratio)
        r_covs  <- rbind(cov1, cov0[r_match[, 2], ])
        r_grp   <- c(rep(1, nrow(cov1)),
                     rep(0, nrow(r_match)))
        rst <- rwe_bal_metric(covs = r_covs, grp = r_grp, ...)
    } else {
        r_strata <- rweCut(ps1, c(ps1, ps0), breaks = n_strata)
        strata_1 <- r_strata[seq_len(length(ps1))]
        strata_0 <- r_strata[-seq_len(length(ps1))]

        rst <- NULL
        for (i in seq_len(n_strata)) {
            cur_cov1 <- cov1[which(i == strata_1), ]
            cur_cov0 <- cov0[which(i == strata_0), , drop = FALSE]

            cur_grp  <- c(rep(1, nrow(cur_cov1)),
                          rep(0, nrow(cur_cov0)))

            cur_rst  <- rwe_bal_metric(covs = rbind(cur_cov1, cur_cov0),
                                       grp  = cur_grp, ...)

            if (1 == i) {
                rst <- cur_rst
            } else {
                rst$value <- rst$value + cur_rst$value
            }
        }
    }

    rst
}


#' Get balance by different metrics
#'
#' @inheritParams rwe_bal_metric
#'
#' @export
#'
get_metric <- function(cov0, cov1,
                       metric = c("std", "abd", "ovl", "ksd",
                                  "ley", "mhb")) {
    metric <- match.arg(metric);
    switch(metric,
           std = {
               s <- sqrt((var(cov1) + var(cov0)) / 2)
               abs(mean(cov1) - mean(cov0)) / s;
           },
           abd = abs(mean(cov0) - mean(cov1)),
           ovl = metric_ovl(cov0, cov1),
           ksd = metric_ksd(cov0, cov1),
           ley = metric_ley(cov0, cov1),
           mhb = metric_mhb(cov0, cov1)
           )
}


## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
##             PRIVATE FUNCTIONS
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

## overlapping coefficient
metric_ovl <- function(cov0, cov1) {
  cov <- c(cov0, cov1);
  if (length(unique(cov)) <= 10) {
    all.x <- c(rep(0, length(cov0)), rep(1, length(cov1)));
    pt    <- apply(prop.table(table(cov, all.x), 2), 1, min);

    ## reversed to measure imbalance
    return(1 - sum(pt))
  }

  mn <- min(cov) * 1.25;
  mx <- max(cov) * 1.25;
  f1 <- approxfun(density(cov1, from = mn, to = mx,
                          bw = "nrd"));
  f0 <- approxfun(density(cov0, from = mn, to = mx,
                          bw = "nrd"));

  fn <- function(x)
    pmin(f1(x), f0(x))

  s <- try(integrate(fn, lower = mn, upper = mx,
                     subdivisions = 500)$value)
  ## Reverse: measure imbalance
  ifelse(inherits(s, "try-error"), NA, 1 - s)
}

## K-S distance
metric_ksd <- function(cov0, cov1) {
    cov    <- c(cov0, cov1);
    cdf_1  <- ecdf(cov1);
    cdf_0  <- ecdf(cov0);
    max(abs(cdf_1(cov) - cdf_0(cov)))
}

## Levy distance
metric_ley <- function(cov0, cov1) {
    cov   <- c(cov0, cov1);
    cdf_1 <- ecdf(cov1);
    cdf_0 <- ecdf(cov0);
    e     <- max(abs(cdf_1(cov) - cdf_0(cov)))

    if (length(unique(cov)) <= 10)
        return(e)

    x     <- seq(min(cov), max(cov), length.out = 1000);
    check <- all(cdf_0(x - e) - e <= cdf_1(x) &
                 cdf_1(x) <= cdf_0(x + e) + e)

    while (check) {
        e <- e - .01
        check <- all(cdf_0(x - e) - e <= cdf_1(x) &
                     cdf_1(x) <= cdf_0(x + e) + e)
    }

    e
}

#' mahalanobis balance
#'
#' covs should be a reduced datset that contains only those covariates
#' that will be used for calculating Mahalanobis balance, for example,
#' covs = dat[,1:6]
#' trt should be the exposure variable,
#' for example, trt=dat$X
#'
metric_mhb <- function(cov0, cov1) {
  cov01 <- rbind(cov0, cov1)
  sinv  <- solve(cov(cov01))
  x0    <- colMeans(cov0)
  x1    <- colMeans(cov1)

  sum((t(x1 - x0) %*% sinv) * (x1 - x0))
}
