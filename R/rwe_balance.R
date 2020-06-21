#' Get unbalance in baseline covariates
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
                         var.group = "Z", cov.pattern = "^V[0-9]+$") {
    if (is.null(pts)) {
        pts <- rweSimuTwoArm(nPat, ...);
    }

    ##unbalance
    inx.0 <- which(0 == pts[, var.group]);
    if (is.null(covs)) {
        c.xy  <- colnames(pts);
        c.xy  <- c.xy[grep(cov.pattern, c.xy)];
    } else {
        c.xy <- covs;
    }

    unb   <- NULL;
    for (i in 1:length(c.xy)) {
        x0     <- sample(pts[inx.0,  c.xy[i]], size = nPat, replace = TRUE);
        x1     <- sample(pts[-inx.0, c.xy[i]], size = nPat, replace = TRUE);

        if (diff) {
            x.diff <- x1 - x0;
            unb    <- rbind(unb, data.frame(V=c.xy[i], Diff=x1 - x0));
        } else {
            unb    <- rbind(unb, data.frame(V=c.xy[i], Z=1, Value=x1));
            unb    <- rbind(unb, data.frame(V=c.xy[i], Z=0, Value=x0));
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
#' @param cov0 covariates from group 0
#' @param cov1 covariates from group 1
#'
#' @export
#'
rweBalMetric <- function(cov0, cov1, metric = c("std", "abd", "ovl", "ksd",
                                                "ley", "mhb")) {
    metric <- match.arg(metric);
    switch(metric,
           std = {
               s <- sqrt((var(cov1) + var(cov0)) / 2)
               abs(mean(cov1) - mean(cov0)) / s;
           },
           abd = abs(mean(cov0) - mean(cov1)),
           ovl = metric.ovl(cov0, cov1),
           ksd = metric.ksd(cov0, cov1),
           ley = metric.ley(cov0, cov1),
           mhb = metric.mhb(cov0, cov1)
           )
}

## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------
##             PRIVATE FUNCTIONS
## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------

## overlapping coefficient
metric.ovl <- function(cov0, cov1) {
  cov <- c(cov0, cov1);
  if (length(unique(cov)) <= 10) {
    all.x <- c(rep(0, length(cov0)), rep(1, length(cov1)));
    pt    <- apply(prop.table(table(cov, all.x), 2), 1, min);

    ## reversed to measure imbalance
    return(1-sum(pt))
  }

  mn <- min(cov) * 1.25;
  mx <- max(cov) * 1.25;
  f1 <- approxfun(density(cov1, from = mn, to = mx, bw = "nrd"));
  f0 <- approxfun(density(cov0, from = mn, to = mx, bw = "nrd"));

  fn <- function(x)
    pmin(f1(x), f0(x))

  s <- try(integrate(fn, lower = mn, upper = mx,
                     subdivisions = 500)$value)
  ## Reverse: measure imbalance
  ifelse(inherits(s, "try-error"), NA, 1-s)
}

## K-S distance
metric.ksd <- function(cov0, cov1) {
    cov <- c(cov0, cov1);
    F1  <- ecdf(cov1);
    F0  <- ecdf(cov0);
    max(abs(F1(cov) - F0(cov)));
}

## Levy distance
metric.ley <- function(cov0, cov1) {
    cov <- c(cov0, cov1);
    F1  <- ecdf(cov1);
    F0  <- ecdf(cov0);
    e   <- max(abs(F1(cov) - F0(cov)));

    if (length(unique(cov)) <= 10)
        return(e)

    x     <- seq(min(cov), max(cov), length.out=1000);
    check <- all(F0(x-e) - e <= F1(x) & F1(x) <= F0(x+e) + e)

    while (check) {
        e <- e-.01
        check <- all(F0(x-e) - e <= F1(x) & F1(x) <= F0(x+e) + e);
    }

    e
}

#' mahalanobis balance
#'
#' covs should be a reduced datset that contains only those covariates
#' that will be used for calculating Mahalanobis balance, for example,
#' covs = dat[,1:6]
#' trt should be the exposure variable, for example, trt=dat$X
metric.mhb <- function(cov0, cov1) {
  S    <- rbind(cov0, cov1)
  Sinv <- solve(cov(S))
  x0   <- colMeans(S[1:length(cov0), ])
  x1   <- colMeans(S[-(1:length(cov0)), ])

  sum((t(x1 - x0) %*% Sinv) * (x1 - x0))
}
