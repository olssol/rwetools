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
               s <- sqrt((var(cov1) + var(cov0))/2);
               abs(mean(cov1) - mean(cov0))/s;
           },
           abd = abs(mean(cov0) - mean(cov1)),
           ovl = metric.ovl(cov0, cov1),
           ksd = metric.ksd(cov0, cov1),
           ley = metric.ley(cov0, cov1),
           mhb = metric.mhb(cov0, cov1)
           )
}

