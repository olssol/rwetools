#' Get propensity scores
#'
#'
#' @export
#'
rwePS <- function(data, formula = NULL, v.grp = "group", v.covs = "V1", d1.grp = 1, delta = 0, nstrata = 5) {

    dnames <- colnames(data);
    stopifnot(v.grp %in% dnames);
    all.ps <- get.ps(data, ps.fml = formula, ps.cov = v.covs, grp = v.grp, delta = delta);

    D1.ps   <- all.ps[which(d1.grp == data[[v.grp]])];
    cuts    <- quantile(D1.ps, seq(0, 1,length=nstrata+1));
    cuts[1] <- cuts[1] - 0.001;

    strata  <- rep(NA, nrow(data));
    for (i in 2:length(cuts)) {
        inx         <- which(all.ps > cuts[i-1] & all.ps <= cuts[i]);
        strata[inx] <- i-1;
    }

    grp <- rep(1, nrow(data));
    grp[which(data[[v.grp]] != d1.grp)] <- 0;

    data[["_ps_"]]     <- all.ps;
    data[["_strata_"]] <- strata;
    data[["_grp_"]]    <- grp;
    class(data)        <- c(class(data), "RWE_DWITHPS");
    data
}

