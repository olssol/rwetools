#' Get propensity scores
#'
#' @param ... parameters to get propensity scores
#'
#' @export
#'
rwePS <- function(data, formula = NULL, v.grp = "group", v.covs = "V1", d1.grp = 1, nstrata = 5, ...) {

    dnames <- colnames(data);
    stopifnot(v.grp %in% dnames);

    all.ps  <- get.ps(data, ps.fml = formula, ps.cov = v.covs, grp = v.grp, ...);
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
    class(data)        <- append(class(data),
                                 get.rwe.class("DWITHPS"));
    data
}

#' Get number of subjects and the KL distance for each PS strata
#'
#'
#' @export
#'
rwePSKL <- function(data.withps, ...) {

    stopifnot(inherits(data.withps,
                       what = get.rwe.class("DWITHPS")));

    nstrata <- max(data.withps[["_strata_"]], na.rm = TRUE);
    rst     <- NULL;
    for (i in 1:nstrata) {
        ps0 <- data.withps[i == data.withps[["_strata_"]] &
                           0 == data.withps[["_grp_"]],
                           "_ps_"];
        ps1 <- data.withps[i == data.withps[["_strata_"]] &
                           1 == data.withps[["_grp_"]],
                           "_ps_"];

        cur.kl <- rweKL(ps0, ps1, ...);
        rst    <- rbind(rst, c(i, cur.kl));
    }

    colnames(rst) <- c("Strata", "N0", "N1", "KL");
    rst           <- data.frame(rst);
    class(rst)    <- append(class(rst), get.rwe.class("PSKL"));

    rst
}

#' Get the actual power term in the power prior
#'
#' @param a           power term
#' @param overall.ess ratio of overall added number of patients to N1
#' @param adjust.size whether adjust for sizes in group 0 and 1 in the power term
#' @param adjust.kl   whether adjust for KL distance in ps scores in the power term
#'
#' @export
#'
rweGetPowerA <- function(pskl, a = NULL, overall.ess = 0.3, adjust.size = TRUE, adjust.kl = TRUE) {

    stopifnot(inherits(pskl, what = get.rwe.class("PSKL")));

    ## compute a
    if (is.null(a)) {
        stopifnot(1 == adjust.size);
        stopifnot(overall.ess >= 0);

        if (1 == adjust.kl) {
            a <- overall.ess * (1 + mean(pskl$KL));
        } else {
            a <- overall.ess;
        }
    }


    ## compute as, power term for each strata
    rst <- rep(a, nrow(pskl));
    if (1 == adjust.size) {
        rst <- rst * pskl$N1 / pskl$N0;
    }

    if (1 == adjust.kl) {
        rst <- rst / (1 + pskl$KL)
    }

    list(a  = a,
         as = rst);
}
