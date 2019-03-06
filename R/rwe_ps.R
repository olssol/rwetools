#' Get propensity scores
#'
#' @param ... parameters to get propensity scores
#' @param d1.arm define strata based on a specific arm. Ignored if NULL.
#'
#' @export
#'
rwePS <- function(data, ps.fml = NULL,
                  v.grp   = "group",
                  v.arm   = "arm",
                  v.covs  = "V1",
                  d1.arm  = NULL,
                  d1.grp  = 1,
                  nstrata = 5, ...) {

    dnames <- colnames(data);
    stopifnot(v.grp %in% dnames);

    ## generate formula
    if (is.null(ps.fml))
        ps.fml <- as.formula(paste(v.grp, "~", paste(v.covs, collapse="+"), sep=""));


    ## d1 index
    d1.inx  <- d1.grp == data[[v.grp]];
    if (!is.null(d1.arm))
        d1.inx <- d1.inx & d1.arm == data[[v.arm]];

    all.ps  <- get.ps(data, ps.fml = ps.fml, ...);
    D1.ps   <- all.ps[which(d1.inx)];

    ## stratification
    strata  <- rweCut(D1.ps, all.ps, breaks = nstrata);
    grp     <- rep(1, nrow(data));
    grp[which(data[[v.grp]] != d1.grp)] <- 0;

    data[["_ps_"]]     <- all.ps;
    data[["_strata_"]] <- strata;
    data[["_grp_"]]    <- grp;
    data[["_arm_"]]    <- data[[v.arm]];

    rst <- list(data    = data,
                ps.fml  = ps.fml,
                nstrata = nstrata);

    class(rst) <- get.rwe.class("DWITHPS");

    rst
}

#' Get number of subjects and the distances of PS distributions for each PS
#' strata
#'
#' @param min.n0 threshold for N0, below which the external data in the
#'     current stratum will be ignored by setting the PS distance to 0
#'
#' @param d1.arm calculate distance based on a specific arm. Ignored if NULL.
#'
#' @export
#'
rwePSDist <- function(data.withps, n.bins = 10, min.n0 = 10, type = c("ovl", "kl"), d1.arm = NULL, ...) {

    f.narm <- function(inx, dataps) {

        if (is.null(dataps[["_arm_"]]))
            return(c(length(inx), 0,0));

        n0 <- length(which(0 == dataps[inx, "_arm_"]));
        n1 <- length(which(1 == dataps[inx, "_arm_"]));

        c(length(inx), n0 ,n1);
    }

    stopifnot(inherits(data.withps,
                       what = get.rwe.class("DWITHPS")));

    type <- match.arg(type);

    dataps   <- data.withps$data;
    nstrata  <- data.withps$nstrata;
    rst     <- NULL;
    for (i in 1:nstrata) {

        inx.ps0 <- i == dataps[["_strata_"]] & 0 == dataps[["_grp_"]];
        inx.ps1 <- i == dataps[["_strata_"]] & 1 == dataps[["_grp_"]];
        n0.01   <- f.narm(which(inx.ps0), dataps);
        n1.01   <- f.narm(which(inx.ps1), dataps);

        if (!is.null(d1.arm) & !is.null(dataps[["_arm_"]])) {
            inx.ps0 <- inx.ps0 & d1.arm == dataps[["_arm_"]];
            inx.ps1 <- inx.ps1 & d1.arm == dataps[["_arm_"]];
        }

        ps0 <- dataps[which(inx.ps0), "_ps_"];
        ps1 <- dataps[which(inx.ps1), "_ps_"];

        if (0 == length(ps0) | 0 == length(ps1))
            warning("No samples in strata");

        if (any(is.na(c(ps0, ps1))))
            warning("NA found in propensity scores in a strata");

        if (length(ps0) < min.n0) {
            warning("Not enough data in the external data in the current stratum. External data ignored.");
            cur.dist <- 0;
        } else {
            cur.dist <- rweDist(ps0, ps1, n.bins = n.bins, type = type, ...);
        }

        rst <- rbind(rst, c(i, n0.01, n1.01, cur.dist));
    }

    ##overall
    inx.tot.ps0 <- which(0 == dataps[["_grp_"]]);
    inx.tot.ps1 <- which(1 == dataps[["_grp_"]]);
    n0.tot.01   <- f.narm(inx.tot.ps0, dataps);
    n1.tot.01   <- f.narm(inx.tot.ps1, dataps);

    ps0         <- dataps[inx.tot.ps0, "_ps_"];
    ps1         <- dataps[inx.tot.ps1, "_ps_"];
    all.dist    <- rweDist(ps0, ps1, n.bins = nstrata*n.bins, type = type, ...);
    rst         <- rbind(rst, c(0, n0.tot.01, n1.tot.01, all.dist));


    colnames(rst) <- c("Strata", "N0", "N00", "N01", "N1", "N10", "N11", "Dist");
    rst           <- data.frame(rst);
    class(rst)    <- append(get.rwe.class("PSDIST"), class(rst));

    rst
}

#' Get the actual power term in the power prior
#'
#' @param psdist      A RWE_PSDIST type object
#' @param a           power term
#' @param overall.ess ratio of overall added number of patients to N1
#' @param adjust.size whether adjust for sizes in group 0 and 1 in the power term
#' @param adjust.dist whether adjust for distance in ps scores in the power term
#'
#' @export
#'
rweGetPowerA <- function(psdist, a = NULL, overall.ess = 0.3, adjust.size = TRUE, adjust.dist = TRUE) {

    stopifnot(inherits(psdist, what = get.rwe.class("PSDIST")));

    ## compute a
    if (is.null(a)) {
        stopifnot(1 == adjust.size);
        stopifnot(overall.ess >= 0);

        if (1 == adjust.dist) {
            a <- overall.ess / mean(psdist$Dist);
        } else {
            a <- overall.ess;
        }
    }


    ## compute as, power term for each strata
    rst <- rep(a, nrow(psdist));
    if (1 == adjust.size) {
        rst <- rst * psdist$N1 / psdist$N0;
    }

    if (1 == adjust.dist) {
        rst <- rst * psdist$Dist;
    }

    list(a  = a,
         as = rst);
}
