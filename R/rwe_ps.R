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
        ps.fml <- as.formula(paste(v.grp, "~",
                                   paste(v.covs, collapse = "+"),
                                   sep = ""))

    ## d1 index will be kept in the results
    d1.inx   <- d1.grp == data[[v.grp]];
    keep.inx <- which(d1.inx);

    ## for 2-arm studies only
    if (!is.null(d1.arm))
        d1.inx <- d1.inx & d1.arm == data[[v.arm]]

    ## get ps
    all.ps  <- rwe_tk_ps(data, ps_fml = ps.fml, d1_grp = d1.grp,
                         ...)
    D1.ps   <- all.ps[which(d1.inx)]

    ## add columns to data
    grp     <- rep(1, nrow(data));
    grp[which(data[[v.grp]] != d1.grp)] <- 0;

    data[["_ps_"]]     <- all.ps;
    data[["_grp_"]]    <- grp;
    data[["_arm_"]]    <- data[[v.arm]];

    ## stratification
    if (nstrata > 0) {
        strata  <- rweCut(D1.ps, all.ps,
                          breaks = nstrata,
                          keep.inx = keep.inx)
        data[["_strata_"]] <- strata;
    }

    ## return
    rst <- list(data    = data,
                ps.fml  = ps.fml,
                nstrata = nstrata);

    class(rst) <- get.rwe.class("DWITHPS");

    rst
}

#' Get generalized propensity scores
#'
#' @param ... parameters to get propensity scores
#'
#' @export
#'
rweGPS <- function(data, ps.fml = NULL,
                    v.grp   = "group",
                    v.covs  = "V1",
                    nstrata = 5, ...) {

    ## likelihood
    f_ll <- function(beta) {
        xbeta <- apply(d_mat, 1,
                       function(x) sum(x * beta[-(1:nd1)]))

        exb <- sapply(xbeta,
                      function(x) {
                          tx <- beta[1:nd1] + x
                          tx <- 1 + sum(exp(tx))
                          log(tx)
                      })

        rst <- sum(beta[1:nd1] * n_j)
        rst <- rst + sum(xbeta[inx_j])
        rst <- rst - sum(exb)
        rst
    }

    f_gradient <- function(beta) {
        xbeta <- apply(d_mat, 1,
                       function(x) sum(x * beta[-(1:nd1)]))

        exb <- sapply(xbeta,
                      function(x) {
                          tx <- beta[1:nd1] + x
                          tx <- exp(tx)
                      })
        s_exb  <- apply(exb, 2, sum)
        s_1exb <- sapply(s_exb, function(x) {1 / (1+x)})

        g <- numeric(length(beta))
        for (i in 1:nd1) {
            g[i] <- n_j[i] - sum(s_1exb * exb[i,])
        }

        for (i in 1:nx) {
            g[i + nd1] <- sum(d_mat[inx_j, i]) - sum(s_1exb * s_exb * d_mat[, i])
        }

        g
    }

    dnames <- colnames(data);
    stopifnot(v.grp %in% dnames);

    nd  <- max(data[[v.grp]])
    nd1 <- nd - 1
    n_j <- NULL
    for (i in 2:nd) {
        n_j <- c(n_j, sum(data[[v.grp]] == i))
    }
    inx_j <- which(data[[v.grp]] > 1)

    ## design matrix
    if (is.null(ps.fml))
        ps.fml <- as.formula(paste(v.grp, "~", paste(v.covs, collapse = "+"), sep = ""))

    d_mat <- model.matrix(ps.fml, data)[, -1]

    ## mle
    ## mle0 <- optim(rep(0, nd1 + NC), f_ll,
    ##               method = "Nelder-Mead",
    ##               control = list(fnscale = -1, ...))
    mle <- optim(rep(0, nd1 + ncol(d_mat)),
                 f_ll, gr = f_gradient,
                 method = "BFGS", control = list(fnscale = -1, ...))

    ps_par <- mle$par

    ## gps and strata
    gps    <- apply(d_mat, 1,
                    function(x) sum(x * mle$par[-(1:nd1)]))

    strata <- rweCut(gps[which(1 == data[[v.grp]])], gps, breaks = nstrata)

    ## return
    data[["_gps_"]]    <- gps
    data[["_strata_"]] <- strata
    data[["_grp_"]]    <- data[[v.grp]]

    rst <- list(data    = data,
                ps.fml  = ps.fml,
                nd      = nd,
                nstrata = nstrata,
                mle     = mle)

    class(rst) <- get.rwe.class("D_GPS");
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
rwePSDist <- function(data.withps, n.bins = 10, min.n0 = 10,
                      type = c("ovl", "kl"), d1.arm = NULL, ...) {

    f.narm <- function(inx, dataps) {

        if (is.null(dataps[["_arm_"]]))
            return(c(length(inx), 0,0));

        n0 <- length(which(0 == dataps[inx, "_arm_"]));
        n1 <- length(which(1 == dataps[inx, "_arm_"]));

        c(length(inx), n0 ,n1);
    }

    stopifnot(inherits(data.withps,
                       what = get.rwe.class("DWITHPS")));

    type     <- match.arg(type)
    dataps   <- data.withps$data
    nstrata  <- data.withps$nstrata
    rst      <- NULL;
    for (i in 1:nstrata) {

        inx.ps0 <- i == dataps[["_strata_"]] & 0 == dataps[["_grp_"]];
        inx.ps1 <- i == dataps[["_strata_"]] & 1 == dataps[["_grp_"]];
        n0.01   <- f.narm(which(inx.ps0), dataps);
        n1.01   <- f.narm(which(inx.ps1), dataps);

        if (!is.null(d1.arm) & !is.null(dataps[["_arm_"]])) {
            inx.ps0 <- inx.ps0 & d1.arm == dataps[["_arm_"]];
            inx.ps1 <- inx.ps1 & d1.arm == dataps[["_arm_"]];
        }

        ps0 <- dataps[which(inx.ps0), "_ps_"]
        ps1 <- dataps[which(inx.ps1), "_ps_"]

        if (0 == length(ps0) | 0 == length(ps1))
            warning("No samples in strata");

        if (any(is.na(c(ps0, ps1))))
            warning("NA found in propensity scores in a strata");

        if (length(ps0) < min.n0) {
            warning("Not enough data in the external data in the current stratum.
                     External data ignored.");
            cur.dist <- 0;
        } else {
            cur.dist <- rweDist(ps0, ps1, n.bins = n.bins, type = type, ...);
        }

        rst <- rbind(rst, c(i, n0.01, n1.01, cur.dist));
    }

    ## overall
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

#' Get number of subjects and the distances of PS distributions for each PS
#' strata for multiple data sources
#'
#' @param min.n0 threshold for N0, below which the external data in the
#'     current stratum will be ignored by setting the PS distance to 0
#'
#' @param d1.arm calculate distance based on a specific arm. Ignored if NULL.
#'
#' @export
#'
rweGpsDist <- function(data.gps, n.bins = 10, min.n0 = 10,
                       type = c("ovl", "kl"),  ...) {

    stopifnot(inherits(data.gps,
                       what = get.rwe.class("D_GPS")));

    type     <- match.arg(type)
    dataps   <- data.gps$data
    nstrata  <- data.gps$nstrata
    nd       <- data.gps$nd

    rst      <- NULL
    for (i in 1:nstrata) {
        inx.ps1 <- i == dataps[["_strata_"]] & 1 == dataps[["_grp_"]]
        ps1     <- dataps[which(inx.ps1), "_gps_"];
        if (0 == length(ps1))
            warning(paste("No samples in strata", i, "in current study"))

        dist <- length(ps1)
        for (j in 2:nd) {
            inx.psj <- i == dataps[["_strata_"]] & j == dataps[["_grp_"]]
            psj     <- dataps[which(inx.psj), "_gps_"]

            if (0 == length(psj))
                warning(paste("No samples in strata", i, "in Study", j))

            if (length(psj) < min.n0) {
                warning("Not enough data in the external data in
                         the current stratum.")
                cur_dist <- 0
            } else {
                cur_dist <- rweDist(psj, ps1, n.bins = n.bins, type = type, ...)
            }

            dist <- c(dist, length(psj), cur_dist)
        }
        rst <- rbind(rst, c(i, dist))
    }

    ##overall
    inx.tot.ps1 <- which(1 == dataps[["_grp_"]]);
    ps1         <- dataps[inx.tot.ps1, "_gps_"];

    dist <- length(ps1)
    for (j in 2:nd) {
        inx.tot.psj <- which(j == dataps[["_grp_"]])
        psj         <- dataps[inx.tot.psj, "_gps_"]
        cur_dist    <- rweDist(psj, ps1, n.bins = nstrata * n.bins,
                               type = type, ...)
        dist        <- c(dist, length(psj), cur_dist)
    }
    rst <- rbind(rst, c(0, dist));


    colnames(rst) <- c("Strata", "N1",
                       paste(rep(c("N", "Dist"), nd - 1),
                             rep(2:nd, each = 2), sep = ""))
    rst           <- data.frame(rst);
    class(rst)    <- append(get.rwe.class("GPSDIST"), class(rst));

    rst
}


#' Get the actual power term in the power prior
#'
#' @param psdist A RWE_PSDIST type object
#' @param a power term
#' @param overall.ess ratio of overall added number of patients to N1
#' @param adjust.size whether adjust for sizes in group 0 and 1 in the power
#'     term
#' @param adjust.dist whether adjust for distance in ps scores in the power term
#'
#' @export
#'
rweGetPowerA <- function(psdist, a = NULL, overall.ess = 0.3,
                         adjust.size = TRUE, adjust.dist = TRUE) {

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

#' Match on PS
#'
#' Match patients in RWD with patients in current study based on PS
#'
#' @param dta_cur current study data
#' @param dta_ext external data source data
#'
#' @export
#'
rwe_match_ps <- function(dta_cur, dta_ext, n_match = 3, ps.fml = NULL,
                         v.covs  = "V1") {

    n_cur <- nrow(dta_cur)
    n_ext <- nrow(dta_ext)
    ratio <- min(n_match, floor(n_ext / n_cur))

    dta_cur$grp_tmp <- 1
    dta_ext$grp_tmp <- 0
    dta_ext         <- dta_ext[, colnames(dta_cur)]

    dta    <- rbind(dta_cur, dta_ext)
    dta_ps <- rwePS(data = dta, v.grp = "grp_tmp", v.covs = v.covs, nstrata = 0)

    ps        <- dta_ps$data[["_ps_"]]
    target    <- ps[1:n_cur]
    candidate <- ps[-(1:n_cur)]

    ## match
    rst <- cMatch(target, candidate, ratio = ratio)[]
    rst <- rst[1:(ratio * n_cur)] + 1

    cbind(pid      = rep(1:n_cur, each = ratio),
          match_id = rst)
}

#' Get Propensity Scores
#'
#' General function for estimating PS
#'
#'
#' @export
#'
rwe_tk_ps <- function(dta, ps_fml = NULL, d1_grp = 1,
                       type = c("logistic", "randomforest", "gmm"),
                       ...,
                       grp = NULL, ps_cov = NULL,
                       ntree = 5000) {

    type <- match.arg(type)

    ## generate formula
    if (is.null(ps_fml))
        ps_fml <- as.formula(paste(grp, "~",
                                   paste(ps_cov, collapse = "+"),
                                   sep = ""))

    ## identify grp if passed from formula
    grp <- all.vars(ps_fml)[1];

    ## fit model
    switch(type,
           logistic = {
               ## d1_grp vs. all other groups
               dta[[grp]] <- d1_grp == dta[[grp]]
               glm_fit    <- glm(ps_fml, family = "binomial", data = dta, ...)
               est_ps     <- glm_fit$fitted
           },

           randomforest = {
               dta[[grp]] <- as.factor(dta[[grp]])
               rf_fit     <- randomForest(ps_fml, data = dta,
                                          ntree = ntree, ...);
               est_ps     <- predict(rf_fit, type = "prob")[, 2]
           },

           gmm = {
               gmm_fit <- rwe_gmm_ps(grp, ps_fml,
                                     dta = dta,
                                     d1_grp = d1_grp,
                                     ...)
               est_ps  <- gmm_fit$fitted
           })

    est_ps
}

#' Get Propensity Scores by GMM
#'
#' Estimating PS by GMM for three groups
#'
#'
#' @export
#'
rwe_gmm_ps <- function(grp, ps_fml, dta, att = 0, d1_grp = 1,
                       method  = "BFGS",
                       crit    = 1e-10,
                       itermax = 5000,
                       control = list(reltol = 1e-10, maxit  = 20000),
                       ...) {

    f_mm <- function(beta, mat_grp_x) {
        rst <- c_ps_gmm_g(beta, mat_grp_x, att = att)
        rst
    }

    f_mm_d <- function(beta, mat_grp_x) {
        rst <- c_ps_gmm_dg(beta, mat_grp_x, att = att)
        rst
    }

    ## check and convert groups
    d_grp <- dta[[grp]]
    lvl_g <- sort(unique(d_grp))
    stopifnot(3 == length(lvl_g))

    lvl_g   <- c(d1_grp,
                 lvl_g[-which(d1_grp == lvl_g)])

    lst_inx <- list()
    for (i in seq_len(length(lvl_g))) {
        lst_inx[[i]] <- which(lvl_g[i] == d_grp)
    }

    for (i in seq_len(length(lst_inx))) {
        d_grp[lst_inx[[i]]] <- i
    }

    ## x and grp
    x         <- model.matrix(ps_fml, dta)
    mat_grp_x <- cbind(d_grp, x)

    ## get initial values
    d_grp_01 <- d_grp
    d_grp_01[which(1 != d_grp_01)] <- 0
    dta[[grp]] <- d_grp_01
    glm_fit    <- glm(ps_fml, family = "binomial", data = dta)

    ## fit gmm
    gmm_fit <- gmm::gmm(f_mm,
                        mat_grp_x,
                        gradv   = NULL, ## f_mm_d,
                        t0      = coefficients(glm_fit),
                        method  = method,
                        crit    = crit,
                        itermax = itermax,
                        control = control,
                        ...)

    ## g values
    ## print(apply(gmm_fit$gt, 2, mean))

    gmm_beta <- coefficients(gmm_fit)
    xbeta    <- x %*% gmm_beta
    exbeta   <- exp(xbeta)
    ps       <- exbeta / (1 + exbeta)

    list(fitted = ps,
         beta   = gmm_beta)
}

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
##             PRIVATE FUNCTIONS
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
