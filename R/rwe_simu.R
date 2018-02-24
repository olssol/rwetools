#' Simulate covariates
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuCov <- function(nPat, muCov, sdCov, corCov) {
    cov.x <- rmvnorm(nPat,
                     mean=muCov,
                     sigma=get.covmat(sdCov, corCov));

    colnames(cov.x) <- paste("V", 1:ncol(cov.x), sep="");
    data.frame(cov.x);
}

#' Simulate X*Beta
#'
#' @inheritParams simupara
#'
#' @export
#'
rweXBeta <- function(nPat, muCov, sdCov, corCov, regCoeff, cov.x = NULL, fmla = NULL) {
    stopifnot(inherits(fmla,"formula") | is.null(fmla));

    if (is.null(cov.x))
        cov.x <- rweSimuCov(nPat, muCov, sdCov, corCov);

    if (is.null(fmla)) {
        ## add intercept
        d.matrix <- cov.x;
    } else {
        d.matrix <- model.matrix(fmla, cov.x);
    }

    xbeta <- get.xbeta(d.matrix, regCoeff);
    xbeta
}


#' Compute standard error of the random error
#'
#' @inheritParams simupara
#'
#' @return standard error or the random error term
#'
#' @export
#'
rweGetYSig <- function(sig2Ratio, nPat=10000, xbeta = NULL, ...) {
    if (is.null(xbeta))
        xbeta   <- rweXBeta(nPat=nPat, ...);
    v.xbeta <- var(xbeta);
    ysig    <- sqrt(v.xbeta * sig2Ratio);
    ysig
}

#' Get intercept for a binary outcome.
#'
#' The binary outcome may be an outcome or a treatment assignment.
#'
#' @inheritParams simupara
#'
#' @return standard error or the random error term
#'
#' @export
#'
rweGetBinInt <- function(bin.mu, nPat=10000, xbeta = NULL, ...) {
    if (is.null(xbeta))
        ey <- rweXBeta(nPat, ...);

    fx <- function(b0) {
        expy <- exp(b0+ey);
        m    <- mean(expy/(1+expy));
        abs(m - bin.mu);
    }

    rst <- optimize(fx, c(-10+min(ey),10+max(ey)))$minimum;
    rst
}

#' Simulate random errors
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuError <- function(nPat, ysig = 1,
                         type = c("normal", "skewed", "mixture"),
                         skew.n = NULL, skew.p = NULL, ...) {
    type <- match.arg(type);
    rst <- switch(type,
                  normal = {rnorm(nPat, 0, ysig)},
                  skewed = {
        mu        <- skew.n * (1-skew.p) / skew.p;
        va        <- skew.n * (1-skew.p) / skew.p^2;
        noise.sig <- 0.1;
        rst       <- rnbinom(nPat, skew.n, skew.p);
        rst       <- rst - mu + rnorm(nPat, 0, noise.sig);
        rst       <- rst/sqrt(va + noise.sig^2)*ysig;
    });

    rst
}

#' Simulate outcomes for a single arm study
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuOneArm <- function(nPat, muCov, sdCov, corCov, regCoeff, fmla = NULL,
                          type = c("continuous", "binary"),
                          ysig = NULL, sig2Ratio=2, bin.mu = 0.5,
                          ...) {

    type  <- match.arg(type);
    COV.X <- rweSimuCov(nPat, muCov, sdCov, corCov);

    ##simulate Y
    if ("continuous" == type) {
        ## epsilon
        if (is.null(ysig)) {
            ysig <- rweGetYSig(sig2Ratio,
                               muCov    = muCov,
                               sdCov    = sdCov,
                               corCov   = corCov,
                               regCoeff = regCoeff,
                               fmla     = fmla);
        }
        XBETA   <- rweXBeta(cov.x    = COV.X,
                            regCoeff = regCoeff,
                            fmla     = fmla);

        EPSILON <- rweSimuError(nPat,
                                ysig = ysig,
                                ...);

        Y <- XBETA + EPSILON;
    } else if ("binary" == type) {
        stopifnot(!is.null(bin.mu));
        b0  <- rweGetBinInt(bin.mu,
                            muCov    = muCov,
                            sdCov    = sdCov,
                            corCov   = corCov,
                            regCoeff = regCoeff,
                            fmla     = fmla);

        regCoeff <- c(b0, regCoeff);
        XBETA    <- rweXBeta(cov.x    = COV.X,
                             regCoeff = regCoeff,
                             fmla     = fmla);

        Y <- rbinom(nPat, 1, expit(XBETA));
    }

    ##return
    Data           <- cbind(1:nPat, Y, COV.X);
    colnames(Data) <- c("pid", "Y", paste("V", 1:ncol(COV.X), sep=""));
    data.frame(Data);
}

#' Simulate continuous outcomes for a two arm study
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuTwoArm <- function(nPat, muCov, sdCov, corCov, regCoeff.y,
                          regCoeff.z=0,
                          fmla.y = NULL, fmla.z = NULL, ysig = NULL, b0 = NULL,
                          sig2Ratio = 2, bin.mu = 0.5, ..., do.simu=TRUE) {

    ##treatment assignment
    if (is.null(b0)) {
        b0  <- rweGetBinInt(bin.mu,
                            nPat     = nPat,
                            muCov    = muCov,
                            sdCov    = sdCov,
                            corCov   = corCov,
                            regCoeff = regCoeff.z,
                            fmla     = fmla.z);
    }

    if (is.null(ysig)) {
        ysig <- rweGetYSig(sig2Ratio,
                           nPat     = nPat,
                           muCov    = muCov,
                           sdCov    = sdCov,
                           corCov   = corCov,
                           regCoeff = regCoeff.y[-1],
                           fmla     = fmla.y);
    }

    if (do.simu) {
        ##covariates
        COV.X   <- rweSimuCov(nPat, muCov, sdCov, corCov);
        xbeta.z <- rweXBeta(cov.x    = COV.X,
                            regCoeff = regCoeff.z,
                            fmla     = fmla.z);
        Z       <- rbinom(nPat, 1, expit(b0+xbeta.z));
        xbeta.y <- rweXBeta(cov.x    = COV.X,
                            regCoeff = regCoeff.y[-1],
                            fmla     = fmla.y);
        epsilon.y <- rweSimuError(nPat, ysig = ysig, ...);
        Y         <- Z * regCoeff.y[1] + xbeta.y + epsilon.y;
        ##return
        rst       <- cbind(pid=1:nPat, Y=Y, Z=Z, COV.X);
    } else {
        rst <- c(b0=b0, ysig=ysig);
    }
    rst
}

#' Get unbalance in baseline covariates
#'
#' @inheritParams simupara
#'
#' @export
rweUnbalance <- function(nPat, ..., pts = NULL) {
    if (is.null(pts)) {
        pts <- rweSimuTwoArm(nPat, ...);
    }

    ##unbalance
    inx.0 <- which(0 == pts[,"Z"]);
    c.xy  <- colnames(pts);
    c.xy  <- c.xy[grep("^V[0-9]+$", c.xy)];
    unb   <- NULL;
    for (i in 1:length(c.xy)) {
        x0     <- sample(pts[inx.0,  c.xy[i]], size = nPat, replace = TRUE);
        x1     <- sample(pts[-inx.0, c.xy[i]], size = nPat, replace = TRUE);
        x.diff <- x1 - x0;
        unb    <- rbind(unb, data.frame(V=c.xy[i], Diff=x1 - x0));
    }
    unb
}

#' Simulate data from an existing dataset
#'
#'
#'
#' @export
#'
rweSimuFromTrial <- function(trial.data, sizes,
                             group = "A", outcome = "Y", keep.group = TRUE,
                             trt.effect = 0,
                             with.replacement = FALSE,
                             seed = NULL) {
    if (!is.null(seed))
        set.seed(seed);

    if (1 == length(sizes)) {
        sizes <- rep(sizes,2);
    }

    if (keep.group) {
        arms   <- unique(trial.data[[group]]);
        rst    <- NULL;
        mean.y <- NULL;
        for (i in 1:length(arms)) {
            cur.n   <- sizes[min(i, length(sizes))];
            cur.d   <- trial.data[arms[i] == trial.data[[group]],];
            mean.y  <- c(mean.y, mean(cur.d[[outcome]]));
            if (nrow(cur.d) > cur.n) {
                wr <- with.replacement;
            } else {
                wr <- TRUE;
            }
            cur.smp <- cur.d[sample(1:nrow(cur.d), cur.n, replace=wr), ];
            cur.smp[[group]] <- i - 1;
            rst <- rbind(rst, cur.smp);
        }

        trt.effect <- mean.y[length(mean.y)] - mean.y[1];
        simu.data  <- rst;
    } else {
        cur.d <- trial.data;
        cur.n <- sum(sizes);
        if (nrow(cur.d) > cur.n) {
            wr <- with.replacement;
        } else {
            wr <- TRUE;
        }
        cur.smp <- cur.d[sample(1:nrow(cur.d), cur.n, replace=wr), ];

        for (i in 1:length(sizes)) {
            if (1 == i) {
                tmp.start <- 1
            } else {
                tmp.start <- sum(sizes[1:(i-1)]);
            }

            tmp.end <- sum(sizes[1:i]);
            cur.smp[tmp.start:tmp.end,group] <- i-1;
        }

        last.arm <- (length(sizes)-1) == cur.smp[[group]];
        cur.smp[last.arm, outcome] <- cur.smp[last.arm, outcome] + trt.effect;

        trt.effect <- trt.effect;
        simu.data  <- cur.smp;
    }

    ## randomize the order of simulated patients
    list(true.effect = trt.effect,
         simu.data   = simu.data[sample(1:nrow(simu.data)),]);
}
