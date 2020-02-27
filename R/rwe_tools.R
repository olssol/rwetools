#' Cut a sequence of numbers into bins with equal numbers in each bin
#'
#' @param x vector of values based on which cut points will be determined
#' @param y vector of values to be cut, default to be the same as \code{x}
#' @param breaks number of cut points
#' @param keep.inx indices of y that will be categorized as 1 or the largest bin
#'     even if their values are out of range of x
#'
#' @export
#'
rweCut <- function(x, y = x, breaks = 5, keep.inx = NULL) {
    cuts    <- quantile(x, seq(0, 1,length = breaks+1));
    cuts[1] <- cuts[1] - 0.001;
    rst     <- rep(NA, length(y));
    for (i in 2:length(cuts)) {
        inx      <- which(y > cuts[i-1] & y <= cuts[i]);
        rst[inx] <- i-1;
    }

    if (!is.null(keep.inx)) {
        inx <- which(y[keep.inx] <= cuts[1]);
        if (0 < length(inx)) {
            rst[keep.inx[inx]] <- 1;
        }

        inx <- which(y[keep.inx] > cuts[length(cuts)]);
        if (0 < length(inx)) {
            rst[keep.inx[inx]] <- length(cuts) - 1;
        }
    }

    rst
}


#' Compute distance from F0 to F1
#'
#' @param type type of distances. ovl: overlapping coefficient, kl:
#'     1/(1+Kullback-Leibler divergence)
#' @param n.bins number of bins for KL computation
#' @param epsilon small integer for Dirichlet smoothing
#'
#' @return a vector with the number of samples in group 0, the number of samples
#'     in group 1, and 1/(1+KL divergence) from group 0 to group 1 when type is
#'     kl, or the overlapping coefficient when type is ovl
#' @export
#'
rweDist <- function(sample.F0, sample.F1, n.bins = 10, type = c("ovl", "kl"), epsilon = 10^-6) {

    type     <- match.arg(type);

    smps     <- c(sample.F0, sample.F1);
    n0       <- length(sample.F0);
    n1       <- length(sample.F1);

    if (0 == n0 | 0 == n1)
        return(c(n0, n1, NA));

    if (1 == length(unique(smps))) {
        cut.smps <- rep(1, n0+n1)
        n.bins   <- 1;
        warning("Distributions for computing distances are degenerate.",
                call. = FALSE);
    } else {
        cut.smps <- rweCut(smps, breaks = n.bins);
    }

    rst <- 0;
    for (j in 1:n.bins) {
        n0.j <- length(which(j == cut.smps[1:n0]));
        n1.j <- length(which(j == cut.smps[(n0+1):(n0+n1)]));

        rst  <- rst + switch(type,
                             kl = {ep0  <- (n0.j+epsilon)/(n0 + epsilon * n.bins);
                                 ep1  <- (n1.j+epsilon)/(n1 + epsilon * n.bins);
                                 ep1 * log(ep1/ep0)},
                             ovl = min(n0.j/n0, n1.j/n1));
    }

    if ("kl" == type)
        rst <- 1/(1+rst);

    rst;
}


#' Generate frequency table for factor columns
#'
#' @return a vector with the number of samples in group 0, the number of samples
#'     in group 1, and the KL divergence from group 0 to group 1
#' @export
#'
rweFreqTbl <- function(data, var.groupby, vars = NULL) {

    if (is.null(vars))
        vars <- colnames(data);

    rst <- NULL;
    for (v in vars) {
        if (!is.factor(data[[v]]))
            next;

        cur.freq <- data %>% count_(c(var.groupby, v)) %>%
            group_by_(.dots = var.groupby) %>%
            mutate(Sum = sum(n), Freq = n/sum(n)) %>%
            mutate_if(is.factor, as.character) %>%
            mutate(Cov = v) %>%
            rename_(Value = v);

        rst <- rbind(rst, data.frame(cur.freq));
    }

    rst
}

#' Split A into Bins
#'
#' Split A into bins to minimize the number difference between two arms
#'
#' @param ns1.trt numbers of treatment arm patients in each stratum
#' @param ns1.ctl numbers of control arm patients in each stratum
#' @param ns0     numbers of external patients in each stratum
#' @param A       number of target patients to borrow
#'
#' @export
#'
rweEvenLmbdS <- function(ns1.trt, ns1.ctl, ns0, A, init.lmbds = NULL) {

    f.target <- function(lmbds) {
        cl  <- c(lmbds, A - sum(lmbds));
        rst <- sum((ns1.trt - ns1.ctl - cl)^2);
        rst
    }

    f.gradient <- function(lmbds) {
        g <- -2*(ns1.trt[-NS] - ns1.ctl[-NS] - lmbds);
    }

    NS <- length(ns1.trt);

    stopifnot(NS == length(ns1.ctl) &
              NS == length(ns0));


    ## there is only one stratum
    if (1 == NS) {
        return(A);
    }

    ## multiple strata

    ## restrictions lmb > 0; sum(lmb) < A; lmb < ns0;
    ui <- rbind(diag(NS-1),
                -1 * diag(NS-1),
                rep(-1, NS-1),
                rep(1, NS-1));

    ci <- c(rep(0, NS-1),
            -ns0[-NS],
            -A,
            A - ns0[NS]);

    if (is.null(init.lmbds))
        init.lmbds <- rep(A/NS, NS-1);

    rst <- constrOptim(theta = init.lmbds,
                       f     = f.target,
                       grad  = f.gradient,
                       ui, ci, mu = 1e-04, control = list(),
                       outer.iterations = 100, outer.eps = 1e-05,
                       hessian = FALSE);

    c(rst$par, A - sum(rst$par))
}

#' Get weights
#'
#' @param A target number of subjects to be borrowed
#' @param m.lambda method to split A. rs: by overlapping coefficient; even: by
#'     minimizing trt and control imbalance in numbers
#'
#' @return power parameter before standardization
#'
#' @export
#'
rweGetLambda <- function(A, rs = NULL, ns1.trt = NULL, ns1.ctl = NULL, ns0,
                         m.lambda = c("rs", "even", "inverse"), ...) {
    m.lambda <- match.arg(m.lambda);

    rst <- switch(m.lambda,
                  rs      = apply(cbind(ns0, A * rs/sum(rs)), 1, min),
                  even    = rweEvenLmbdS(ns1.trt, ns1.ctl, ns0, A, ...),
                  inverse = {mrs <- 1/(1-rs);
                             apply(cbind(ns0, A * mrs/sum(mrs)), 1, min)})
    rst
}


#' Get weights
#'
#' @param A target number of subjects to be borrowed
#' @param m.lambda method to split A. rs: by overlapping coefficient; even: by
#'     minimizing trt and control imbalance in numbers
#'
#' @return power parameter before standardization
#'
#' @export
#'
rweGpsLambda <- function(A, ps_dist) {
    ps_dist <- ps_dist[which(ps_dist$Strata > 0), ]
    inx     <- seq(4, ncol(ps_dist), by = 2)

    ## average distance
    avg   <- apply(ps_dist[, inx], 1, mean)
    avg_a <- A * avg / sum(avg)

    ## lambdas
    lambdas <- apply(cbind(avg_a, ps_dist[, inx]),
                     1,
                     function(x) {
                         x[-1] / sum(x[-1]) * x[1]
                     })

    ## minimum
    rst <- apply(cbind(as.vector(t(lambdas)),
                       as.vector(data.matrix(ps_dist[, inx-1]))
                      ),
                 1,
                 min)

    matrix(rst, nrow = nrow(ps_dist))
}


#' Summary statistics
#' l
#'
#'
#'
#'
#' @export
#'
rweSummary <- function(cur.m, cur.var, true.theta,  cur.ci = NULL) {
    rst      <- c(thetahat = cur.m,
                  thetavar = cur.var,
                  bias     = cur.m - true.theta,
                  mse      = (cur.m - true.theta)^2)

    if (!is.null(cur.ci)) {
        range.ci <- range(cur.ci)
        rst <- c(rst,
                 width    = range.ci[2] - range.ci[1],
                 cover    = true.theta >= range.ci[1] & true.theta <= range.ci[2],
                 lb       = as.numeric(range.ci[1]),
                 ub       = as.numeric(range.ci[2]))
    }
    rst
}
