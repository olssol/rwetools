#' Call STAN models
#'
#'
#' @param chains STAN parameter. Number of Markov chainsm
#' @param iter STAN parameter. Number of iterations
#' @param warmup STAN parameter. Number of burnin.
#' @param control STAN parameter. See \code{rstan::stan} for details.
#' @param ... other options to call STAN sampling such as \code{thin},
#'     \code{algorithm}. See \code{rstan::sampling} for details.#'
#'
#'
#' @export
#'
rweSTAN <- function(lst.data, stan.mdl = "powerp",
                    chains = 4, iter = 5000, warmup = 1000,
                    control = list(adapt_delta=0.95), ...) {

    stan.rst <- rstan::sampling(stanmodels[[stan.mdl]],
                                data    = lst.data,
                                chains  = chains,
                                iter    = iter,
                                warmup  = warmup,
                                control = control,
                                ...);

    stan.rst;
}



#' Get Posterior for all stratum
#'
#'
#' @param data class DWITHPS data frame
#' @param As   power of the power prior for each strata
#' @param ...  extra parameters for calling function \code{\link{rweSTAN}}
#'
#' @export
#'
rweDrawPost <- function(data, v.outcome = "Y", As=0, type = c("continuous", "binary"), ...) {
    stopifnot("RWE_DWITHPS" %in% class(data));
    stopifnot(v.outcome %in% colnames(data));

    type     <- match.arg(type);
    stan.mdl <- switch(type,
                       continuous = "powerp",
                       binary     = "powerpbinary");

    data      <- data[!is.na(data[["_strata_"]]),];
    nstrata   <- max(data[["_strata_"]]);
    rst.theta <- NULL;
    for (i in 1:nstrata) {
        cur.d1 <- data[data[["_strata_"]] == i & data[["_grp_"]] == 1, v.outcome];
        cur.d0 <- data[data[["_strata_"]] == i & data[["_grp_"]] == 0, v.outcome];

        if (0 == length(cur.d1))
            next;

        Y1 <- cur.d1;
        N1 <- length(cur.d1);
        N0 <- length(cur.d0);
        if (0 == N0) {
            ##dummy value
            Y0 <- 0;
        } else {
            Y0 <- cur.d0;
        }

        cur.A     <- As[min(i, length(As))];
        cur.post  <- rweSTAN(lst.data = list(A=cur.A, N0=N0, N1=N1, Y0=Y0, Y1=Y1),
                             stan.mdl = stan.mdl, ...);
        cur.theta <- rstan::extract(cur.post, pars = "theta")$theta;
        rst.theta <- rbind(rst.theta, cur.theta);
    }

    rst.theta
}

#' Summary Posterior theta
#'
#'
#' @param post.theta posterior samples from STAN
#' @param true.theta true value of theta
#' @param quants     quantiles
#'
#' @export
#'
rweSummaryPost <- function(post.theta, true.theta, quants = c(0.025, 0.975)) {

    cur.post        <- apply(post.theta, 2, mean);
    cur.m           <- mean(cur.post);
    cur.ci          <- quantile(cur.post, quants);
    post.var        <- apply(post.theta, 1, var);
    names(post.var) <- paste("var", 1:length(post.var), sep = "");
    range.ci        <- range(cur.ci);

    rst <- c(thetahat = cur.m,
             thetavar = var(cur.post),
             bias     = cur.m - true.theta,
             mse      = (cur.m - true.theta)^2,
             width    = range.ci[2] - range.ci[1],
             cover    = true.theta >= range.ci[1] & true.theta <= range.ci[2],
             post.var,
             cur.ci);
    rst

}
