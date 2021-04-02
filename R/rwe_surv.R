#' PS-Integrated Survival Curve
#'
#' For all stratum. Variance estimated by Jack-Knife.
#'
#' @param data class DWITHPS data frame
#' @param ... parameters for \code{rweWL}
#' @param lambdas power parameter without standardization by ns0
#' @param m_var method to get variance: jackknife or bootstrap
#' @param seed random seed
#'
#' @export
#'
rwe_ps_surv <- function(data,
                        lambdas,
                        v_censored = "censored",
                        v_event    = "time",
                        m_var      = c("jk"),
                        pred_tps   = 1,
                        seed       = NULL, ...) {

    f_surv <- function(cur_data, cur_weights) {
        cur_surv <- survfit(fml_surv,
                            data    = cur_data,
                            weights = cur_weights)

        rst <- summary(cur_surv, time = pred_tps)$surv
        rst
    }

    stopifnot(all(c(v_event, v_censored) %in%
                  colnames(data)))

    m_var <- match.arg(m_var)

    ## number of thetas
    n_est <- length(pred_tps)

    ## set random seed
    if (!is.null(seed))
        old_seed <- set.seed(seed)

    ## prepare data
    data               <- data[!is.na(data[["_strata_"]]), ]
    data[[v_censored]] <- 1 - data[[v_censored]]
    S                  <- max(data[["_strata_"]])

    ## time points
    if (is.null(pred_tps)) {
        pred_tps <- data[which(1 == data[[v_censored]]), v_event]
        pred_tps <- sort(unique(pred_tps))
    }

    ## formula
    fml_surv <- paste("Surv(", v_event, ",",
                      v_censored, ") ~ 1", sep = "")
    fml_surv <- as.formula(fml_surv)

    ## find survival estimate
    rst_est <- NULL;
    for (i in 1:S) {
        cur_d1 <- data[data[["_strata_"]] == i & data[["_grp_"]] == 1, ]
        cur_d0 <- data[data[["_strata_"]] == i & data[["_grp_"]] == 0, ]

        ns1 <- nrow(cur_d1)
        ns0 <- nrow(cur_d0)
        if (0 == ns1) {
            stop(paste("Stratum ", i,
                       " contains no subjects from group 1",
                       sep = ""))
        }

        cur_data    <- rbind(cur_d1, cur_d0)
        cur_weights <- c(rep(1, ns1), rep(lambdas[i] / ns0, ns0))
        cur_est     <- f_surv(cur_data, cur_weights)

        cur_jk  <- NULL;
        ##jackknife
        if ("jk" == m_var) {
            cur_weights <- c(rep(1, ns1 - 1), rep(lambdas[i] / ns0, ns0))
            for (j in 1:ns1) {
                cur_data  <- rbind(cur_d1[-j, ], cur_d0)
                jk_est    <- f_surv(cur_data, cur_weights)
                cur_jk    <- rbind(cur_jk, (jk_est - cur_est)^2)
            }

            if (ns0 > 0 & lambdas[i] > 0) {
                cur_weights <- c(rep(1, ns1), rep(lambdas[i] / ns0, ns0 - 1))
                for (j in 1:ns0) {
                    cur_data <- rbind(cur_d1, cur_d0[-j, ])
                    jk_est   <- f_surv(cur_data, cur_weights)
                    cur_jk   <- rbind(cur_jk, (jk_est - cur_est)^2)
                }
            }
        }
        cur_n      <- nrow(cur_jk)
        cur_jk_var <- apply(cur_jk, 2, sum) * (ns1 + ns0 - 1) / (ns1 + ns0)
        rst_est    <- rbind(rst_est, c(ns1, ns0, cur_est, cur_jk_var))
    }

    ## summarize the results
    ws          <- rst_est[, 1] / sum(rst_est[,1])
    est_strata  <- rst_est[, 2 + (1:n_est),         drop = F]
    var_strata  <- rst_est[, 2 + n_est + (1:n_est), drop = F]
    est_overall <- apply(est_strata, 2, function(x) sum(ws * x))
    est_var     <- apply(var_strata, 2, function(x) sum(ws^2 * x))

    ## reset seed
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    list(est         = est_overall,
         var         = est_var,
         est_strata  = est_strata,
         var_strata  = var_strata,
         ns1         = rst_est[, 1],
         ns0         = rst_est[, 2])
}
