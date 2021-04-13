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
                        m_var      = c("jk", "bs"),
                        pred_tps   = 1,
                        nbs        = 100,
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
    rst_est <- NULL
    rst_jk  <- NULL
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

        ## estimate
        cur_est     <- f_surv(cur_data, cur_weights)

        ## variances
        if ("jk" == m_var) {
            ##jackknife
            cur_jk      <- NULL;
            cur_weights <- c(rep(1, ns1 - 1), rep(lambdas[i] / ns0, ns0))
            for (j in 1:ns1) {
                cur_data  <- rbind(cur_d1[-j, ], cur_d0)
                jk_est    <- f_surv(cur_data, cur_weights)
                cur_jk    <- rbind(cur_jk, (jk_est - cur_est)^2)
                rst_jk    <- rbind(rst_jk, c(i, jk_est))
            }

            if (ns0 > 0 & lambdas[i] > 0) {
                cur_weights <- c(rep(1, ns1), rep(lambdas[i] / ns0, ns0 - 1))
                for (j in 1:ns0) {
                    cur_data <- rbind(cur_d1, cur_d0[-j, ])
                    jk_est   <- f_surv(cur_data, cur_weights)
                    cur_jk   <- rbind(cur_jk, (jk_est - cur_est)^2)
                    rst_jk   <- rbind(rst_jk, c(i, jk_est))
                }
            }

            cur_n   <- nrow(cur_jk)
            cur_var <- apply(cur_jk, 2, sum)
            cur_var <- cur_var * (cur_n - 1) / cur_n

        } else if ("bs" == m_var) {
            ## bootstrap
            cur_bs <- NULL
            for (j in 1:nbs) {
                cur_data <- rbind(cur_d1[sample(1:ns1, replace = T), ],
                                  cur_d0[sample(1:ns0, replace = T), ])
                bs_est   <- f_surv(cur_data, cur_weights)
                cur_bs   <- rbind(cur_bs, bs_est)
            }
            cur_var <- apply(cur_bs, 2, var)
        }

        ## append
        rst_est <- rbind(rst_est,
                         c(ns1, ns0, cur_est, cur_var))
    }

    ## summarize the results
    ws          <- rst_est[, 1] / sum(rst_est[, 1])
    est_strata  <- rst_est[, 2 + (1:n_est),         drop = F]
    var_strata  <- rst_est[, 2 + n_est + (1:n_est), drop = F]
    est_overall <- apply(est_strata, 2, function(x) sum(ws * x))
    est_var     <- apply(var_strata, 2, function(x) sum(ws^2 * x))

    est_var_2 <- NULL
    if ("jk" == m_var) {
        ## overall JK variance
        jk_var <- apply(rst_jk, 1, function(x) {
            est         <- est_strata
            est[x[1], ] <- x[-1]
            est         <- apply(est, 2, function(t) sum(ws * t))
            (est - est_overall)^2
        })

        ## switch variance to overall JK approach
        est_var_2 <- est_var

        est_var   <- apply(jk_var, 1, sum)
        est_var   <- est_var_2 * (nrow(rst_jk) - 1) / nrow(rst_jk)
    }

    ## reset seed
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    list(est         = est_overall,
         var         = est_var,
         var_jk      = est_var_2,
         est_strata  = est_strata,
         var_strata  = var_strata,
         ns1         = rst_est[, 1],
         ns0         = rst_est[, 2])
}
