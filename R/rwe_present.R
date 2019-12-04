#' Plot unbalance in covariates
#'
#'
#' @export
#'
rwePlotUnbalance <- function(data.unb,
                             var.x     = "Diff",
                             var.group = NULL,
                             xlim      = NULL,
                             ylim      = NULL,
                             title     = "",
                             f.grid    = formula("V~Study"),
                             adjust    = 1) {

    if (is.null(var.group)) {
        rst <- ggplot(data.unb, aes_string(x=var.x)) +
            geom_density(alpha = 0.4, fill = "gray", na.rm = TRUE, adjust = adjust) +
            geom_vline(xintercept = 0, linetype="dashed", col = "red");
    } else {
        rst <- ggplot(data.unb,
                      aes_string(x     = var.x,
                                 group = var.group,
                                 color = var.group,
                                 linetype = var.group)) +
            geom_density(alpha = 0.2, fill = "gray", na.rm = TRUE, adjust = adjust, );
    }

    rst <- rst + labs(x = "", y="", title=title) +
        theme_bw() +
        theme(strip.background = element_blank(),
              panel.grid       = element_blank(),
              panel.border     = element_blank(),
              axis.text.y      = element_blank(),
              axis.text.x      = element_text(size=9),
              axis.ticks.y     = element_blank(),
              panel.spacing    = unit(0, "lines"))+
        facet_grid(f.grid);

    if (!is.null(xlim))
        rst <- rst + scale_x_continuous(limits = xlim, breaks = 0);
    if (!is.null(ylim))
        rst <- rst + scale_y_continuous(limits = ylim);

    rst
}


#' Plot PS distributions
#'
#'
#' @method plot RWE_DWITHPS
#'
#' @export
#'
plot.RWE_DWITHPS <- function(x, type = c("ps", "balance"), ...) {
    type <- match.arg(type);

    switch(type,
           ps = plotRwePs(x, ...),
           balance = plotRweBalance(x, ...));
}


#' add a method summary2
#'
#'
#' @export
#'
summary2 <- function (x, ...) {
    UseMethod("summary2", x)
}

#' Summary unbalance metrics
#'
#'
#' @export
#'
#' @method summary2 RWE_DWITHPS
#'
summary2.RWE_DWITHPS <- function(x, v.cov = NULL, label.cov = v.cov, ...) {
    if (is.null(v.cov))
        v.cov <- all.vars(x$ps.fml)[-1];

    if (is.null(label.cov))
        label.cov <- v.cov;

    nstrata      <- x$nstrata;
    dtaps        <- x$data;
    dtaps$Strata <- dtaps[["_strata_"]];
    dtaps$Group  <- dtaps[["_grp_"]];

    rst <- NULL;
    for (iv in 1:length(v.cov)) {
        v        <- v.cov[iv];
        cur.cov  <- dtaps[[v]];
        if (is.factor(cur.cov)) {
            cur.cov <- as.numeric(cur.cov);
        }

        for (stra in 1:max(dtaps$Strata, na.rm=TRUE)) {
            cur.cov0 <- cur.cov[which(0 == dtaps$Group & stra == dtaps$Strata)];
            cur.cov1 <- cur.cov[which(1 == dtaps$Group & stra == dtaps$Strata)];

            cur.met  <- rweBalMetric(cur.cov0, cur.cov1, ...);
            rst      <- rbind(rst,
                              data.frame(Covariate = label.cov[iv],
                                         Strata    = stra,
                                         Metric    = cur.met)
                              );
        }

        ##overall
        cur.cov0 <- cur.cov[which(0 == dtaps$Group & !is.na(dtaps$Strata))];
        cur.cov1 <- cur.cov[which(1 == dtaps$Group & !is.na(dtaps$Strata))];
        cur.met  <- rweBalMetric(cur.cov0, cur.cov1, ...);
        rst      <- rbind(rst,
                          data.frame(Covariate = label.cov[iv],
                                     Strata    = 0,
                                     Metric    = cur.met)
                          );

    }

    rst
}
