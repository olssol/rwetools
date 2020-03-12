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

## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------
##             PRIVATE FUNCTIONS
## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------

## plot density of propensity score
plotRwePs <- function(data.withps, overall.inc = TRUE, add.text = TRUE,
                      facet.scales = "free_y", ...) {
    stopifnot(inherits(data.withps,
                       what = get.rwe.class("DWITHPS")));

    pskl        <- rwePSDist(data.withps, ...);
    nstrata     <- data.withps$nstrata;
    dtaps       <- data.withps$data;

    pskl$Strata <- as.factor(c(paste("Stratum ", 1:nstrata, sep = ""), "Overall"));
    xlim        <- range(dtaps[which(!is.na(dtaps[["_strata_"]])),"_ps_"], na.rm = TRUE);

    all.data <- NULL;
    for (i in 1:nstrata) {
        cur.sub  <- dtaps[which(i == dtaps[["_strata_"]]),];
        cur.data <- data.frame(Strata = paste("Stratum ", i, sep = ""),
                               Ps     = cur.sub[["_ps_"]],
                               Group  = cur.sub[["_grp_"]]);
        all.data <- rbind(all.data, cur.data);
    }

    if (overall.inc) {
        cur.data <-  data.frame(Strata = "Overall",
                                Ps     = dtaps[["_ps_"]],
                                Group  = dtaps[["_grp_"]]);

        all.data <- rbind(all.data, cur.data);
    } else {
        pskl <- pskl %>% filter(Strata != "Overall");
    }

    all.data$Group <- as.factor(all.data$Group);
    rst <- ggplot(data = all.data, aes(x = Ps)) +
        geom_density(alpha = 0.2,
                       aes(group = Group,
                           fill  = Group,
                           linetype = Group),
                     trim  = TRUE,
                     na.rm = TRUE) +
        labs(x = "Propensity Score", y = "Density") +
        scale_y_continuous(breaks = NULL) +
        scale_x_continuous(limits = xlim) +
        scale_fill_manual(values=c("gray20", "gray80")) +
        theme_bw() +
        theme(strip.background = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_rect(colour = "black"),
              panel.spacing = unit(0, "lines")) +
        facet_grid(Strata ~ ., scales = facet.scales);

    if (add.text) {
        rst <- rst +
            geom_text(x = Inf, y = Inf, hjust = 1, vjust = 1,
                      aes(label = paste('n0=', N0,
                                        ", n1=", N1,
                                        ", OVL=", format(Dist, digits = 3),
                                        sep = "")),
                      data = pskl, size = 4);
    }
    rst
}

plot.balance.fac <- function(dtaps, v, overall.inc = TRUE) {
    cur.d <- rweFreqTbl(dtaps, var.groupby = c("Strata", "Group"), vars = v);
    cur.d <- cur.d %>% dplyr::filter(!is.na(Strata));
    cur.d$Strata <- paste("Stratum ", cur.d$Strata, sep = "")
    if (overall.inc) {
        cur.overall <- rweFreqTbl(dtaps, var.groupby = "Group", vars = v);
        cur.overall$Strata <- "Overall";
        cur.d <- rbind(cur.d, cur.overall);
    }
    cur.d$Group <- as.factor(cur.d$Group);
    cur.d$Value <- as.factor(cur.d$Value);

    rst <- ggplot(data = cur.d, aes(x = Value, y = Freq)) +
        geom_bar(alpha = 0.4,
                 stat = "identity",
                 position = "dodge",
                 color = "black",
                 aes(group = Group,
                     fill  = Group)) +
        scale_fill_manual(values=c("gray20", "gray80")) +
        scale_y_continuous(breaks = NULL, limits = c(0,1)) +
        labs(x = "", y = "") +
        facet_grid(Strata ~ .);
    rst
}

plot.balance.cont <- function(dtaps, v, nstrata, overall.inc = TRUE, facet.scales = "free_y") {
    cur.d <- NULL;
    for (i in 1:nstrata) {
        cur.sub      <- dtaps[which(i == dtaps[["_strata_"]]),];
        cur.v        <- data.frame(Cov    = v,
                                   Value  = cur.sub[[v]],
                                   Group  = cur.sub[["_grp_"]]);
        cur.v$Strata <- paste("Stratum ", i, sep = "");
        cur.d        <- rbind(cur.d, cur.v);
    }

    if (overall.inc) {
        cur.sub      <- dtaps;
        cur.v        <- data.frame(Cov    = v,
                                   Value  = cur.sub[[v]],
                                   Group  = cur.sub[["_grp_"]]);
        cur.v$Strata <- paste("Overall");
        cur.d        <- rbind(cur.d, cur.v);
    }
    cur.d$Group <- as.factor(cur.d$Group);

    rst <- ggplot(data = cur.d, aes(x = Value)) +
        geom_density(alpha = 0.2,
                     aes(group = Group,
                         fill  = Group,
                         linetype = Group),
                     na.rm = TRUE) +
        scale_y_continuous(breaks = NULL) +
        scale_fill_manual(values=c("gray20", "white")) +
        labs(x = "", y = "") +
        facet_grid(Strata ~ ., scales = facet.scales);
    rst
}


plotRweBalance <- function(data.withps, overall.inc = TRUE, v.cov = NULL,
                           facet.scales = "free_y", label.cov = v.cov, legend.width = 0.08,
                           ...) {

    if (is.null(v.cov))
        v.cov <- all.vars(data.withps$ps.fml)[-1];

    if (is.null(label.cov))
        label.cov <- v.cov;

    nstrata      <- data.withps$nstrata;
    dtaps        <- data.withps$data;
    dtaps$Strata <- dtaps[["_strata_"]];
    dtaps$Group  <- dtaps[["_grp_"]];

    rst <- list();
    for (v in v.cov) {
        if (is.factor(dtaps[[v]])) {
            cur.p <- plot.balance.fac(dtaps, v, overall.inc = overall.inc);
        } else {
            cur.p <- plot.balance.cont(dtaps, v, nstrata = nstrata,
                                       overall.inc = overall.inc, facet.scales = facet.scales);
        }
        cur.p <- cur.p +
            labs(title = label.cov[v == v.cov]) +
            theme_bw() +
            theme(strip.background = element_blank(),
                  strip.placement  = "right",
                  strip.text       = element_blank(),
                  panel.grid       = element_blank(),
                  panel.border     = element_blank(),
                  panel.spacing    = unit(0, "lines"),
                  plot.title = element_text(hjust = 0.5),
                  legend.position  = "none",
                  plot.margin      = unit(c(1,0,1,-0.5), "lines"));

        rst[[v]] <- cur.p;
    }

    rst[[length(rst)]] <- rst[[length(rst)]] +
      theme(strip.text = element_text(size=8),
            legend.position = "right")

    rst$nrow        <- 1;
    rst$rel_widths <- c(rep(1, length(v.cov)-1),
                        1+legend.width*length(v.cov));
    do.call(plot_grid, rst);
}
