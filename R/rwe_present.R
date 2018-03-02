
#' Export results into a template file
#'
#' @param numbers    vector of results
#' @param template.f template file name
#' @param out.f      output file name
#' @param sub.str    pattern of string to be replaced
#'
#'
#' @export
#'
rweExpRst <- function(numbers, template.f,  out.f="rst.txt", sub.str="AA") {
    if (!file.exists(template.f)) {
        return;
    }
    ##read
    tpla <- readChar(template.f, file.info(template.f)$size);

    ##substitute
    for (i in 1:length(numbers)) {
        tpla <- sub(sub.str, numbers[i], tpla);
    }

    ##write out
    write(tpla, file=out.f);
}


#' Plot unbalance in covariates
#'
#'
#' @export
#'
rwePlotUnbalance <- function(data.unb,
                             var.x = "Diff",
                             var.group = NULL,
                             xlim = NULL,
                             ylim = NULL,
                             title = "",
                             f.grid = formula("V~Study"),
                             adjust = 1) {

    if (is.null(var.group)) {
        rst <- ggplot(data.unb, aes_string(x=var.x)) +
            geom_density(alpha = 0.4, fill = "gray", na.rm = TRUE, adjust = adjust) +
            geom_vline(xintercept = 0, linetype="dashed", col = "red");
    } else {
        rst <- ggplot(data.unb,
                      aes_string(x=var.x,
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
