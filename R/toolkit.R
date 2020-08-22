
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
tkExpRst <- function(numbers, template.f,  out.f="rst.txt", sub.str="AA",
                     append = TRUE) {
    if (!file.exists(template.f)) {
        return;
    }
    ##read
    tpla <- readChar(template.f, file.info(template.f)$size);

    ##substitute
    for (i in 1:length(numbers)) {
        tpla <- sub(sub.str, numbers[i], tpla)
    }

    ##write out
    write(tpla, file=out.f, append = append)
}


#' Import objects in a list into a designated environment
#'
#' @param alist list of objects
#' @param dest.env designated environment
#'
#' @export
#'
tkMakeLocal <- function(alist, dest.env='.GlobalEnv') {
    for (i in 1:length(alist)) {
        assign(names(alist[i]), alist[[i]], dest.env);
    }
}


#' Call function by its name organized as a vector
#'
#' @param vec function names as a vector
#'
#' @export
#'
tkCallFun <- function(vec, ...) {
    eval(parse(text=paste("rst <- ",
                          paste(vec, collapse = ""),
                          "(...)",
                          sep = "")
               )
         );
    rst
}


#' Extract summary statistics from a data frame
#'
#' @param dta data frame
#' @param quants quantiles for continuous variables
#'
#' @return a list of summary statistics
#'
#' @export
#'
tk_extract_stats <- function(x, quants = c(0.5), xlev = NULL,
                             type = NULL, weights  = NULL) {

    f_msd <- function(x, weights) {
        x_mean  <- sum(x * weights)
        x2_mean <- sum(x^2 * weights)

        x_sd    <- x2_mean - (x_mean)^2
        x_sd    <- sqrt(x_sd)

        c(x_mean, x2_mean, x_sd)
    }

    f_quantile <- function(x, quants) {
        ## x_quants <- Hmisc::wtd.quantile(x, weights = weights,
        ##                               probs = quants, normwt = TRUE)

        x_quants <- quantile(x, probs = quants)
        cbind(quants, x_quants)
    }

    f_factor <- function(x, xlev, weights) {
        x <- as.factor(x)
        if (!is.null(xlev))
            levels(x) <- xlev

        probs <- NULL
        for (i in levels(x)) {
            c_prob <- sum(weights[which(i == x)])
            probs  <- c(probs, c_prob)
        }

        rst <- list(type   = "discrete",
                    values = levels(x),
                    probs  = probs)
    }

    f_continous <- function(x, quants, weights) {
        x    <- as.numeric(x)
        m_sd <- f_msd(x, weights)
        rst  <- list(type   = "continuous",
                    range  = range(x),
                    mean   = m_sd[1],
                    ex2    = m_sd[2],
                    sd     = m_sd[3],
                    quants = f_quantile(x, quants, weights))
    }

    if (is.null(x))
        return(NULL)

    if (is.null(weights))
        weights <- rep(1, length(x))

    ## standardize weights
    weights <- weights / sum(weights)

    if (is.null(type)) {
        if (is.factor(x)) {
            rst <- f_factor(x, xlev, weights)
        } else {
            rst <- f_continous(x, quants, weights)
        }
    } else {
        rst <- switch(type,
                      "discrete" = f_factor(x, xlev, weights),
                      {
                          rst      <- f_continous(x, quants, weights)
                          rst$type <- type
                          rst
                      })
    }

    rst
}
