
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
tkExtractStats <- function(x, quants = c(0.5), xlev = NULL, type = NULL) {

    f_factor <- function(x, xlev) {
        x <- as.factor(x)
        if (!is.null(xlev))
            levels(x) <- xlev

        probs <- as.numeric(table(x))
        probs <- probs / sum(probs)

        rst <- list(type   = "discrete",
                    values = levels(x),
                    probs  = probs)
    }

    f_continous <- function(x, quants) {
        x   <- as.numeric(x)
        rst <- list(type   = "continuous",
                    range  = range(x),
                    mean   = mean(x),
                    sd     = sd(x),
                    quants = cbind(quants, quantile(x, quants)))
    }

    if (is.null(x))
        return(NULL)

    if (is.null(type)) {
        if (is.factor(x)) {
            rst <- f_factor(x, xlev)
        } else {
            rst <- f_continous(x, quants)
        }
    } else {
        rst <- switch(type,
                      "discrete" = f_factor(x, xlev),
                      {
                          rst      <- f_continous(x, quants)
                          rst$type <- type
                          rst
                      })
    }

    rst
}

