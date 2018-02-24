
#' Run Web-Based \code{rwetools} application
#'
#' Call Shiny to run \code{rwetools} as a web-based application. A web browser will
#' be brought up.
#'
#' @examples
#' \dontrun{
#' rweShiny()}
#'
#'
#' @export
#'
rweShiny <- function() {

    if (!requireNamespace("shiny", quietly = TRUE)) {
        stop("Shiny needed for this function to work. Please install it.",
             call. = FALSE)
    }

    if (!requireNamespace("shinythemes", quietly = TRUE)) {
        stop("shinythemes needed for this function to work. Please install it.",
             call. = FALSE)
    }

    if (!requireNamespace("DT", quietly = TRUE)) {
        stop("DT needed for this function to work. Please install it.",
             call. = FALSE)
    }

    appDir <- system.file("shiny", package = "etpi")
    if (appDir == "") {
        stop("Could not find Shiny directory. Try re-installing `rwetools`.",
             call. = FALSE)
    }


    shiny::runApp(appDir, display.mode = "normal");
}
