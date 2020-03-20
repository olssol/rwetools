
#' Run Web-Based \code{rwetools} application
#'
#' Call Shiny to run \code{rwetools} as a web-based application. A web browser
#' will be brought up.
#'
#' @examples
#' \dontrun{
#' rweShiny()}
#'
#'
#' @export
#'
rwe_shiny <- function() {

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

    app_dir <- system.file("shiny", package = "etpi")
    if (app_dir == "") {
        stop("Could not find Shiny directory. Try re-installing `rwetools`.",
             call. = FALSE)
    }


    shiny::runApp(app_dir, display.mode = "normal");
}
