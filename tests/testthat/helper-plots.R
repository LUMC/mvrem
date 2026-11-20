# Helpers of the plot tests (test-plots.R): widget checks and a node runner of the
# widget's JavaScript (tests/testthat/js/).

hasWidgets <- function() requireNamespace("htmlwidgets", quietly = TRUE)
skipNoWidgets <- function() testthat::skip_if_not(hasWidgets(), "htmlwidgets not installed")

# number of non-finite numbers (NA, NaN, Inf) anywhere in a nested list
nonFinite <- function(x) {
  if (is.list(x)) return(sum(vapply(x, nonFinite, 0)))
  if (is.numeric(x)) return(sum(!is.finite(x)))
  0
}

# a widget whose spec holds only finite numbers and serialises to JSON without NA / Inf
expectCleanWidget <- function(w) {
  testthat::skip_if_not_installed("jsonlite")
  testthat::expect_s3_class(w, "mvremPlot")
  testthat::expect_s3_class(w, "htmlwidget")
  testthat::expect_equal(nonFinite(w$x), 0)
  j <- jsonlite::toJSON(w$x, auto_unbox = TRUE)
  testthat::expect_false(grepl("\"(NA|NaN|Inf|-Inf)\"", j))
  invisible(w)
}

# ids of the lines of a widget panel (of one `kind` if given), without duplicates
lineIds <- function(panel, kind = NULL) {
  ls <- panel$lines
  if (!is.null(kind)) ls <- Filter(function(l) l$kind == kind, ls)
  unique(vapply(ls, `[[`, "", "id"))
}

lineKinds <- function(panel) vapply(panel$lines, `[[`, "", "kind")

subjectIds <- function(panel) vapply(panel$subjects, `[[`, "", "id")

panelTitles <- function(w) vapply(w$x$panels, `[[`, "", "title")

legendKinds <- function(w) vapply(w$x$legend, `[[`, "", "kind")

# all subject ids of the simulated data (sorted, as the plots pick them)
allIds <- function() sort(unique(fxData()$long$id))

# the text size plotModel() chose on a device of the given size (inches)
plotModelCex <- function(m, width, height, ...) {
  grDevices::pdf(NULL, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  plotModel(m, ...)
}

# ---- node runner of inst/htmlwidgets/mvremPlot.js --------------------------------------------

hasNode <- function() nzchar(Sys.which("node"))

# Render the specs of the widgets (a named list) with the widget's JavaScript in node, on the
# DOM stub of tests/testthat/js/domStub.js, and exercise the search box, the legend and the
# mouse / keyboard handlers (tests/testthat/js/renderSpecs.js). Returns one result per widget
# (parsed JSON: `ok`, `error`, element counts and the states after the interactions).
runWidgetJs <- function(widgets) {
  testthat::skip_if_not_installed("jsonlite")
  testthat::skip_if_not(hasNode(), "node not available")
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  files <- vapply(names(widgets), function(nm) {
    f <- file.path(dir, paste0(nm, ".json"))
    # as htmlwidgets serialises the spec: scalars unboxed, NULL / NA as null
    jsonlite::write_json(widgets[[nm]]$x, f, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
    f
  }, "")
  js <- system.file("htmlwidgets", "mvremPlot.js", package = "mvrem")
  driver <- testthat::test_path("js", "renderSpecs.js")
  out <- suppressWarnings(system2("node", shQuote(c(driver, js, files)), stdout = TRUE, stderr = TRUE))
  res <- tryCatch(jsonlite::fromJSON(paste(out, collapse = "\n"), simplifyVector = FALSE),
                  error = function(e) stop("node output is not JSON:\n", paste(out, collapse = "\n"), call. = FALSE))
  stats::setNames(res, names(widgets))
}
