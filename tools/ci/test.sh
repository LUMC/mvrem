#!/usr/bin/env bash
# The package's tests (testthat), with a JUnit report for GitLab (test-results.xml).
# Fails when a test fails, or when INLA, htmlwidgets or node is missing (their tests would be skipped).
set -euo pipefail
cd "$(dirname "$0")/../.."
Rscript -e 'stopifnot(requireNamespace("INLA", quietly = TRUE), requireNamespace("htmlwidgets", quietly = TRUE),
                      nzchar(Sys.which("node")))'
# Install this version: the jsmApply(backend = "psock") tests start separate R sessions, which load
# the INSTALLED mvrem (they are skipped without one, and test a stale one silently).
log="${TMPDIR:-/tmp}/mvrem-install.log"
R CMD INSTALL --no-docs --no-multiarch . > "$log" 2>&1 || { cat "$log"; exit 1; }
Rscript -e 'rep <- testthat::MultiReporter$new(list(testthat::SummaryReporter$new(show_praise = FALSE),
                                                    testthat::JunitReporter$new(file = file.path(getwd(), "test-results.xml"))))
            devtools::test(reporter = rep, stop_on_failure = TRUE)'
