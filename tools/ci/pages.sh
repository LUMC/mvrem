#!/usr/bin/env bash
# The web page of the package: the vignette "Getting started with mvrem" (vignettes/mvrem.Rmd),
# rendered with the package installed from this commit, as public/index.html (GitLab Pages /
# GitHub Pages). The HTML is self-contained (KaTeX, widgets and figures embedded).
set -euo pipefail
cd "$(dirname "$0")/../.."
R CMD INSTALL --no-test-load .
mkdir -p public
Rscript -e 'rmarkdown::render("vignettes/mvrem.Rmd", output_file = "index.html", output_dir = "public")'
