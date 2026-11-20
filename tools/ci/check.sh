#!/usr/bin/env bash
# R CMD check of the package (tests are run by tools/ci/test.sh), with the vignette built and its
# code run (needs INLA). Fails on errors and warnings.
set -euo pipefail
cd "$(dirname "$0")/../.."
Rscript -e 'rcmdcheck::rcmdcheck(args = c("--no-manual", "--no-tests"), error_on = "warning", check_dir = "check")'
