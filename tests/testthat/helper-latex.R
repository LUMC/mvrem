# Helpers for the LaTeX tests (test-latex.R); reusable by other files:
#
#   texLines(x)               lines of the aligned body of an mvremLatex (no \begin / \end)
#   texLine(lines, pattern)   the ONE line matching `pattern` (fixed), else a test failure
#   bracesBalanced(x)         TRUE if { } are balanced (escaped \{ \} and \\ ignored)
#   expectWellFormed(x)       braces and \left / \right balanced
#   unwrapCommand(x, cmd)     strip every `cmd{...}{sym}` wrapper (\htmlData, \underset), keeping sym
#   htmlDataWrappers(x)       data frame of the \htmlData{arg}{sym} wrappers: `arg`, `sym`
#   expectValidHtmlData(arg)  a \htmlData key=value string is well formed
#   skipWithoutKatex()        skip when node or the vendored KaTeX is missing
#   katexCheck(texs)          "OK" or the KaTeX error per string (ONE node call)
#   famTex(f, fam)            latexEquation() of one parsed formula with a family attached
#   withEstState(est, expr)   run expr with latexState$est set (restored afterwards)

texLines <- function(x) {
  body <- sub("^\\\\begin\\{aligned\\}\n", "", sub("\n\\\\end\\{aligned\\}$", "", unclass(x)))
  strsplit(body, " \\\\\n", fixed = TRUE)[[1]]
}

texLine <- function(lines, pattern) {
  hit <- grep(pattern, lines, fixed = TRUE)
  testthat::expect_length(hit, 1L)
  if (length(hit) == 1L) lines[hit] else NA_character_
}

bracesBalanced <- function(x) {
  x <- gsub(r"(\\\\)", "", unclass(x))
  x <- gsub(r"(\\[{}])", "", x)
  ch <- strsplit(x, "")[[1]]
  depth <- cumsum((ch == "{") - (ch == "}"))
  !length(depth) || (all(depth >= 0) && depth[length(depth)] == 0)
}

leftRightBalanced <- function(x) {
  count <- function(p) lengths(regmatches(x, gregexpr(p, x, perl = TRUE)))
  count(r"(\\left(?![a-zA-Z]))") == count(r"(\\right(?![a-zA-Z]))")
}

expectWellFormed <- function(x) {
  testthat::expect_true(bracesBalanced(x), info = unclass(x))
  testthat::expect_true(leftRightBalanced(unclass(x)), info = unclass(x))
}

# Index of the brace matching the "{" at position i of the character vector ch.
matchBrace <- function(ch, i) {
  depth <- 0L
  j <- i
  while (j <= length(ch)) {
    if (ch[j] == "\\") {
      j <- j + 2L
      next
    }
    if (ch[j] == "{") depth <- depth + 1L
    if (ch[j] == "}") {
      depth <- depth - 1L
      if (depth == 0L) return(j)
    }
    j <- j + 1L
  }
  stop("unbalanced braces")
}

# Remove every `cmd{drop args}{sym}` wrapper (innermost last), keeping `sym`.
unwrapCommand <- function(x, cmd, drop = 1L) {
  x <- unclass(x)
  repeat {
    pos <- gregexpr(cmd, x, fixed = TRUE)[[1]]
    if (pos[1] < 0) return(x)
    p <- pos[length(pos)]
    ch <- strsplit(x, "")[[1]]
    i <- p + nchar(cmd)
    for (k in seq_len(drop)) {
      stopifnot(ch[i] == "{")
      i <- matchBrace(ch, i) + 1L
    }
    stopifnot(ch[i] == "{")
    j <- matchBrace(ch, i)
    x <- paste0(substr(x, 1L, p - 1L), substr(x, i + 1L, j - 1L), substring(x, j + 1L))
  }
}

# The \htmlData{arg}{sym} wrappers of x, outermost first: `arg` (the key=value string)
# and `sym` (the wrapped LaTeX, which may itself contain wrappers).
htmlDataWrappers <- function(x) {
  x <- unclass(x)
  ch <- strsplit(x, "")[[1]]
  cmd <- "\\htmlData"
  pos <- gregexpr(cmd, x, fixed = TRUE)[[1]]
  if (pos[1] < 0) return(data.frame(arg = character(), sym = character(), stringsAsFactors = FALSE))
  res <- lapply(pos, function(p) {
    i <- p + nchar(cmd)
    j <- matchBrace(ch, i)
    k <- matchBrace(ch, j + 1L)
    c(arg = substr(x, i + 1L, j - 1L), sym = substr(x, j + 2L, k - 1L))
  })
  data.frame(arg = vapply(res, `[[`, "", "arg"), sym = vapply(res, `[[`, "", "sym"),
             stringsAsFactors = FALSE)
}

# Checks a \htmlData key=value string: `par` first, known keys only, no duplicated keys,
# values without = { }, numbers (or NA) in mean / sd / lo / hi (vectors joined by ";").
expectValidHtmlData <- function(arg) {
  testthat::expect_false(grepl("[{}]", arg), info = arg)
  pairs <- strsplit(arg, ",", fixed = TRUE)[[1]]
  testthat::expect_match(pairs[1], "^par=", info = arg)
  testthat::expect_true(all(grepl("^(par|item|mean|sd|lo|hi|band|note)=[^=]*$", pairs)), info = arg)
  keys <- sub("=.*$", "", pairs)
  testthat::expect_false(anyDuplicated(keys) > 0, info = arg)
  vals <- sub("^[^=]*=", "", pairs)
  for (v in vals[keys %in% c("mean", "sd", "lo", "hi")])
    testthat::expect_true(all(grepl("^(NA|-?[0-9.]+(e[-+][0-9]+)?)$", strsplit(v, ";", fixed = TRUE)[[1]])),
                          info = arg)
}

katexJs <- function() testthat::test_path("..", "..", "vignettes", "katex", "katex.min.js")

skipWithoutKatex <- function() {
  testthat::skip_if_not(file.exists(katexJs()), "vendored KaTeX not available")
  testthat::skip_if(Sys.which("node") == "", "node not available")
}

# KaTeX parse of each string (node + vendored KaTeX; trust for \htmlData): "OK" or the error,
# named as `texs`. Call skipWithoutKatex() first; node runs once for all strings.
katexCheck <- function(texs) {
  script <- tempfile(fileext = ".js")
  writeLines(c("const katex = require(process.argv[2]);",
               "const fs = require('fs');",
               "for (const f of process.argv.slice(3)) {",
               "  try {",
               "    katex.renderToString(fs.readFileSync(f, 'utf8'),",
               "      {throwOnError: true, trust: true, strict: 'ignore', displayMode: true});",
               "    console.log('OK');",
               "  } catch (e) { console.log('ERR ' + e.message.split('\\n')[0]); }",
               "}"), script)
  files <- vapply(texs, function(t) {
    f <- tempfile(fileext = ".tex")
    writeLines(unclass(t), f, useBytes = TRUE)
    f
  }, "")
  on.exit(unlink(c(script, files)), add = TRUE)
  out <- system2("node", c(shQuote(script), shQuote(normalizePath(katexJs())), shQuote(files)), stdout = TRUE)
  stats::setNames(as.character(out), names(texs)[seq_along(out)])
}

famTex <- function(f, fam, part = "out", data = NULL) {
  p <- parseJsmFormula(f, part = part)
  if (!missing(fam)) p$family <- fam
  latexEquation(p, data)
}

withEstState <- function(est, expr) {
  old <- latexState$est
  latexState$est <- est
  on.exit(latexState$est <- old, add = TRUE)
  force(expr)
}
