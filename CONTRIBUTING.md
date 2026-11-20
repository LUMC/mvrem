# Contributing to mvrem

Thank you for helping to improve **mvrem**.

If anything in this guide does not work for you, open an issue (step 1 below) and say so. That
is a useful contribution too.

---

## Contents

1. [What you can contribute](#1-what-you-can-contribute)
2. [Words used in this guide](#2-words-used-in-this-guide)
3. [One-time setup](#3-one-time-setup)
4. [Making a contribution, step by step](#4-making-a-contribution-step-by-step)
5. [What a good contribution contains](#5-what-a-good-contribution-contains)
6. [After you opened the pull request](#6-after-you-opened-the-pull-request)
7. [Small changes without installing anything](#7-small-changes-without-installing-anything)
8. [When something goes wrong](#8-when-something-goes-wrong)
9. [Command summary](#9-command-summary)

---

## 1. What you can contribute

- **A report**: a wrong estimate, an error message you did not understand, a model mvrem
  should handle, a sentence in the documentation that is unclear. Only a GitHub account is
  needed. See "Reporting a problem or an idea" below.
- **An example**: simulated or public data with a model written in mvrem notation, ideally with
  the true values or published estimates to compare with.
- **A statistical check**: a simulation that shows whether an estimate is biased or an interval
  covers.
- **A change to the code or documentation**: a fix, a new family, a new prior, a better
  default.

**Before you write code, open an issue first.** Describe the change in a few sentences.
The maintainers can then tell you whether it fits the package, and whether someone is already
working on it. This saves you work.

### Reporting a problem or an idea

1. Go to <https://github.com/LUMC/mvrem/issues> and click **New issue**.
2. Give it a short title, for example "Beta outcome: gamma attenuated with 20 subjects".
3. In the text, write:
   - what you did: the R code, as short as possible, with simulated data if you can (use
     `set.seed()` so that others get the same numbers);
   - what you expected, and what happened instead (copy the output or the error message);
   - the output of `sessionInfo()` and `packageVersion("INLA")`.
4. Put R code between lines of three backticks, so that GitHub shows it as code:

   ````
   ```r
   library(mvrem)
   set.seed(1)
   ...
   ```
   ````

5. Click **Submit new issue**.

---

## 2. Words used in this guide

| Word | Meaning |
|---|---|
| **repository** ("repo") | The project folder with its whole history. The official one is `LUMC/mvrem` on GitHub. |
| **fork** | Your own copy of the repository on GitHub (`<your-name>/mvrem`). You may change it freely. |
| **clone** | A copy of a repository on your own computer. |
| **commit** | A saved step in the history: a set of changes with a short message saying what and why. |
| **branch** | A separate line of work. Each contribution gets its own branch, so that unrelated changes do not mix. |
| **`main`** | The branch of mvrem on GitHub: one commit per release. All contributions start from it and go back into it. |
| **push** | Send your commits from your computer to your fork on GitHub. |
| **pull request** ("PR") | A request to the maintainers to take the commits of your branch into `main`. It is where the change is discussed and reviewed. |
| **upstream** | The official repository `LUMC/mvrem`, as seen from your clone. |
| **origin** | Your fork, as seen from your clone. |

The picture to keep in mind:

```
 LUMC/mvrem (upstream)  --- fork --->  <you>/mvrem (origin)  --- clone --->  your computer
        ^                                     ^                                   |
        |                                     +------------- push ----------------+
        +------------- pull request (reviewed, then merged by a maintainer) -------+
```

---

## 3. One-time setup

You do this once per computer. It takes about an hour, mostly waiting for installations.

### 3.1 Create a GitHub account

Sign up at <https://github.com>. Use an email address you read; GitHub sends review comments
there.

### 3.2 Install the software

1. **R** (version 4.3 or newer): <https://cran.r-project.org>.
2. **RStudio Desktop**: <https://posit.co/download/rstudio-desktop/>. This guide uses
   RStudio, but any editor works.
3. **Git**:
   - Windows: <https://git-scm.com/download/win>. Keep all default options.
   - macOS: open the Terminal app and type `git --version`. If Git is missing, macOS offers
     to install it; accept.
   - Linux: install the `git` package of your distribution, e.g. `sudo apt install git`.
4. Restart RStudio, so that it finds Git.
5. **INLA and the R packages needed for development.** In the RStudio *Console*:

   ```r
   install.packages("INLA",
     repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
     dep = TRUE)
   install.packages(c("devtools", "usethis", "gitcreds", "lme4", "testthat",
                      "htmlwidgets", "jsonlite", "tibble", "rmarkdown"))
   ```

   Check that INLA works:

   ```r
   INLA::inla(y ~ 1, data = data.frame(y = rnorm(20)))$summary.fixed
   ```

   You should see a small table with the estimated intercept.

> **Alternative: the ready-made container.** The repository contains a development container
> (`.devcontainer/`) with R, INLA and everything else installed. If you know Docker, or use
> VS Code with the "Dev Containers" extension, you can skip step 5 and open the project in that
> container instead.

### 3.3 Tell Git who you are

Git writes your name and email into every commit. In the RStudio *Console*:

```r
usethis::use_git_config(user.name = "Jane Doe", user.email = "jane.doe@example.org")
```

Use the same email address as on GitHub.

### 3.4 Let your computer log in to GitHub

GitHub does not accept your password from R or Git; it needs a **token** (a long random
password that only your computer knows).

1. In the RStudio *Console*:

   ```r
   usethis::create_github_token()
   ```

   A GitHub page opens. Sign in if asked. Leave the pre-selected options, set an
   expiration date (for example 90 days), and click **Generate token** at the bottom of the
   page.
2. Copy the token that GitHub shows (it starts with `ghp_`). GitHub shows it only once.
3. Back in the *Console*:

   ```r
   gitcreds::gitcreds_set()
   ```

   Paste the token when asked and press Enter.
4. Check:

   ```r
   usethis::gh_token_help()
   ```

   It should say that a token is found. When the token expires, repeat these steps.

### 3.5 Fork and clone mvrem

One command creates your fork on GitHub, copies it to your computer and opens it as an
RStudio project:

```r
usethis::create_from_github("LUMC/mvrem", fork = TRUE, destdir = "~/projects")
```

(Choose any folder for `destdir`.) A new RStudio window opens with the project. In this window,
open the *Terminal* tab (next to *Console*; if it is missing: menu *Tools > Terminal > New
Terminal*) and type:

```sh
git remote -v
```

You should see four lines: `origin` pointing to `github.com/<your-name>/mvrem` and `upstream`
pointing to `github.com/LUMC/mvrem`. That is exactly the picture of section 2.

Finally, make sure your `main` follows the official one:

```sh
git fetch upstream
git checkout main
git branch -u upstream/main
```

### 3.6 Check that everything works

In the RStudio *Console* of the project:

```r
devtools::load_all()   # loads the package from your copy of the code
devtools::test()       # runs all tests: takes about 5 minutes
```

At the end, the test summary should say `[ FAIL 0 | ... ]`. One test is skipped on purpose
(it marks a known bug). If tests fail on a fresh copy, please open an issue with the output.

The setup is done.

---

## 4. Making a contribution, step by step

Repeat these steps for every contribution. Commands starting with `git` are typed in the
RStudio **Terminal** tab; R commands in the **Console** tab.

### Step 1: Agree on the change in an issue

See [section 1](#1-what-you-can-contribute). Note the issue number (e.g. `#42`).

### Step 2: Bring your copy up to date

Others may have changed `main` since you last worked on it:

```sh
git checkout main
git pull upstream main
```

The first line switches to the development branch; the second fetches the newest version of it
from `LUMC/mvrem`.

### Step 3: Create a branch for your change

Choose a short name without spaces that says what the change is, e.g. `fix-beta-gamma` or
`example-pbc`:

```sh
git checkout -b fix-beta-gamma
```

You are now on your own branch. Nothing you do here affects `main` until a maintainer
merges your pull request.

To see which branch you are on at any time: `git status` (first line), or the branch name at
the top right of RStudio's *Git* tab.

### Step 4: Make the change

Edit the files in RStudio as you would edit any R script. Where things are:

| Folder / file | Contains |
|---|---|
| `R/` | the code of the package, one file per topic (`families.R`, `priors.R`, `plot.R`, ...) |
| `tests/testthat/` | the tests, one file per topic (`test-families-priors.R`, `test-fit.R`, ...) |
| `man/` | help pages: **generated**, do not edit by hand (see below) |
| `tools/` | simulation scripts used to check the statistical behaviour, the CI scripts (`ci/`) |
| `vignettes/` | the vignette "Getting started with mvrem" (also the web page of the package) and the KaTeX copy it and the LaTeX tests use (`katex/`) |
| `CLAUDE.md` | a detailed technical description of the package's design and decisions |

While working, after each change of the code:

```r
devtools::load_all()
```

loads your version of the package, so you can try it out right away in the *Console* (no need
to install it).

The help pages are written as comments above each function (lines starting with `#'`,
"roxygen"). After changing such comments, regenerate the help pages:

```r
devtools::document()
```

### Step 5: Add or update tests

A test is a small piece of R code that checks one behaviour and complains if it changes. Tests
protect your contribution from being broken by later changes. Add them to the test file of the
topic, in the same style as the tests around them, for example:

```r
test_that("betaFamily(): rejects a response outside (0, 1)", {
  d <- data.frame(id = rep(1:5, each = 2), t = 0:1, x = rnorm(10), y = 1.5)
  expect_error(jsm(list(y ~ 1 + scaled(x), x ~ 1 + t + (1 | id)), data = d,
                   family = list(y = betaFamily()), fit = FALSE),
               "must be in \\(0, 1\\)")
})
```

For a statistical change, a test can fit a small simulated model and check that the estimate
is near the true value (with a generous tolerance: INLA's results differ slightly between
machines and versions).

Run the tests of one file while working (faster), and all tests before you continue:

```r
devtools::test_active_file()   # the test file open in the editor
devtools::test()               # everything (about 5 minutes)
```

### Step 6: Commit

A commit saves your changes in the history. Make one commit per meaningful step, e.g. "the
fix" and "the test", or one for both.

**With RStudio's Git tab** (top right pane):

1. Tick the *Staged* box of each file that belongs to this step.
2. Click **Commit**. A window shows the changes (red = removed, green = added): check them.
3. Write a message: a first line of at most about 70 characters saying *what* changed, then an
   empty line, then *why*, if it is not obvious. Example:

   ```
   betaFamily(): check that the response lies in (0, 1)

   A response of exactly 0 or 1 made INLA fail with an unclear message (#42).
   ```

4. Click **Commit**, then close the window.

**Or in the Terminal:**

```sh
git status                      # which files changed
git diff                        # what changed in them
git add R/families.R tests/testthat/test-families-priors.R
git commit -m "betaFamily(): check that the response lies in (0, 1)"
```

Do not commit large data files, personal files, or files with passwords.

### Step 7: Push your branch to your fork

The first time for a branch:

```sh
git push -u origin fix-beta-gamma
```

After that, `git push` (or RStudio's green **Push** arrow) is enough.

### Step 8: Open the pull request

1. Go to your fork, `https://github.com/<your-name>/mvrem`. GitHub usually shows a yellow bar
   "fix-beta-gamma had recent pushes" with a button **Compare & pull request**: click it. (If
   not: click **Pull requests > New pull request**.)
2. **Important: check the two branches at the top.** It must read

   ```
   base repository: LUMC/mvrem   base: main   <-   head repository: <your-name>/mvrem   compare: fix-beta-gamma
   ```

   GitHub suggests exactly this by default.
3. Title: the same as your main commit message.
4. Description: explain
   - what the change does and why (refer to the issue: `Fixes #42`; GitHub then closes the issue
     when the pull request is merged);
   - for a statistical change: the reasoning or reference, and the evidence (simulation
     results with the true values, number of replicates, bias / coverage);
   - how you tested it (`devtools::test()` passed; which new tests);
   - anything you are unsure about. Questions are welcome in a pull request.
5. Click **Create pull request**.

If the work is not finished but you would like early feedback, choose **Create draft pull
request** (arrow next to the button) instead.

---

## 5. What a good contribution contains

- **One topic per pull request.** Two unrelated fixes = two branches and two pull requests.
- **Tests** for the new or changed behaviour (step 5); all tests pass.
- **Documentation**: roxygen comments for new or changed arguments and return values, then
  `devtools::document()`. For a user-visible feature, a short example in the help page
  (roxygen `@examples`) or a description in the pull request is very welcome; the maintainers
  add it to the worked examples.
- **Statistical evidence** for anything that changes estimates: a simulation with known true
  values (several seeds), as done for the existing priors (see `tools/simDataPriors.R` and the
  "Evidence" notes in `CLAUDE.md`). Report bias and interval coverage, and the settings.
- **Code style**: follow the code around your change. Names are in `camelCase`
  (`fittedValues`, not `fitted_values`); indent with two spaces; keep functions short. Error
  messages say what is wrong and what to do (`"... Use na.action = na.omit to drop such rows."`).
- **No unrelated changes**: do not reformat files or rename things outside your topic; it makes
  the review hard.

You are not expected to get all of this right the first time. The maintainers will help in the
review.

---

## 6. After you opened the pull request

- A maintainer reviews the change and may comment on specific lines, ask questions or request
  changes. GitHub emails you about every comment.
- **To change something:** stay on the same branch on your computer, edit, commit and push
  again (steps 4-7). The pull request updates itself; do not open a new one.

  ```sh
  git checkout fix-beta-gamma
  # ... edit, test ...
  git add ...
  git commit -m "Address review: clearer error message"
  git push
  ```

- Answer each comment (a short reply such as "done" is fine), so that the reviewer knows it was
  handled.
- **If `main` moved on** and the maintainer asks you to update your branch:

  ```sh
  git checkout fix-beta-gamma
  git pull upstream main
  git push
  ```

  If Git reports a *conflict*, see [section 8](#8-when-something-goes-wrong).
- When the reviewer is satisfied, a maintainer **merges** the pull request. Your change is now
  part of mvrem.
- **Clean up** (optional): switch back and delete the branch:

  ```sh
  git checkout main
  git pull upstream main
  git branch -d fix-beta-gamma
  git push origin --delete fix-beta-gamma
  ```

---

## 7. Small changes without installing anything

For a typo or a sentence in the documentation (for example in `README.md` or a help page's roxygen comment),
GitHub's website is enough:

1. Open the file on GitHub, e.g. <https://github.com/LUMC/mvrem/blob/main/README.md>.
2. Click the pencil icon (**Edit this file**). GitHub creates a fork for you if you have none.
3. Make the change. Click **Commit changes...**, write a short message, choose
   **Create a new branch** and **Propose changes**.
4. On the next page, check that the base is `LUMC/mvrem` / `main` and click **Create pull request**.

Do not use this for code changes: the tests cannot run this way.

---

## 8. When something goes wrong

**"I made changes on `main` instead of on a new branch."**
Before committing: just create the branch now, the changes come with you:
`git checkout -b my-change`. After committing, ask for help in your issue; nothing is lost.

**"`git pull` says: Your local changes would be overwritten."**
You have uncommitted changes. Commit them first (step 6), or set them aside with `git stash`,
pull, and bring them back with `git stash pop`.

**"Merge conflict" / "CONFLICT (content)".**
You and someone else changed the same lines. Git marks them in the file like this:

```
<<<<<<< HEAD
your version of the lines
=======
the other version of the lines
>>>>>>> upstream/main
```

Open each file Git names (`git status` lists them as "both modified"), decide what the lines
should be, delete the three marker lines, and save. Then:

```sh
git add <the files>
git commit -m "Merge main"
git push
```

If you are unsure, run `git merge --abort` to go back to before the pull, and ask in the pull
request.

**"Push rejected" / "Authentication failed".**
Your token has expired or is missing: repeat [section 3.4](#34-let-your-computer-log-in-to-github).
If the message says *non-fast-forward*, run `git pull` first, then push again.

**"The tests fail, but not because of my change."**
Run `devtools::test()` on an up-to-date `main` without your change (`git checkout main`).
If they fail there too, report it in your pull request; it is not your problem to fix.

**"I am lost."**
`git status` always tells you which branch you are on and what has changed. Copy its output into
a comment in your issue or pull request and ask. Nothing you do in your fork or branch can
damage the official repository.

---

## 9. Command summary

| When | Where | Command |
|---|---|---|
| start a contribution | Terminal | `git checkout main` |
| | Terminal | `git pull upstream main` |
| | Terminal | `git checkout -b my-change` |
| while working | Console | `devtools::load_all()` |
| | Console | `devtools::document()` (after changing `#'` comments) |
| | Console | `devtools::test()` |
| save a step | Terminal | `git add <files>` then `git commit -m "message"` (or the Git tab) |
| publish | Terminal | `git push -u origin my-change` (first time), then `git push` |
| propose | GitHub | **Compare & pull request**, base = `LUMC/mvrem` / **`main`** |
| update after review | Terminal | edit, commit, `git push` |
| where am I? | Terminal | `git status` |

Thank you for contributing!
