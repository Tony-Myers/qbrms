## Resubmission
This is a resubmission (1.0.1) to address feedback from Uwe Ligges.
* Removed Bioconductor from 'Additional_repositories' field in DESCRIPTION, as requested.

## Test environments
* Local OS X install, R 4.3.2
* GitHub Actions (via rhub):
  * Ubuntu Linux 20.04 (R-release, R-devel)
  * Windows Server 2022 (R-release)
  * macOS (R-release)

## R CMD check results
0 errors | 0 warnings | 1 NOTE

## Explanation of NOTE
* **Suggests or Enhances not in mainstream repositories (INLA):**
  INLA is available from https://inla.r-inla-download.org/R/stable/, which is specified in the 'Additional_repositories' field.
* **New submission:** This is the first submission of this package.
* **Possibly misspelled words:**
  These are surnames of researchers cited in the Description field. They have been added to inst/WORDLIST.

## Reverse dependencies
This is a new package, so there are no reverse dependencies.