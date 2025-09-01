pkgname <- "qbrms"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "qbrms-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('qbrms')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("qbrms")
### * qbrms

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: qbrms
### Title: Main q-brms model fitting interface (simplified example)
### Aliases: qbrms

### ** Examples

## Not run: 
##D # Simple linear regression
##D fit1 <- qbrms(y ~ x, data = data, family = gaussian())
##D 
##D # Mixed effects model
##D fit2 <- qbrms(y ~ x + (1 | group), data = data, family = gaussian())
##D 
##D # Quantile regression
##D fit3 <- qbrms(y ~ x, data = data, family = asymmetric_laplace(), quantile = 0.9)
##D 
##D # Prior predictive check
##D fit_prior <- qbrms(y ~ x, data = data, sample_prior = "only")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("qbrms", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("qmbs")
### * qmbs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: qmbs
### Title: Fit Bayesian Models using qbrms (Alternative Interface)
### Aliases: qmbs

### ** Examples

## Not run: 
##D # Simple linear regression
##D fit <- qmbs(y ~ x, data = data, family = gaussian())
##D 
##D # Prior predictive check
##D fit_prior <- qmbs(y ~ x, data = data, sample_prior = "only")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("qmbs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
