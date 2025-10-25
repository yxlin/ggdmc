# q(save = "no")
cat("\n\n----------------- Testing model ecpe --------------------\n")
rm(list = ls())
pkg <- c("ggdmc")
pkg_ok <- sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
# data_path <- file.path(home_dir, "Group9_gen_cdm/data/dina_ecpe.rda")
# save_path <- file.path(home_dir, "Group10_cdm_subjects/data/dina_ecpe.rda")
source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))

#############################################################################
## EXAMPLE 1: simulate DINA/DINO data according to a tetrachoric correlation
#############################################################################

# define Q-matrix for 4 items and 2 attributes
q.matrix <- matrix(c(1, 0, 0, 1, 1, 1, 1, 1), ncol = 2, nrow = 4)
q.matrix

# Slipping parameters
slip <- c(0.2, 0.3, 0.4, 0.3)

# Guessing parameters
guess <- c(0, 0.1, 0.05, 0.2)

set.seed(1567) # fix random numbers
dat1 <- CDM::sim.din(
    N = 200, q.matrix, slip = slip, guess = guess,
    # Possession of the attributes with high probability
    mean = c(0.5, 0.2),
    # Possession of the attributes is weakly correlated
    Sigma = matrix(c(1, 0.2, 0.2, 1), ncol = 2), rule = "DINA"
)
names(dat1)
tibble::as_tibble(dat1$dat)
head(dat1$alpha)

set.seed(15367) # fix random numbers
res <- CDM::sim.din(
    N = 200, q.matrix, slip = slip, guess = guess, mean = c(0.5, 0.2),
    Sigma = matrix(c(1, 0.2, 0.2, 1), ncol = 2), rule = "DINO"
)

# extract simulated data
dat2 <- res$dat
# extract attribute patterns
head(res$alpha)
##        [,1] [,2]
##   [1,]    1    1
##   [2,]    1    1
##   [3,]    1    1
##   [4,]    1    1
##   [5,]    1    1
##   [6,]    1    0

#  simulate data based on given attributes
#          -> 5 persons with 2 attributes -> see the Q-matrix above
alpha <- matrix(c(1, 0, 1, 0, 1, 1, 0, 1, 1, 1),
    nrow = 5, ncol = 2, byrow = TRUE
)
# alpha
CDM::sim.din(q.matrix = q.matrix, alpha = alpha)

## Not run:
#############################################################################
# EXAMPLE 2: Simulation based on attribute vectors
#############################################################################
set.seed(76)
# define Q-matrix
Qmatrix <- matrix(c(1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1), 8, 2, byrow = TRUE)
colnames(Qmatrix) <- c("Attr1", "Attr2")
# define skill patterns
alpha.patt <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), 4, 2, byrow = TRUE)
AP <- nrow(alpha.patt)
# define pattern probabilities
alpha.prob <- c(.20, .40, .10, .30)
# simulate alpha latent responses
N <- 1000 # number of persons
ind <- sample(x = 1:AP, size = N, replace = TRUE, prob = alpha.prob)
alpha <- alpha.patt[ind, ] # (true) latent responses
# define guessing and slipping parameters
guess <- c(.26, .3, .07, .23, .24, .34, .05, .1)
slip <- c(.05, .16, .19, .03, .03, .19, .15, .05)
# simulation of the DINA model
dat <- CDM::sim.din(
    N = 0, q.matrix = Qmatrix, guess = guess,
    slip = slip, alpha = alpha
)$dat
# estimate model
res <- CDM::din(dat, q.matrix = Qmatrix)
# extract maximum likelihood estimates for individual classifications
est <- paste(res$pattern$mle.est)
# calculate classification accuracy
mean(est == apply(alpha, 1, FUN = function(ll) {
    paste0(ll[1], ll[2])
}))
##   [1] 0.935

#############################################################################
# EXAMPLE 3: Simulation based on already estimated DINA model for data.ecpe
#############################################################################

dat <- CDM::data.ecpe$data
q.matrix <- CDM::data.ecpe$q.matrix

#***
# (1) estimate DINA model
mod <- CDM::din(data = dat[, -1], q.matrix = q.matrix, rule = "DINA")

#***
# (2) simulate data according to DINA model
set.seed(977)
# number of subjects to be simulated
n <- 3000
# simulate attribute patterns
probs <- mod$attribute.patt$class.prob # probabilities
patt <- mod$attribute.patt.splitted # response patterns
alpha <- patt[sample(1:(length(probs)), n, prob = probs, replace = TRUE), ]
# simulate data using estimated item parameters
res <- CDM::sim.din(
    N = n, q.matrix = q.matrix, guess = mod$guess$est, slip = mod$slip$est,
    rule = "DINA", alpha = alpha
)
# extract data
dat <- res$dat

## End(Not run)
