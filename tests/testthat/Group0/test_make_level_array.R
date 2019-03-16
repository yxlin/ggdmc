cat("\n-------------------- Testing make_level_array --------------------")

rm(list = ls())
## Create an array of all factorial combinations of factor levels
##
## Take a list storing one or multiple string vectors and glues the strings in
## these vectors with a dot symbol.
##
## @param x a list storing one or multiple factor vectors. Each vector
## can have different numbers of level.  For example,
## \code{f <- list(E = c("nc", "wc"), S = c("n", "w", "pn", "pw"))}
## @return a table showing the combinations of factor levels.
##
## @examples
## ## Example 1
## factors <- list(S = c("s1", "s2"))
## make_level_array(factors)
## ##   S
## ##  s1   s2
## ## "s1" "s2"
##
## factors <- list(S = c("left", "right"))
## make_level_array(factors)
## ##      S
## ##   left   right
## ##  "left" "right"
##
## ## Example 2
## factors <- list(A = c("00", "45", "90", "135", "180"),
##                 S = c("mirror", "normal"))
## make_level_array(factors)
## ##       S
## ## ------------------------------40
## ##   A        mirror       normal
## ## ------------------------------40
## ##  00   "00.mirror"  "00.normal"
## ##  45   "45.mirror"  "45.normal"
## ##  90   "90.mirror"  "90.normal"
## ## 135  "135.mirror" "135.normal"
## ## 180  "180.mirror" "180.normal"
##
## ## Example 3
## factors <- list(E = c("nc", "wc"),
##                 S = c("n", "w", "pn", "pw"))
## make_level_array(factors)
##
## ##         S
## ## -------------------------------40
## ## E       n      w      pn      pw
## ## -------------------------------40
## ## nc "nc.n" "nc.w" "nc.pn" "nc.pw"
## ## wc "wc.n" "wc.w" "wc.pn" "wc.pw"
##


## Example 1
factors <- list(S = c("s1", "s2"))
ggdmc:::make_level_array(factors)
##   S
##  s1   s2
## "s1" "s2"

factors <- list(S = c("left", "right"))
ggdmc:::make_level_array(factors)
##      S
##   left   right
##  "left" "right"

## Example 2
factors <- list(A = c("00", "45", "90", "135", "180"),
                S = c("mirror", "normal"))
ggdmc:::make_level_array(factors)
##       S
## ------------------------------40
##   A        mirror       normal
## ------------------------------40
##  00   "00.mirror"  "00.normal"
##  45   "45.mirror"  "45.normal"
##  90   "90.mirror"  "90.normal"
## 135  "135.mirror" "135.normal"
## 180  "180.mirror" "180.normal"

## Example 3
factors <- list(E = c("nc", "wc"),
                S = c("n", "w", "pn", "pw"))
ggdmc:::make_level_array(factors)

##         S
## -------------------------------40
## E       n      w      pn      pw
## -------------------------------40
## nc "nc.n" "nc.w" "nc.pn" "nc.pw"
## wc "wc.n" "wc.w" "wc.pn" "wc.pw"
