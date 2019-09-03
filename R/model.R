######### CORE ---------------------------------------------------
##' Bayeisan computation of response time models
##'
##' \pkg{ggdmc} uses the population-based Markov chain Monte Carlo to
##' conduct Bayesian computation on cognitive models.
##'
##' @keywords package
##' @name ggdmc
##' @docType package
##' @author  Yi-Shin Lin <yishinlin001@gmail.com> \cr
##' Andrew Heathcote <andrew.heathcote@utas.edu.au>
##' @references
##' Heathcote, A., Lin, Y.-S., Reynolds, A., Strickland, L., Gretton, M. &
##' Matzke, D., (2018). Dynamic model of choice.
##' \emph{Behavior Research Methods}.
##' https://doi.org/10.3758/s13428-018-1067-y. \cr
##'
##' Turner, B. M., & Sederberg P. B. (2012). Approximate Bayesian computation
##' with differential evolution, \emph{Journal of Mathematical Psychology}, 56,
##' 375--385. \cr
##'
##' Ter Braak (2006). A Markov Chain Monte Carlo version of the genetic
##' algorithm Differential Evolution: easy Bayesian computing for real
##' parameter spaces. \emph{Statistics and Computing}, 16, 239-249.
##'
##' @importFrom Rcpp evalCpp
##' @useDynLib ggdmc
NULL

######### Model Generic -----------------------------------

grepl_exact <- function(xdot, pattern) {
  xvec <- strsplit(xdot, ".", fixed = TRUE)[[1]]
  any(pattern == xvec)
}

grepl_dot <- function(pattern, x) {
  ## Splits up pattern at ".", to get n parts. Matches each part to x and
  ## returens TRUE for each x where exactly n matches in any order.

  ps <- strsplit(pattern, ".", fixed = TRUE)[[1]]
  out <- sapply(x, grepl_exact, pattern = ps[1])
  if (length(ps)>1) {
    for (i in ps[-1])
      out <- rbind(out,sapply(x, grepl_exact, pattern = i))
    apply(out, 2, sum) == length(ps)
  } else out
}

##' Create a model object
##'
##' A model object consists of arraies with model attributes.
##'
##' @param p.map parameter map. This option maps a particular factorial design
##' to model parameters
##' @param match.map match map. This option matches stimuli and responses
##' @param factors specifying a list of factors and their levels
##' @param constants specifying the parameters with fixed values
##' @param responses specifying the response names and levels
##' @param type specifying model type, either "rd" or "norm".
##' @param posdrift a Boolean, switching between enforcing strict postive drift
##' rates by using truncated normal distribution. This option is only useful in
##' "norm" model type.
##' @param verbose Print p.vector, constants and model type
##' @param x a model object
##' @param p.vector parameter vector
##' @param ... other arguments
##' @importFrom utils glob2rx
##' @export
##' @examples
##' model <- BuildModel(
##'         p.map     = list(a = "1", v = "1", z = "1", d = "1", t0 = "1",
##'                      sv = "1", sz = "1", st0 = "1"),
##'         constants = c(st0 = 0, d = 0, sz = 0, sv = 0),
##'         match.map = list(M = list(s1 = "r1", s2 = "r2")),
##'         factors   = list(S = c("s1", "s2")),
##'         responses = c("r1", "r2"),
##'         type      = "rd")
BuildModel <- function(
  p.map,                          # list factors and constants for parameters
  responses,                      # Response (accumulator) names
  factors    = list(A = "1"),     # Factor names and levels
  match.map  = NULL,              # Scores responses
  constants  = numeric(0),        # Parameters set to constant value
  type       = "norm",            # model type
  posdrift   = TRUE,              # only used by norm
  verbose    = TRUE)              # Print p.vector, constants and type
{
  mapinfo <- check_BuildModel(p.map, responses, factors, match.map,
                                      constants, type)
  map.names  <- mapinfo[[1]]
  map.levels <- mapinfo[[2]] ## usually empty named list
  match.map  <- mapinfo[[3]] ## match.map fixed numeric sequences
  factors_long <- check_factors(responses, factors, match.map, map.names)
  dim2 <- make_pnames(p.map, factors_long) ## names.par
  npar <- length(dim2)

  ## Make level array for manifest design and accumulators
  level_array <- make_level_array(factors_long[1:(length(factors) + 1)])
  check_matchmap(match.map, level_array)
  dim1       <- as.vector(level_array)
  ncondition <- length(dim1)
  nresponse  <- length(responses)

  n1order   <- make_n1order(responses, level_array) ## norm type node 1
  matchcell <- make_match_cell(match.map, level_array) ## matching cells Boolean
  rdinfo    <- flipz(responses, match.map, type, level_array)

  MRmap  <- check_keywords(p.map, match.map, type, map.names)
  is.M   <- MRmap[[1]]
  is.R   <- MRmap[[2]]
  is.map <- MRmap[[3]]

  if ( any(is.map) ) {
    p.map.name <- lapply(p.map, function(x) {
      unlist(strsplit(x, "[.]"))[
        unlist(strsplit(x, "[.]")) %in% map.names]
    })
    maps_dim <- c(ncondition / nresponse, nresponse, nresponse)
    map.shuffle <- matrix(aperm(array(1:ncondition, dim=maps_dim),
                                c(1,3,2)), ncol = nresponse)
  }

  out <- array(NA, dim = c(ncondition, npar, nresponse))
  dimnames(out) <- list(dim1, dim2, responses)

  ## col.par = column parameter type (1st name)
  if ( is.null(match.map) ) {
    col.par.levels <- responses
    # cat("match.map is NULL")
    # print(col.par.levels)
  } else {
    ## col.par.levels because map.levels is an empty named list
    col.par.levels <- c(responses, "true", "false", map.levels)
    # cat("match.map found")
    # print(col.par.levels)
  }

  col.par <- strsplit(dim2, "[.]")
  col.fac <- lapply(col.par, function(x){x[-1]})
  col.par <- sapply(col.par, function(x){x[1]})

  ## split into fac and resp
  col.fac  <- lapply(col.fac, function(x){
    if ( length(x) == 0 ) out <- c(NA, NA)
    if ( length(x) == 1 ) {
      if ( x %in% col.par.levels )
        out <- c(NA, x) else out <- c(x,NA)
    }
    if ( length(x)>1 )
      if ( x[length(x)] %in% col.par.levels )
        out <- c(paste(x[-length(x)], collapse="."),x[length(x)]) else
          out <- paste(x, collapse=".")
        out
  })

  col.resp <- sapply(col.fac, function(x){x[2]})
  col.fac  <- sapply(col.fac, function(x){x[1]})

  row.fac <- strsplit(dim1, "[.]")
  ## row.resp <- unlist(lapply(row.fac,function(x){x[length(x)]}))
  row.fac <- sapply( row.fac, function(x){ paste(x[-length(x)], collapse=".") } )

  # Fill out array
  for ( p in unique(col.par) )
  { # parameters
    is.col <- p==col.par
    ncols <- sum(is.col)
    if ( ncols==1 ) out[,is.col,] <- TRUE else
    { # there are parameter subtypes
      for ( i in 1:ncols )
      { # each parameter subtype
        # select rows based on factors
        tmp <- col.fac[is.col][i]
        is.fac.row <- rep(TRUE, ncondition)
        if ( !is.na(tmp) ) is.fac.row[!grepl_dot(tmp, row.fac)] <- FALSE
        # set not applicable rows to false
        out[!is.fac.row,is.col,][,i,] <- FALSE
        if ( is.M[p] )
        { # has a match factor
          for ( j in names(match.map$M) )
          { # response cell
            correct.response <- match.map$M[[j]]
            is.rcell <- is.fac.row & grepl_dot(j, row.fac)
            for ( k in responses )
            { # responses
              if ( k==correct.response )
              {
                if ( grepl("true", col.resp[is.col][i]) )
                  out[,is.col,][is.rcell,i,k] <- TRUE else
                    out[,is.col,][is.rcell,i,k] <- FALSE

              } else {
                if ( grepl("false",col.resp[is.col][i]) )
                  out[,is.col,][is.rcell,i,k] <- TRUE else
                    out[,is.col,][is.rcell,i,k] <- FALSE
              }
            }
          }
        } else if ( is.R[p] ) {
          for ( k in responses )
            out[is.fac.row,is.col,k][,i] <- k==col.resp[is.col][i]
        }  else if ( is.map[p] ) {
          out[is.fac.row,is.col,][,i,] <-
            match.map[[ p.map.name[[p]] ]] [map.shuffle[is.fac.row,]]==col.resp[is.col][i]
        } else out[is.fac.row,is.col,][,i,] <- TRUE
      }
    }
  }

  if (any(is.na(out))) stop("Some cells of the map were not assigned!")

  # add in constants
  all.par <- out[1,,1]
  all.par[1:length(all.par)] <- NA
  if ( length(constants) > 0 ) {
    if ( !all(names(constants) %in% names(all.par)) )
      stop("Name(s) in constants not in p.map")
    all.par[names(constants)] <- constants
  }

  attr(out, "all.par")    <- all.par
  attr(out, "p.vector")   <- all.par[is.na(all.par)]
  attr(out, "par.names")  <- unique(col.par)
  attr(out, "type")       <- type
  attr(out, "factors")    <- factors
  attr(out, "responses")  <- responses
  attr(out, "constants")  <- constants
  attr(out, "posdrift")   <- posdrift
  attr(out, "n1.order")   <- n1order
  attr(out, "match.cell") <- matchcell
  if ( !is.null(match.map) ) attr(out, "match.map") <- match.map
  attr(out, "is.r1") <- rdinfo[[1]]
  attr(out, "bound") <- rdinfo[[2]]

  if (verbose) {
    cat("\nParameter vector names are: ( see attr(,\"p.vector\") )\n")
    print(names(all.par[is.na(all.par)]))
    cat("\nConstants are (see attr(,\"constants\") ):\n")
    print(constants)
    mod <- paste("\nModel type =", type)
    if (type == "norm") mod <- paste(mod, "(posdrift =", posdrift,")")
    cat(paste(mod, "\n\n"))
  }

  class(out) <- "model"
  return(out)
}

##' @rdname BuildModel
##' @export
print.model <- function(x, p.vector = NULL, ...) {

  if (!is.array(unclass(x)))
  {
    message("model is not an array. ")
    stop("Do you attempt to print posterior samples or multiple models?")
  }


  if (is.null(p.vector)) {
    nr <- length(attr(x, "response"))
    for (i in 1:nr) {
      dim3 <- dimnames(x)[[3]]
      cat(dim3[i], "\n")
      print(x[,, i])
    }
    message("Attributes: ")
    print(names(attributes(x)))
    return(invisible(x))

  } else {
    dim1 <- dimnames(x)[[1]]

    out <- lapply(dim1, function(xx) {
      print(xx)
      print(TableParameters(p.vector, xx, x, TRUE))
    })
    return(invisible(out))
  }
}

##' @rdname BuildModel
##' @importFrom utils head
##' @export
print.dmi <- function(x, ...) {
  model <- attr(x, "model")
  print.model(model)
  print(head(as.data.frame(x)))
  return(invisible(x))
}


##' Bind data and models
##'
##' Binding a data set with a model object. The function also checks whether
##' they are compatible and adds attributes on a data model instance.
##'
##' @param x data as in data frame
##' @param model a model object
##' @return a data model instance
##' @export
BuildDMI <- function(x, model) {

  res <- check_BuildDMI(x, model)
  subject_models <- res$issm
  modeli <- res$model
  fnams <- names(attr(modeli, "factors"))

  if ( any(names(x) == "s") ) { # more than one subject

    s    <- levels(x$s)
    nsub <- length(s)
    dat  <- vector("list", nsub)
    names(dat) <- s
    if (subject_models) names(model) <- s

    # is.sim <- !is.null(attr(dat, "parameters"))

    for (i in 1:nsub)
    {
      k <- s[i]
      if (subject_models) modeli <- model[[k]] else modeli <- model

      ## Remove s column and extract ith subject
      dat[[i]] <- x[x$s == k, names(x) != "s"]

      # add model and index attribute to data
      cells <- apply(dat[[k]][, c(fnams, "R")], 1, paste, collapse = ".")
      cell.index <- vector("list", dim(modeli)[1])
      names(cell.index) <- row.names(modeli)

      for ( j in names(cell.index) ) cell.index[[j]] <- cells %in% j

      attr(dat[[k]], "cell.index") <- cell.index
      attr(dat[[k]], "cell.empty") <- sapply(cell.index, function(xx){sum(xx)}) == 0
      attr(dat[[k]], "model") <- modeli

      # if (is.sim) attr(data[[s]], "parameters") <- attr(dat, "parameters")[s,]

      class(dat) <- c("dmi", "model", "list")
    }

    return(dat)

  } else { # one subject
    ## add model and index attribute to data
    ## Presume only one dependent variable RT
    ## x[, c(fnams, "R")] extract all columns except the dependent variable (RT)
    cells <- apply(x[, c(fnams, "R")], 1, paste, collapse = ".")
    ## ncell == nobservation

    cell.index <- vector("list", dim(model)[1])
    names(cell.index) <- row.names(model)
    ## scan trial-by-trial (every observation)
    for ( j in names(cell.index) ) cell.index[[j]] <- cells %in% j

    attr(x, "cell.index") <- cell.index
    attr(x, "cell.empty") <- sapply(cell.index, function(xx){sum(xx)}) == 0
    attr(x, "model") <- model
    class(x) <- c("dmi", "model", "data.frame")
    return(x)
  }
}

make_level_array <- function(x = NA) {
  if (all(is.na(x))) stop("Found no factors!")
  out <- x[[1]]
  nf  <- length(x)

  if ( nf > 1 ) {
    for (i in 2:nf) out <- outer(out, x[[i]], "paste", sep = ".")
    dimnames(out) <- x
  } else {
    out <- array(out, dim = length(x[[1]]), dimnames = x)
  }
  return(out)
}

##' Table response and parameter
##'
##' \code{TableParameters} arranges the values in a parameter
##' vector and creates a response x parameter matrix. The matrix is used
##' by the likelihood function, assigning a trial to a cell for calculating
##' probability densities.
##'
##' @param p.vector a parameter vector
##' @param cell a string or an integer indicating a design cell, e.g.,
##' \code{s1.f1.r1} or 1. Note the integer cannot exceed the number of cell.
##' One can check this by entering \code{length(dimnames(model))}.
##' @param model a model object
##' @param n1order a Boolean switch, indicating using node 1 ordering. This is
##' only for LBA-like models and its n1PDF likelihood function.
##' @return each row corresponding to the model parameter for a response.
##' When \code{n1.order} is FALSE, TableParameters returns a martix without
##' rearranging into node 1 order.  For example, this is used in
##' the \code{simulate} function. By default \code{n1.order} is TRUE.
##' @export
##' @examples
##' m1 <- BuildModel(
##'   p.map     = list(a = "1", v = "F", z = "1", d = "1", sz = "1", sv = "F",
##'                    t0 = "1", st0 = "1"),
##'   match.map = list(M = list(s1 = "r1", s2 = "r2")),
##'   factors   = list(S = c("s1", "s2"), F = c("f1","f2")),
##'   constants = c(st0 = 0, d = 0),
##'   responses = c("r1","r2"),
##'   type      = "rd")
##'
##' m2 <- BuildModel(
##'   p.map = list(A = "1", B = "1", mean_v = "M", sd_v = "1",
##'     t0 = "1", st0 = "1"),
##'   constants = c(st0 = 0, sd_v = 1),
##'   match.map = list(M = list(s1 = 1, s2 = 2)),
##'   factors   = list(S = c("s1", "s2")),
##'   responses = c("r1", "r2"),
##'   type      = "norm")
##'
##' pvec1 <- c(a = 1.15, v.f1 = -0.10, v.f2 = 3, z = 0.74, sz = 1.23,
##'            sv.f1 = 0.11, sv.f2 = 0.21, t0 = 0.87)
##' pvec2 <- c(A = .75, B = .25, mean_v.true = 2.5, mean_v.false = 1.5,
##'            t0 = .2)
##'
##' print(m1, pvec1)
##' print(m2, pvec2)
##'
##' accMat1 <- TableParameters(pvec1, "s1.f1.r1", m1, FALSE)
##' accMat2 <- TableParameters(pvec2, "s1.r1",    m2, FALSE)
##'
##' ##    a    v   t0    z d   sz   sv st0
##' ## 1.15 -0.1 0.87 0.26 0 1.23 0.11   0
##' ## 1.15 -0.1 0.87 0.26 0 1.23 0.11   0
##'
##' ##    A b  t0 mean_v sd_v st0
##' ## 0.75 1 0.2    2.5    1   0
##' ## 0.75 1 0.2    1.5    1   0
TableParameters <- function(p.vector, cell, model, n1order)
{
  pnames   <- names(attr(model, "p.vector"))
  allpar   <- attr(model, "all.par")
  parnames <- attr(model, "par.names")
  type     <- attr(model, "type")
  n1idx    <- attr(model, "n1.order")
  resp     <- attr(model, "responses")
  cell     <- check_cell(cell, model)
  isr1     <- check_rd(type, model)

  dim0 <- dimnames(model)[[1]]
  dim1 <- dimnames(model)[[2]]
  dim2 <- dimnames(model)[[3]]

  parmeter_matrix <- p_df(p.vector, cell, type, pnames, parnames, dim0, dim1,
                          dim2, allpar, model,  isr1, n1idx, n1order)

  out <- as.data.frame(parmeter_matrix)

  if(type == "rd")
  {
    names(out) <- c("a","v","z","d","sz","sv","t0","st0")
    rownames(out) <- attr(model, "response")
  }

  if(type %in% c("norm", "norm_pda", "norm_pda_gpu"))
  {
    if (dim(out)[[2]] != 6)   ## Prospective memory?
    {
      names(out) <- c("A", "b", "t0", "mean_v", "sd_v", "st0", "nacc")
    } else {
      names(out) <- c("A", "b", "t0", "mean_v", "sd_v", "st0")
    }
  }

  return(out)
}


######### Model checks  -----------------------------------
check_BuildModel <- function(p_map, responses, factors, match_map, constants,
                          type) {
  ## Check requried inputs supplied
  if (is.null(p_map)) stop("Must supply p.map")
  if (is.null(responses)) stop("Must supply responses")
  ## if (is.null(factors)) stop("Must supply factors")

  ## Check factors
  keywords <- c("1", "s", "R")
  if ( length(unlist(factors)) != length(unique(unlist(factors))) )
    stop("All factors levels must be unqiue")
  if (any(names(factors) %in% keywords))
    stop("Do not use 1, s, or R as a factor name")

  ## Check no parameter names have a dot
  has_dot <- sapply(strsplit(names(p_map), "[.]"), length) > 1
  if (any(has_dot)) {
    stop(paste("Dots not allowed in p.map names, fix:",
               paste(names(p_map)[has_dot]), "\n"))
  }

  ## Check R last if used
  if (any(sapply(p_map, function(x){any(x=="R") && x[length(x)]!="R"})))
    stop("R factors must always be last")

  ## Check responnses
  if (type == "rd") {
    if (is.null(match_map)) stop("Must specify supply a match.map for the DDM")
    if (length(responses) != 2) stop("DDM only applicable for two responses")
  }

  ## Check match.map (if supplied)
  if (!is.null(match_map)) {
    # Check structure
    if ( length(match_map) < 1 || class(match_map[[1]]) != "list" )
      stop("match.map must be a list of lists")
    # Check match.map contains at least name M
    if ( !any(names(match_map) %in% "M") )
      stop("match.map must have a list named M")
    map_names <- names(match_map)[names(match_map) != "M"]
    map_levels <- sapply(match_map[names(match_map) != "M"], levels)

    # convert match.map$M to responses and check
    if ( is.numeric(unlist(match_map$M)) )
      match_map$M <- lapply(match_map$M, function(x){responses[x]})

    if ( !all(unlist(match_map$M) %in% responses) )
      stop("match.map$M has index or name not in response names")
    if ( !(all(sort(responses)==sort(unique(unlist(match_map$M))))) )
      stop("Not all response names are scored by match.map$M")
    if ( length(match_map) > 1 &&
         !all(lapply(match_map[names(match_map)!="M"], class) == "factor") )
      stop("Entries in match.map besides M must be factors")
    if ( length(unlist(map_levels)) != length(unique(unlist(map_levels))) )
      stop("All match.map levels must be unqiue")
    # Check factors
    if ( any(names(factors) == "M") )
      stop("Do not use M as a factor name")
    if ( any(names(factors) %in% names(match_map)) )
      stop(paste(match_map, "used in match.map, can not use as a factor name"))
    if ( any(unlist(factors) %in% c("true","false")) )
      stop("\"true\" and \"false\" cannot be used as factor levels")
    if ( any(map_levels %in% c("true","false")) )
      stop("\"true\" and \"false\" cannot be used as match.map levels")
    if ( length(unlist(c(factors, map_levels))) !=
         length(unique(unlist(c(factors, map_levels)))) )
      stop("Factor levels cannot overlap match.map levels")

    # Check M and R are last
    if (any(sapply(p_map, function(x){any(x=="M") && x[length(x)]!="M"})))
      stop("M factors must always be last")

    ## Must return match_map. This is because the user is allowed to enter
    ## numberic sequences (1, 2, ...) to represent string response types (r1, r2).
    ## Internally, string (response) types are favoured. (should change this)
    out <- list(map_names, map_levels, match_map)

  } else {
    hasM <- sapply(p_map, function(x){ any(x == "M")})
    if(any(hasM)) stop("match.map is NULL, but found M factor in p.map")
    out <- NULL
  }

  return(out)
}

check_factors <- function(responses, factors,  match_map, map_names) {
  factors_long <- factors
  factors_long$R <- responses
  if (!is.null(match_map)) factors_long$M <- c("true", "false")

  # protect againt grep problems
  for (i in unlist(factors_long)) if ( length(grep(i,unlist(factors_long)))!=1 )
    stop("Factor, response or map level is not unique or is substring of another
      level or of \"true\" or \"false\"!" )

  # Add in extra match.map factors (if any)
  if ( !is.null(match_map) ) for (i in map_names)
    factors_long[[i]] <- levels(match_map[[i]])

  return(factors_long)
}

make_pnames <- function(p_map, factors_long) {
  ## Make parameter names
  names.par <- character(0)
  for( i in names(p_map) ) {

    if ( length(p_map[[i]])==1 && p_map[[i]] == "1" ) {
      ## If parameters are not affected by any factors, just use parameter names
      new.names <- i

    } else {
      ## Stick the factor levels onto parameter names, separated by a dot
      new.names <- paste(i, factors_long[[ p_map[[i]][1] ]], sep = ".")

      if ( length(p_map[[i]]) > 1 ) {
        ## more than 1 factor affect parameters
        for (j in 2:length(p_map[[i]])) {
            new.names <- as.vector(outer(
              new.names, factors_long[[p_map[[i]][j]]], "paste", sep = "."
            ))
        }
      }
    }
    ## Return parameter names + factor levels
    names.par <- c(names.par, new.names)
  }

  return(names.par)
}

check_matchmap <- function(match_map, level_array) {
  # # Check match.map
  if ( !is.null(match_map) ) for ( i in names(match_map$M) ) {
    match.position <- grep(i, level_array)
    if ( length(match.position) == 0 )
      stop(paste(i, "in match.map is not in the design"))
  }
  invisible(NULL)
}

check_keywords <- function(p.map, match.map, type, map.names) {
  ## Does the par use the match factor?
  is.M <- sapply(p.map, function(x){
    any(unlist(strsplit(x, "[.]") == "M"))
  })

  ## Does the par use a response factor
  is.R <- sapply(p.map, function(x){
    any(unlist(strsplit(x, "[.]") == "R"))
  })

  if (type == "rd" & (any(is.M) | any(is.R))) stop("Cannot use M or R in DDM p.map")

  # ## Does the par use a map factor (i.e., in match map but not M)
  if (!is.null(match.map)) {
    is.map <- sapply(p.map, function(x){
      any(unlist(strsplit(x, "[.]") %in% map.names))
    })
  } else {
    is.map <- logical(length(p.map))
    names(is.map) <- names(p.map)
  }

  if ( any(apply(cbind(is.M,is.R,is.map),1,function(x){sum(x)>1})) )
    stop("Parameters cannot have more than one of match.map and R factors")

  return(list(is.M, is.R, is.map))
}

make_match_cell <- function(match.map, level_array) {
  dim1 <- as.vector(level_array)
  ncondition <- length(dim1)
  out <- logical(ncondition)
  names(out) <- dim1

  if (!is.null(match.map)) {

    for (i in 1:length(match.map$M)) {
      match.num <- grep(match.map$M[i], dim1)
      idx1 <- match.num %in% grep(names(match.map$M[i]), dim1)
      idx2 <- match.num[idx1]
      out[idx2] <- TRUE
    }

  } else {
    message("match.map is NULL")
  }

  return(out)
}

make_n1order <- function(responses, level_array) {
  nresponse  <- length(responses)
  dim1       <- as.vector(level_array)
  resp <- sapply( strsplit(dim1, "[.]"), function(x){ x[length(x)] } )
  out  <- matrix(nrow = length(resp), ncol = nresponse)
  row.names(out) <- dim1

  for (i in 1:length(resp)) {
    out[i,] <- c(
      c(1:nresponse)[resp[i] == responses],
      c(1:nresponse)[resp[i] != responses]
    )
  }
  return(out)
}

flipz <- function(responses, match.map, type, level_array) {
  dim1       <- as.vector(level_array)
  ncondition <- length(dim1)

  if (type == "rd") # r1 cells have z flipped
  {
    is.r1 <- rep(FALSE, ncondition)

    idx1 <- sapply(lapply(
      as.list(names(match.map$M)[match.map$M == responses[1]]),
      function(x)grepl(x, dim1)
      ), function(x) which(x == TRUE)
    )

    is.r1[idx1] <- TRUE

    # add bound attributes
    bound <- rep("lower", ncondition)
    idx2 <- as.vector(sapply(
      paste("", names(match.map$M), match.map$M, sep = "*"),
      function(x){grep(glob2rx(x), dim1)}))

    bound[idx2] <- "upper"

    names(is.r1) <- dim1
    names(bound) <- dim1

  } else {
    is.r1 <- 0
    bound <- NULL
  }

  return(list(is.r1, bound))
}

check_BuildDMI <- function(data, model) {
  message1 <- "Model list is to match multiple subjects - models. No s column was found in data frame"
  message2 <- "Mostl list must be same length as the number of subjects"
  message3 <- "data must be a data frame"
  type <- attr(model, "type")

  if (is.list(model)) {
    if (!any(names(data) == "s")) stop(message1)
    if (length(model) != length(levels(data$s))) stop(message2)
    subject_models <- TRUE
    modeli <- model[[1]]
  } else {
    subject_models <- FALSE
    modeli <- model
  }

  fnams   <- names(attr(modeli, "factors"))
  factors <- attr(modeli, "factors")
  resps   <- attr(modeli, "responses")
  message4 <- paste("data must have columns named:",
                    paste(fnams, collapse = " "), "R", "RT")

  ## check data
  if ( !is.data.frame(data) ) stop(message3)
  if (!type %in% c("glm", "logit")) {
    if (!all(c(fnams, "R", "RT") %in% names(data)) ) stop(message4)
    if ( !all(sort(resps) == sort(levels(data[, "R"]))) )
      stop(paste("R must have levels:", paste(resps, collapse=" ")))
    if ( !is.numeric(data$RT) ) stop("RT must be of type numeric")
  }

  # i <- 1
  # data <- d
  for ( i in fnams ) {
    factori <- factors[[i]]
    if ( !all(sort(factori) == sort(levels(data[,i]))) )
      stop(paste("Factor", i, "must have levels:", paste(factori, collapse=" ")))
  }

  list(issm = subject_models, model=modeli)
}

checkddm1 <- function(p.map, responses,  type) {
  nr       <- length(responses)
  message1 <- "DDM only applicable for two responses"
  message2 <- "Cannot use M or R in DDM p.map"
  if (type == "rd" & (nr != 2)) stop(message1)

  # Does the par use a match factor or a response factor?
  is.M <- unlist(lapply(p.map, function(x){
    any(unlist(strsplit(x, "[.]") == "M"))
  }))
  is.R <- unlist(lapply(p.map, function(x){
    any(unlist(strsplit(x, "[.]") == "R"))
  }))
  if (type =="rd"  & ( any(is.M) | any(is.R) )) stop(message2)
}



checkdesign  <- function(match.map, levelarray) {
  ## Check match.map and expand
  for (i in names(match.map$M)) {
    match.position <- grep(i, levelarray)
    if ( length(match.position)==0 )
      stop(paste(i, "in match.map is not in the design"))
  }
}


##' Does a model object specify a correct p.vector
##'
##' Check a parameter vector
##'
##' @param p.vector parameter vector
##' @param model a model object
##' @export
check_pvec <- function(p.vector, model)
{
  modpvec <- names(attr(model, "p.vector"))
  ism1 <- modpvec %in% names(p.vector)
  ism2 <- names(p.vector) %in% modpvec
  bad  <- any(!ism1)
  if (bad) warning(paste("Parameter", modpvec[!ism1],
    "in model not present in p.vector\n"))
  bad <- bad | any(!ism2)
  if (any(!ism2)) warning(paste("Parameter",
    names(p.vector)[!ism2], "in p.vector not present in model\n"))
  invisible(bad)
}

check_cell <- function(cell, model) {

  dim1 <- dimnames(model)[[1]]
  if(is.numeric(cell)) {
    if ( cell > length(dim1) ) stop("cell out-of-range!")
    cell <- dim1[cell]
  }
  cell
}

check_rd <- function(type, model) {
  ## depreciate this
  if(type != "rd") { isr1 <- 0 } else { isr1 <- attr(model, "is.r1") }
  isr1
}




######### Simulation checks -----------------------------------------------------
createfacsdf <- function(model) {
  dim1 <- dimnames(model)[[1]]
  responses <- attr(model, "responses")
  levs      <- attr(model, "factors")
  nr <-  length(responses)

  facs  <- lapply(strsplit(dim1, "[.]"), function(x){x[-length(x)]})
  nf    <- length(facs)
  facs  <- facs[1:( nf/nr )]
  fnams <- names(levs)
  facs  <- data.frame(t(matrix(unlist(facs), length(fnams))))
  names(facs) <- fnams
  return(facs)
}

##' @importFrom data.table is.data.table
check_n <- function(n, facs) {
  if ( length(n) == 1 ) n <- rep(n, dim(facs)[1])
  test_n <- ifelse(is.data.frame(n), (nrow(n) != dim(facs)[1]),
    (length(n) != dim(facs)[1]))
  if (test_n) stop(paste("n must either be length 1 or", dim(facs)[1], "for this model."))

  if (data.table::is.data.table(n)) {
    n <- n$N
  } else {
    if ( !is.null(dimnames(n)) ) n <- merge(facs, data.frame(n), sort = F)$Freq
    for (i in 1:dim(facs)[1]) {if (n[i] < 0) warning("n is less than 0")}
  }

  return(n)
}


nadf <- function(n, facs, levs, type) {
  # create a data-frame container

  out <- data.frame(lapply(facs, rep, n))
  if (type == "logit") {
    for (i in names(levs)) {
      out[[i]] <- as.numeric(out[[i]]) - 1
    }
  } else {
    for (i in names(levs)) {
      out[[i]] <- factor(as.character(out[[i]]), levels = levs[[i]])
    }
  }


  out$R <- NA
  if (type == "glm") {
    out$X <- NA
    out$Y <- NA
  } else if (type == "logit") {
    out$Y <- NA
    out$N <- NA
  } else {
    out$RT <- NA
  }
  return(out)
}

FlipResponse_rd <- function(model, data, facs) {
  cell.names <- apply(data[, 1:length(names(facs)), drop = F], 1, paste, collapse=".")
  M <- attr(model, "match.map")$M
  R <- attr(model, "responses")
  for ( i in names(M) )
    if ( M[[i]] == R[1] )
      data[grep(i, cell.names),"R"] <- as.character(
        factor(as.character(data[grep(i, cell.names), "R"]),
          levels=R, labels=R[2:1]))
  return(data)
}


######### CUSTOM MAP MAKING for PM -------------------------------
MakeEmptyMap <- function(FR, levels) {
  ## This derives from DMC's empty.map
  level_array <- make_level_array(FR)
  map <- rep("", length = length(level_array))
  names(map) <- level_array
  factor(map, levels = levels)
}

AssignMap <- function(map, value = "", eq.list = list(), funs = NULL,
  include = NA, match_values = NA) {

  if (any(is.na(include))) ok <- rep(TRUE,length(map)) else
  {
    ok <- grepl(include[1],names(map))
    if (length(include)>1) for (i in 2:length(include))
      ok <- ok | grepl(include[i],names(map))
  }
  if ( length(eq.list)==0 ) # Fill in empty elements if included
    map[is.na(map) & ok] <- value else if (all(unlist(lapply(eq.list, length))==1))
  {
    if (all(is.na(match_values)) || length(match_values)!=length(eq.list))
      stop("If eq.list only has length one entries match_value must contain target for each entry")
    match.mat <- matrix(unlist(lapply(eq.list,function(x){
      (unlist(lapply(strsplit(names(map),".", fixed = T),function(y){
        y[x]})))})),ncol=length(eq.list))
    map[apply(match.mat,1,function(x){all(x == match_values)})] <- value
  } else {
    if ( is.null(funs) ) funs <- lapply(eq.list,function(x){
      list("identity","identity")
    })
    if (length(funs)!=length(eq.list))
      stop("Must specify one function pair in funs for each entry in eq.list")
    if ( class(funs)!="list" || !all(unlist(lapply(funs,class))=="list") )
      stop("funs must be  list of lists")
    if ( !all(unlist(lapply(funs,length))==2) )
      stop("Each entry in funs must be length 2")
    funs <- lapply(funs,function(x){lapply(x,function(y){
      if ( is.character(y) && y=="" ) "identity" else y
    })})
    pair.list <- lapply(eq.list,function(x){
      matrix(unlist(lapply(strsplit(names(map),".",fixed=T),function(y){
        y[x]})),nrow=2)})
    map.mat <- matrix(nrow=length(map),ncol=length(eq.list))
    for ( i in 1:length(eq.list) )
      map.mat[,i] <- apply(pair.list[[i]],2,function(x){
        do.call(funs[[i]][[1]],list(x[1])) ==
        do.call(funs[[i]][[2]],list(x[2]))
      })
      map[apply(map.mat,1,all) & ok] <- value
  }
  map
}

