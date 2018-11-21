######### CORE ------------------------------------------------------
##' Bayeisan Computation for Cognitive Models
##'
##' \pkg{ggdmc} evolves from Dynamic Models of Choice (DMC), using graphic
##' styles of ggplot2, highly efficient computations of Armadillo C++ and
##' Rcpp and connecting GPU parallel computation to \pkg{ppda} to make fitting
##' complex cognitive models feasible. \pkg{ggdmc} uses the sampling technique
##' of population-based Monte Chain Monte Carlo.
##'
##' @keywords package
##' @name ggdmc
##' @docType package
##' @author  Yi-Shin Lin <yishin.lin@utas.edu.au> \cr
##' Andrew Heathcote <andrew.heathcote@utas.edu.au>
##' @references
##' Heathcote, A., Lin, Y.-S., Reynolds, A., Strickland, L., Gretton, M. &
##' Matzke, D., (2018). Dynamic model of choice.
##' \emph{Behavior Research Methods}.
##' https://doi.org/10.3758/s13428-018-1067-y. \cr
##'
##' Turner, B. M., & Sederberg P. B. (2012). Approximate
##' Bayesian computation with differential evolution,
##' \emph{Journal of Mathematical Psychology}, 56, 375--385. \cr
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
##' Create a model array and attach many model attributes. These attributes
##' specify a particular model and parameterisation.
##'
##' @param p.map mapping factorial design to model parameters
##' @param match.map matching stimuli and responses
##' @param factors specifying factors and factor levels
##' @param constants setting parameters as constant value
##' @param responses Response (accumulator) names
##' @param type using character string to specifying model type.
##' @param posdrift enforce postive drift rate, using truncated normal or
##' just using normal distribution. This is used only by norm type (any
##' LBA variants and extensions)
##' @param verbose Print parameter vector, constants and model type
##' @param x a model object
##' @param p.vector parameter vector (for printing model)
##' @param ... other arguments
##' @importFrom utils glob2rx
##' @export
##' @examples
##' model <- BuildModel(
##'         p.map     = list(a = "1", v = "1", z = "1", d = "1", t0 = "1", sv = "1",
##'                          sz = "1", st0 = "1"),
##'         constants = c(st0 = 0, d = 0, sz = 0, sv = 0),
##'         match.map = list(M = list(s1 = "r1", s2 = "r2")),
##'         factors   = list(S = c("s1", "s2")),
##'         responses = c("r1", "r2"),
##'         type      = "rd")
BuildModel <- function(
  p.map,                          # list factors and constants for parameters
  responses,                      # Response (accumulator) names
  factors    = list(dummy = "1"), # Factor names and levels
  match.map  = NULL,              # Scores responses
  constants  = numeric(0),        # Parameters set to constant value
  type       = "norm",            # model type
  posdrift   = TRUE,              # only used by norm
  verbose    = TRUE               # Print p.vector, constants and type
  ) {

  # Check requried inputs supplied
  if (is.null(p.map)) stop("Must supply p.map")
  if (is.null(factors)) stop("Must supply factors")
  if (is.null(responses)) stop("Must supply responses")

  # Check factors
  if ( length(unlist(factors)) != length(unique(unlist(factors))) )
    stop("All factors levels must be unqiue")
  if ( any(names(factors) %in% c("1","R", "s")) )
    stop("Do not use 1, s, or R as a factor name")
  # Check no parameter names have a dot
  has.dot <- unlist(lapply(strsplit(names(p.map),".",fixed=TRUE),length))>1
  if ( any(has.dot) )
    stop(paste("Dots not allowed in p.map names, fix:",paste(names(p.map)[has.dot])))
  # Check R last if used
  if (any(unlist(lapply(p.map,function(x){any(x=="R") && x[length(x)]!="R"}))))
    stop("R factors must always be last")

  # Check responnses
  if ( type =="rd" ) {
    if (is.null(match.map))
      stop("Must specify supply a match.map for the DDM")
    if ( length(responses)!=2 )
      stop("DDM only applicable for two responses")
  }

  # Check match.map (if supplied)
  if (!is.null(match.map)) {
    # Check structure
    if ( length(match.map)<1 || class(match.map[[1]]) != "list" )
      stop("match.map must be a list of lists")
    # Check match.map contains at least name M
    if ( !any(names(match.map) %in% "M") )
      stop("match.map must have a list named M")
    map.names <- names(match.map)[names(match.map)!="M"]
    map.levels <- unlist(lapply(match.map[names(match.map)!="M"],levels))
    # convert match.map$M to responses and check
    if ( is.numeric(unlist(match.map$M)) )
      match.map$M <- lapply(match.map$M,function(x){responses[x]})
    if ( !all(unlist(match.map$M) %in% responses) )
      stop("match.map$M has index or name not in response names")
    if ( !(all(sort(responses)==sort(unique(unlist(match.map$M))))) )
      stop("Not all response names are scored by match.map$M")
    if ( length(match.map)>1 &&
        !all(lapply(match.map[names(match.map)!="M"],class)=="factor") )
      stop("Entries in match.map besides M must be factors")
    if ( length(unlist(map.levels)) != length(unique(unlist(map.levels))) )
      stop("All match.map levels must be unqiue")
    # Check factors
    if ( any(names(factors) == "M") )
      stop("Do not use M as a factor name")
    if ( any(names(factors) %in% names(match.map)) )
      stop(paste(map.names,"used in match.map, can not use as a factor name"))
    if ( any(unlist(factors) %in% c("true","false")) )
      stop("\"true\" and \"false\" cannot be used as factor levels")
    if ( any(map.levels %in% c("true","false")) )
      stop("\"true\" and \"false\" cannot be used as match.map levels")
    if ( length(unlist(c(factors,map.levels))) !=
        length(unique(unlist(c(factors,map.levels)))) )
      stop("Factor levels cannot overlap match.map levels")
    # Check M and R are last
    if (any(unlist(lapply(p.map,function(x){any(x=="M") && x[length(x)]!="M"}))))
      stop("M factors must always be last")

  }

  factors.short <- factors
  factors$R <- responses
  if (!is.null(match.map)) factors$M <- c("true","false")

  # protect againt grep problems
  for (i in unlist(factors)) if ( length(grep(i,unlist(factors)))!=1 )
    stop("Factor, response or map level is not unique or is substring of another
      level or of \"true\" or \"false\"!" )

  # Add in extra match.map factors (if any)
  if ( !is.null(match.map) ) for (i in map.names)
    factors[[i]] <- levels(match.map[[i]])


  # Make parameter names
  names.par <- character(0)
  for ( i in names(p.map) )
  {
    if ( length(p.map[[i]])==1 && p.map[[i]] == "1" ) new.names <- i else
    {
      new.names <- paste(i,factors[[p.map[[i]][1]]],sep=".")
      if ( length(p.map[[i]])>1 ) for ( j in 2:length(p.map[[i]]) )
        new.names <- as.vector(outer(
          new.names,factors[[p.map[[i]][j]]],"paste",sep="."
        ))
    }
    names.par <- c(names.par,new.names)
  }

  # Make level array for manifest design and accumulators
  level.array <- MakeLevelArray(factors[1:(length(factors.short)+1)])

  # Check match.map
  if ( !is.null(match.map) ) for ( i in names(match.map$M) ) {
    match.position <- grep(i,level.array)
    if ( length(match.position)==0 )
      stop(paste(i,"in match.map is not in the design"))
  }

  # Does the par use the match factor?
  is.M <- unlist(lapply(p.map,function(x){
    any(unlist(strsplit(x,".",fixed=T)=="M"))
  }))

  # Does the par use a response factor
  is.R <- unlist(lapply(p.map,function(x){
    any(unlist(strsplit(x,".",fixed=T)=="R"))
  }))

  if ( type =="rd"  & ( any(is.M) | any(is.R) ) )
    stop("Cannot use M or R in DDM p.map")

  # Does the par use a map factor (i.e., in match map but not M)
  if ( !is.null(match.map) ) {
    is.map <- unlist(lapply(p.map,function(x){
      any(unlist(strsplit(x,".",fixed=T) %in% map.names))
    }))
  } else {
    is.map <- logical(length(p.map))
    names(is.map) <- names(p.map)
  }

  if ( any(is.map) ) {
    p.map.name <- lapply(p.map,function(x){
      unlist(strsplit(x,".",fixed=T))[
        unlist(strsplit(x,".",fixed=T)) %in% map.names]
    })
    nr <- length(responses)
    n <- length(level.array)
    map.shuffle <- matrix(aperm(array(1:n,dim=c(n/nr,nr,nr)),c(1,3,2)),ncol=nr)
  }

  if ( any(apply(cbind(is.M,is.R,is.map),1,function(x){sum(x)>1})) )
    stop("Parameters cannot have more than one of match.map and R factors")

  # use.par = boolean matrix for parameter use, cells x pars x resposnes
  use.par <- array(NA,
    dim=c(length(level.array),length(names.par),length(responses)))
  dimnames(use.par) <-
    list(as.vector(level.array),names.par,responses)

  # col.par = column parameter type (1st name)
  if ( is.null(match.map) )
    col.par.levels <- responses else
      col.par.levels <- c(responses,"true","false",map.levels)

  col.par <- strsplit(dimnames(use.par)[[2]],".",fixed=T)
  col.fac <- lapply(col.par,function(x){x[-1]})
  col.par <- unlist(lapply(col.par,function(x){x[1]}))
  # split into fac and resp
  col.fac <- lapply(col.fac,function(x){
    if ( length(x)==0 ) out <- c(NA,NA)
    if ( length(x)==1 ) {
      if ( x %in% col.par.levels )
        out <- c(NA,x) else out <- c(x,NA)
    }
    if ( length(x)>1 )
      if ( x[length(x)] %in% col.par.levels )
        out <- c(paste(x[-length(x)],collapse="."),x[length(x)]) else
          out <- paste(x,collapse=".")
        out
  })
  col.resp <- unlist(lapply(col.fac,function(x){x[2]}))
  col.fac <- unlist(lapply(col.fac,function(x){x[1]}))

  row.fac <- strsplit(dimnames(use.par)[[1]],".",fixed=T)
  #  row.resp <- unlist(lapply(row.fac,function(x){x[length(x)]}))
  row.fac <- unlist(lapply(row.fac,function(x){
    paste(x[-length(x)],collapse=".")}))

  # Fill use.par array
  for ( p in unique(col.par) )
  { # parameters
    is.col <- p==col.par
    ncols <- sum(is.col)
    if ( ncols==1 ) use.par[,is.col,] <- TRUE else
    { # there are parameter subtypes
      for ( i in 1:ncols )
      { # each parameter subtype
        # select rows based on factors
        tmp <- col.fac[is.col][i]
        is.fac.row <- rep(TRUE,dim(use.par)[1])
        if ( !is.na(tmp) ) is.fac.row[!grepl_dot(tmp, row.fac)] <- FALSE
        # set not applicable rows to false
        use.par[!is.fac.row,is.col,][,i,] <- FALSE
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
                  use.par[,is.col,][is.rcell,i,k] <- TRUE else
                    use.par[,is.col,][is.rcell,i,k] <- FALSE

              } else {
                if ( grepl("false",col.resp[is.col][i]) )
                  use.par[,is.col,][is.rcell,i,k] <- TRUE else
                    use.par[,is.col,][is.rcell,i,k] <- FALSE
              }
            }
          }
        } else if ( is.R[p] ) {
          for ( k in responses )
            use.par[is.fac.row,is.col,k][,i] <- k==col.resp[is.col][i]
        }  else if ( is.map[p] ) {
          use.par[is.fac.row,is.col,][,i,] <-
            match.map[[ p.map.name[[p]] ]] [map.shuffle[is.fac.row,]]==col.resp[is.col][i]
        } else use.par[is.fac.row,is.col,][,i,] <- TRUE
      }
    }
  }

  if ( any(is.na(use.par)) )
    stop("Some cells of the map were not assigned!")

  # add in constants
  all.par <- use.par[1,,1]
  all.par[1:length(all.par)] <- NA
  if ( length(constants)>0 ) {
    if ( !all(names(constants) %in% names(all.par)) )
      stop("Name(s) in constants not in p.map")
    all.par[names(constants)] <- constants
  }

  attr(use.par, "all.par")  <- all.par
  attr(use.par, "p.vector") <- all.par[is.na(all.par)]

  if (length(attr(use.par,"p.vector"))<2)
    stop("DMC cant handle models with less than two parameters")

  if (verbose) {
    cat("\nParameter vector names are: ( see attr(,\"p.vector\") )\n")
    print(names(all.par[is.na(all.par)]))
    cat("\nConstants are (see attr(,\"constants\") ):\n")
    print(constants)
    mod <- paste("\nModel type =", type)
    if (type == "norm") mod <- paste(mod, "(posdrift =", posdrift,")")
    cat(paste(mod,"\n\n"))
  }


  # Names of parameter types (cannot have a period)
  attr(use.par, "par.names") <- unique(col.par)
  attr(use.par, "type") <- type
  attr(use.par, "factors") <- factors.short
  attr(use.par, "responses") <- responses
  attr(use.par, "constants") <- constants
  attr(use.par, "posdrift") <- posdrift

  par.df <- matrix(nrow=0,ncol=length(p.map))
  dimnames(par.df) <- list(NULL,names(p.map))
  par.df <- data.frame(par.df)

  # save "n1" orders
  resp <- unlist(lapply(strsplit(level.array,".",fixed=TRUE),function(x){
    x[length(x)]}))
  nr <- length(responses)
  n1.order <- matrix(nrow=length(resp),ncol=nr)
  for (i in 1:length(resp))
    n1.order[i,] <- c(c(1:nr)[resp[i]==responses],c(1:nr)[resp[i]!=responses])
  row.names(n1.order) <- row.names(use.par)

  # Boolean for matching cells
  match.cell <- logical(length(row.names(n1.order)))
  names(match.cell) <- row.names(n1.order)
  if ( !is.null(match.map) ) for (i in 1:length(match.map$M)) {
    match.num <- grep(match.map$M[i],names(match.cell))
    match.cell[match.num[match.num %in%
        grep(names(match.map$M[i]),names(match.cell))]] <- TRUE
  }

  attr(use.par, "n1.order") <- n1.order
  attr(use.par, "match.cell") <- match.cell
  if ( !is.null(match.map) ) attr(use.par, "match.map") <- match.map

  if (type=="rd") # r1 cells have z flipped
  {
    is.r1 <- rep(FALSE,length(row.names(use.par)))
    names(is.r1) <- row.names(use.par)
    is.r1[unlist(lapply(
      lapply(
        as.list(names(match.map$M)[match.map$M==responses[1]]),
        function(x)grepl(x,row.names(use.par))
      ),
      function(x) which(x==TRUE)))] <- TRUE
    attr(use.par,"is.r1") <- is.r1

    # add bound attributes
    bound <- rep("lower",dim(use.par)[1])
    bound[as.vector(sapply(
      paste("",names(attr(use.par,"match.map")$M),
        attr(use.par,"match.map")$M,sep="*"),
      function(x){grep(glob2rx(x),row.names(use.par))}))] <- "upper"
    names(bound) <- row.names(use.par)
    attr(use.par,"bound") <- bound
  }

  class(use.par) <- "model"
  return(use.par)
}

##' @rdname BuildModel
##' @export
print.model <- function(x, p.vector = NULL, ...) {

  if (!is.array(unclass(x))) stop("model is not an array. Do you attempt to print posterior samples?")

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
      print(TableParameters(p.vector, xx, x))
    })
    return(invisible(out))
  }
}

##' Bind data and models
##'
##' Binding a data with a model object. This function also checks
##' if they are compatible and adds a \code{cell.index} and many other
##' attributes to the data frame.
##'
##' @param x data in a data frame format
##' @param model a model object
##' @param npda number of model simulations
##' @param bw kernel bandwidth
##' @param gpuid GPU ID, indicating using which GPU cards on a machine with
##' multieple GPUs.
##' @param debug debugging switch
##' @export
BuildDMI <- function(x, model, npda = 16384, bw = 0.01, gpuid = 0,
  debug = FALSE) {

  res <- check_dmi(x, model)
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
      attr(dat[[k]], "cell.empty") <-
        unlist(lapply(cell.index, function(xx){sum(xx)})) == 0
      attr(dat[[k]], "model") <- modeli

      # if (is.sim) attr(data[[s]], "parameters") <- attr(dat, "parameters")[s,]

      attr(dat[[k]], "n.pda") <- npda
      attr(dat[[k]], "bw")    <- bw
      attr(dat[[k]], "gpuid") <- gpuid
      attr(dat[[k]], "debug") <- debug
      class(dat) <- c("model", "list")
    }

    return(dat)

  } else { # one subject
    # add model and index attribute to data
    cells      <- apply(x[, c(fnams, "R")], 1, paste, collapse = ".")
    cell.index <- vector("list", dim(model)[1])
    names(cell.index) <- row.names(model)
    for ( j in names(cell.index) ) cell.index[[j]] <- cells %in% j
    attr(x, "cell.index") <- cell.index
    attr(x, "cell.empty") <-
      unlist(lapply(cell.index, function(xx){sum(xx)})) == 0
    attr(x, "model") <- model
    attr(x, "n.pda") <- npda
    attr(x, "bw")    <- bw
    attr(x, "gpuid") <- gpuid
    attr(x, "debug") <- debug
    class(x) <- c("model", "data.frame")
    return(x)
  }
}


##' Create an array of all factorial combinations of factor levels
##'
##' Take a list storing one or multiple string vectors and glues the strings in
##' these vectors with a dot symbol.
##'
##' @param x a list storing one or multiple factor vectors. Each vector
##' can have different numbers of level.  For example,
##' \code{f <- list(E = c("nc", "wc"), S = c("n", "w", "pn", "pw"))}
##' @return a table showing the combinations of factor levels.
##'
##' @examples
##' ## Example 1
##' factors <- list(S = c("s1", "s2"))
##' MakeLevelArray(factors)
##' ##   S
##' ##  s1   s2
##' ## "s1" "s2"
##'
##' factors <- list(S = c("left", "right"))
##' MakeLevelArray(factors)
##' ##      S
##' ##   left   right
##' ##  "left" "right"
##'
##' ## Example 2
##' factors <- list(A = c("00", "45", "90", "135", "180"),
##'                 S = c("mirror", "normal"))
##' MakeLevelArray(factors)
##' ##       S
##' ## ------------------------------40
##' ##   A        mirror       normal
##' ## ------------------------------40
##' ##  00   "00.mirror"  "00.normal"
##' ##  45   "45.mirror"  "45.normal"
##' ##  90   "90.mirror"  "90.normal"
##' ## 135  "135.mirror" "135.normal"
##' ## 180  "180.mirror" "180.normal"
##'
##' ## Example 3
##' factors <- list(E = c("nc", "wc"),
##'                 S = c("n", "w", "pn", "pw"))
##' MakeLevelArray(factors)
##'
##' ##         S
##' ## -------------------------------40
##' ## E       n      w      pn      pw
##' ## -------------------------------40
##' ## nc "nc.n" "nc.w" "nc.pn" "nc.pw"
##' ## wc "wc.n" "wc.w" "wc.pn" "wc.pw"
##'
##' @export
MakeLevelArray <- function(x = NA) {
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
##' vector to a factorial response x parameter matrix. The matrix is used
##' by likelihood functions, assigning a trial to a cell for calculating
##' probability densities.
##'
##' @param x a parameter vector
##' @param cell a string or an integer indicating a design cell, e.g.,
##' \code{s1.f1.r1} or 1. Note the integer cannot exceed the number of cell.
##' use \code{length(dimnames(model))} to check the upper bound.
##' @param model a model ibhect
##' @param n1order a Boolean switch, indicating using node 1 ordering. This is
##' only for LBA-like models and its n1PDF likelihood function.
##' @return each row corresponding to the model parameter for a response.
##' When \code{n1.order} is FALSE, TableParameters returns a martix in natural
##' order, which is used by \code{simulate}. By default \code{n1.order} is TRUE,
##' the returned matrix, used by n1PDF-like functions.
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
TableParameters <- function(x, cell, model, n1order = TRUE) {
  pnames   <- names(attr(model, "p.vector"))
  allpar   <- attr(model, "all.par")
  parnames <- attr(model, "par.names")
  type     <- attr(model, "type")
  n1       <- attr(model, "n1.order")
  resp     <- attr(model, "responses")
  cell     <- check_cell(cell, model)
  isr1     <- check_rd(type, model)

  out <- as.data.frame(p_df(x, cell, pnames, allpar, parnames,
    model, type, dimnames(model)[[1]], dimnames(model)[[2]],
    dimnames(model)[[3]], isr1, n1, n1order))

  if(type == "rd") {
    names(out) <- c("a","v","z","d","sz","sv","t0","st0")
    rownames(out) <- attr(model, "response")
  }

  if(type %in% c("norm", "norm_pda", "norm_pda_gpu")) {
    if (dim(out)[[2]] != 6) {
      ## Prospective memory?
      names(out) <- c("A", "b", "t0", "mean_v", "sd_v", "st0", "nacc")
    } else {
      names(out) <- c("A", "b", "t0", "mean_v", "sd_v", "st0")
    }
  }

  if(type %in% c("plba0_gpu")) {
    names(out) <- c("A", "b", "mean_v", "sd_v", "mean_w", "rD", "t0","swt")
  }

  if(type %in% c("plba1", "plba1_gpu")) {
    names(out) <- c("A", "b", "mean_v", "sd_v", "mean_w", "rD", "t0","swt")
  }

  if(type %in% c("plba2")) {
    names(out) <- c("A","b", "mean_v","mean_w","sd_v", "sd_w", "rD","t0","swt")
  }

  if(type %in% c("plba3")) {
    names(out) <- c("A", "B", "C","mean_v","mean_w","sd_v","sd_w", "rD", "tD",
      "t0","swt")
  }

  if(type %in% c("lnr")) {
    names(out) <- c("meanlog", "sdlog", "t0", "st0")
    rownames(out) <- attr(model, "response")
  }

  if(type %in% c("cnorm")) {
    names(out) <- c("A", "b", "t0", "mean_v", "sd_v", "corr_v", "st0")
    rownames(out) <- attr(model, "response")
  }

  return(out)
}

##' Scoring RT data
##'
##' A convenient function to calculate mean, interquantile range, standard
##' deviation for correct and error RTs
##'
##' @param x a parameter vector
##' @param digits printing digits
##' @importFrom stats IQR
##' @importFrom stats sd
##' @export
score <- function(x, digits = 2) {
  correct <- tolower(x$S) == tolower(x$R)
  cat("Accuracy: \n");
  print(round(mean(correct), digits))
  mrt <- tapply(x$RT, list(correct), mean)
  iqr <- tapply(x$RT, list(correct), IQR)
  SD  <- tapply(x$RT, list(correct), sd)
  out <- data.frame(rbind(mrt, iqr, SD))
  names(out) <- c("error", "correct")
  print(round(out, digits))
  invisible(out)
}



######### Model checks  -----------------------------------
checknull <- function(p.map, factors, responses) {
  if (is.null(p.map)) stop("Must supply p.map")
  if (is.null(factors)) stop("Must supply factors")
  if (is.null(responses)) stop("Must supply responses")
}

checkmap <- function(match.map, responses) {

  ## Extract essential information
  nmap      <- length(match.map)
  mapnames  <- names(match.map)
  message1 <- "match.map must be a list of lists"
  message2 <- "match.map has no list named M"
  message3 <- "match.map$M has index or name not in response names. If you use
  integers, it must be continuous."
  message4 <- "Not all response names are scored by match.map$M"
  message5 <- "Entries in match.map besides M must be factors"

  if (nmap < 1 || class(match.map[[1]]) != "list") stop(message1)
  if (!any(names(match.map) %in% "M") ) stop(message2)
  map.names  <- mapnames[mapnames != "M"]
  umap       <- unlist(match.map$M) ## unlisted match map

  ## We allow the user to enter integers (1, 2, ..., n) to represent
  ## "response 1", "response 2", ... "response n". If "is.numeric" evaluates to
  ## TRUE, it means the user enters integers as response types. Hence, we use
  ## lapply to pick sequentially 1st stimulus type matching to 1st response
  ## type etc. Response types are stored in another variable "responses"
  if (is.numeric(umap)) {
    match.map$M <- lapply(match.map$M, function(x){responses[x]})
    umap <- unlist(match.map$M)
  }
  if ( !all(umap %in% responses) ) stop(message3)
  if ( !(all(sort(responses) == sort(unique(umap)))) ) stop(message4)
  if ( nmap > 1 && !all(lapply(match.map[mapnames != "M"], class)=="factor") )
    stop(message5)
  if ( length(match.map)>1 &&
      !all(lapply(match.map[names(match.map)!="M"], class)=="factor") )
    stop("Entries in match.map besides M must be factors")

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

checkddm3 <- function(pmat, i, model, p.prior, debug = FALSE) {
  ## RT = .3, precision = 2.5 are dummies
  if (anyNA(p.prior)) stop("p.prior not found in checkddm3")
  ps <- c(pmat$a[1], pmat$v[1], pmat$z[1], pmat$d[1], pmat$sz[1], pmat$sv[1],
    pmat$t0[1], pmat$st0[1], .3, 2.5)
  isbadparameter <- checkddm2(ps)
  j <- 0

  repeat {
    if (isbadparameter) {
      j <- j+1
      ps <- rprior(p.prior)
      pmat <- TableParameters(ps, i, model, FALSE)
      ps <- c(pmat$a[1], pmat$v[1], pmat$z[1], pmat$d[1], pmat$sz[1], pmat$sv[1],
        pmat$t0[1], pmat$st0[1], .3, 2.5)
      isbadparameter <- checkddm2(ps)
      if (debug) message("ps changed")
    } else if (j >= 1e3) {
      stop("Fail to simulate DDM data");
      break
    } else {
      break
    }

  }
  return(pmat)
}

checkfacresp <- function(p.map, match.map, factors, responses) {


  mapnames   <- names(match.map)
  map.names  <- mapnames[mapnames != "M"]
  map.levels <- unlist(map.names, levels)

  ufac <- unlist(factors)    ## unlisted factors
  uml  <- unlist(map.levels) ## unlisted map.levels
  reservednames <- c("1", "R", "R2", "M", "s")
  tfnames       <- c("true", "false")
  ufacmaplv     <- unlist(c(factors, map.levels))

  message1 <- "Do not use s, M, R, R2 or 1 as a factor name"
  message2 <- paste(map.names,"used in match.map, can not use as a factor name")
  message3 <- "All factors levels must be unqiue"
  message4 <- "All match.map levels must be unqiue"
  message5 <- "\"true\" and \"false\" cannot be used as factor levels"
  message6 <- "\"true\" and \"false\" cannot be used as match.map levels"
  message7 <- "Factor levels cannot overlap match.map levels"

  # Check factors and add responses
  if ( any(names(factors) %in% reservednames) ) { stop(message1) }
  if ( any(names(factors) %in% mapnames) ) {      stop(message2) }
  if ( length(ufac) != length(unique(ufac)) ) {   stop(message3) }
  if ( length(uml) != length(unique(uml)) ) {     stop(message4) }
  if ( any(ufac %in% tfnames) ) {                 stop(message5) }
  if ( any(map.levels %in% tfnames) ) {           stop(message6) }
  if ( length(ufacmaplv) != length(unique(ufacmaplv)) ) { stop(message7) }
}

checkpmap    <- function(p.map) {
  pmapnames <- names(p.map)
  has.dot   <- unlist(lapply(strsplit(pmapnames, "[.]"), length)) > 1
  message1 <- paste("Dots not allowed in p.map names, fix:", pmapnames[has.dot])
  if (any(has.dot)) stop(message1)

  # Check M and R are last
  if (any(unlist(lapply(p.map, function(x){any(x == "M") && x[length(x)]!="M"}))))
    stop("M factors must always be last")
  if (any(unlist(lapply(p.map, function(x){any(x == "R") && x[length(x)]!="R"}))))
    stop("R factors must always be last")
}

addmap2facs  <- function(match.map, factors, responses) {
  umap <- unlist(match.map$M)
  if (is.numeric(umap)) match.map$M <- lapply(match.map$M, function(x){responses[x]})
  mapnames   <- names(match.map)
  map.names  <- mapnames[mapnames != "M"]

  factors$R     <- responses
  if (!is.null(match.map)) {
    factors$M <- c("true", "false")
    for (i in map.names) factors[[i]] <- levels(match.map[[i]])
  }
  return(factors)
}

getparnames  <- function(p.map, factors) {
  ## Make parameter names
  out <- character(0)
  for ( i in names(p.map) )
  {
    if ( length(p.map[[i]]) == 1 && p.map[[i]] == "1" ) {
      new.names <- i
    } else {
      new.names <- paste(i, factors[[p.map[[i]][1]]], sep=".")
      if ( length(p.map[[i]]) > 1 ) {
        for (j in 2:length(p.map[[i]])) {
          new.names <- as.vector(outer(
            new.names, factors[[p.map[[i]][j]]], "paste", sep="."))
        }
      }
    }
    out <- c(out, new.names)
  }

  return(out)
}

checkdesign  <- function(match.map, levelarray) {
  ## Check match.map and expand
  for (i in names(match.map$M)) {
    match.position <- grep(i, levelarray)
    if ( length(match.position)==0 )
      stop(paste(i, "in match.map is not in the design"))
  }
}

splitfacresp <- function(colfac, responses, maplevels) {
  res <- lapply(colfac, function(x){
    if ( length(x) == 0 ) out <- c(NA, NA)
    if ( length(x) == 1 ) {
      if ( x %in% c(responses, "true", "false", maplevels) )
        out <- c(NA, x) else out <- c(x, NA)
    }
    if ( length(x) > 1 )
      if ( x[length(x)] %in% c(responses, "true", "false", maplevels) )
        out <- c(paste(x[-length(x)], collapse="."), x[length(x)]) else
          out <- paste(x, collapse=".")
        out
  })
  return(res)
}

## Check reserved factor name (M, R, and other user indicated names in match.map)
checkrn <- function(p.map, match.map) {
  mapnames   <- names(match.map)
  map.names  <- mapnames[mapnames != "M"]
  ## Reserved uppercase letters we don't want the user to use
  is.M <- unlist(lapply(p.map, function(x) {
    any(unlist(strsplit(x, "[.]") == "M"))
  }))

  ## Does the par use a response factor
  is.R <- unlist(lapply(p.map,function(x){
    any(unlist(strsplit(x, "[.]") == "R"))
  }))

  ## Does the par use a map factor, other than M
  is.map <- unlist(lapply(p.map,function(x) {
    any(unlist(strsplit(x, "[.]") %in% map.names))
  }))

  out <- cbind(is.M, is.R, is.map)

  if ( any(apply(out, 1, function(x){ sum(x) > 1 })) ) {
    message("M and R are reserved names. Also if you have used a name in")
    message("match.map other than M. It should not be used to name a paramter.")
    message("I will stop here. Please check your match.map or")
    stop("consult the package maintainer.")
  }
  # if ( any(is.map) ) {
  #   p.map.name <- lapply(p.map, function(x){
  #     unlist(strsplit(x,"[.]" ))[
  #       unlist(strsplit(x,"[.]")) %in% map.names]
  #   })
  #   nr <- length(responses)
  #   n  <- length(level.array)
  #   map.shuffle <- matrix(aperm(array(1:n, dim = c(n/nr, nr, nr)), c(1,3,2)), ncol=nr)
  # }

  if (any(is.map)) { stop("map.shuffle is pending.") }
  return(out)
}

# MakeModelArray <- function(match.map, responses, col.par, rnamemat, allpar,
#   level.array) {
#
#   umap <- unlist(match.map$M)
#   if (is.numeric(umap)) match.map$M <- lapply(match.map$M, function(x){responses[x]})
#   mapnames   <- names(match.map)
#   map.names  <- mapnames[mapnames != "M"]
#   map.levels <- unlist(map.names, levels)
#
#
#   modeldn  <- list(as.vector(level.array), allpar, responses)
#   modeldim <- c(length(level.array), length(allpar), length(responses))
#   use.par  <- array(NA, modeldim, modeldn)
#   is.M   <- rnamemat[,1]
#   is.R   <- rnamemat[,2]
#   is.map <- rnamemat[,3]
#
#   ## col.par = column parameter type (1st name)
#   col.par <- strsplit(modeldn[[2]], "[.]")
#   col.fac <- lapply(col.par, function(x){x[-1]})
#   col.par <- unlist(lapply(col.par, function(x){x[1]}))
#   ## split into fac and resp
#   col.fac <- splitfacresp(col.fac, responses, map.levels)
#   col.resp <- unlist(lapply(col.fac, function(x){x[2]}))
#   col.fac  <- unlist(lapply(col.fac, function(x){x[1]}))
#
#   row.fac <- strsplit(modeldn[[1]], "[.]")
#   row.fac <- unlist(lapply(row.fac, function(x){
#     paste(x[-length(x)], collapse=".")}))
#
#
#   for ( p in unique(col.par) ) {
#     is.col <- p==col.par
#     ncols  <- sum(is.col)
#     if ( ncols == 1 ) {
#       use.par[, is.col,] <- TRUE
#     } else {
#       # there are parameter subtypes
#       for ( i in 1:ncols ) {
#         ## each parameter subtype
#         ## select rows based on factors
#         tmp <- col.fac[is.col][i]
#         is.fac.row <- rep(TRUE, dim(use.par)[1])
#         if ( !is.na(tmp) ) is.fac.row[!grepl(tmp, row.fac)] <- FALSE
#         ## set not applicable rows to false
#         use.par[!is.fac.row, is.col, ][,i,] <- FALSE
#         if ( is.M[p] ) {
#           # has a match factor
#           for ( j in names(match.map$M) ) {
#             # response cell
#             correct.response <- match.map$M[[j]]
#             is.rcell <- is.fac.row & grepl(j, row.fac)
#             for ( k in responses ) {
#               # responses
#               if ( k==correct.response ) {
#                 if ( grepl("true", col.resp[is.col][i]) )
#                   use.par[,is.col,][is.rcell,i,k] <- TRUE else
#                     use.par[,is.col,][is.rcell,i,k] <- FALSE
#               } else {
#                 if ( grepl("false", col.resp[is.col][i]) )
#                   use.par[,is.col,][is.rcell,i,k] <- TRUE else
#                     use.par[,is.col,][is.rcell,i,k] <- FALSE
#               }
#             }
#           }
#         } else if ( is.R[p] ) {
#           for ( k in responses )
#             use.par[is.fac.row,is.col,k][,i] <- k==col.resp[is.col][i]
#         }  else if ( is.map[p] ) {
#           use.par[is.fac.row,is.col,][,i,] <-
#             match.map[[ p.map.name[[p]] ]] [map.shuffle[is.fac.row,]]==col.resp[is.col][i]
#         } else use.par[is.fac.row,is.col,][,i,] <- TRUE
#       }
#     }
#   }
#
#   if ( any(is.na(use.par)) ) stop("Some cells of the map were not assigned!")
#   return(use.par)
# }

getcolpar <- function(names.par) {
  col.par <- strsplit(names.par, "[.]")   ## dim2
  col.fac <- lapply(col.par, function(x){x[-1]})
  col.par <- unlist(lapply(col.par, function(x){x[1]}))
  return(col.par)
}

checkconst <- function(use.par, constants) {
  all.par <- use.par[1,,1]  ## First array, first row
  all.par[1:length(all.par)] <- NA
  if (length(constants)>0) {
    if ( !any(names(constants) %in% names(all.par)) )
      stop("Name(s) in constants not in p.map")
    all.par[names(constants)] <- constants
  }
  return(all.par)
}

addisr1 <- function(match.map, responses, use.par) {
  is.r1 <- rep(FALSE, length(row.names(use.par)))
  names(is.r1) <- row.names(use.par)
  is.r1[unlist(lapply(
    lapply(
      as.list(names(match.map$M)[match.map$M==responses[1]]),
      function(x)grepl(x, row.names(use.par))
    ),
    function(x)which(x==TRUE)))] <- TRUE
  attr(use.par, "is.r1") <- is.r1
  return(use.par)
}

getn1order <- function(responses, level.array, use.par) {
  resp <- unlist(lapply(strsplit(level.array,"[.]"),function(x){
    x[length(x)]}))
  nr <- length(responses)
  n1.order <- matrix(nrow=length(resp), ncol=nr)
  for (i in 1:length(resp))
    n1.order[i,] <- c(c(1:nr)[resp[i]==responses], c(1:nr)[resp[i]!=responses])
  row.names(n1.order) <- row.names(use.par)
  return(n1.order)

}

getmatchcell <- function(match.map, n1.order) {
  # Boolean for matching cells
  match.cell <- logical(length(row.names(n1.order)))
  names(match.cell) <- row.names(n1.order)
  for (i in 1:length(match.map$M)) {
    match.num <- grep(match.map$M[i],names(match.cell))
    match.cell[match.num[match.num %in%
        grep(names(match.map$M[i]),names(match.cell))]] <- TRUE
  }
  return(match.cell)
}

##' Check if parameter vector compatible with model object?
##'
##' Check if the user supplies a parameter vector, compatiable with
##' the model object.
##'
##' @param x parameter vector
##' @param model a model object
##' @export
check_pvec <- function(x, model) {
  modpvec <- names(attr(model, "p.vector"))
  ism1 <- modpvec %in% names(x)
  ism2 <- names(x) %in% modpvec
  bad  <- any(!ism1)
  if (bad) warning(paste("Parameter", modpvec[!ism1],
    "in model not present in p.vector\n"))
  bad <- bad | any(!ism2)
  if (any(!ism2)) warning(paste("Parameter",
    names(x)[!ism2], "in p.vector not present in model\n"))
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


check_dmi <- function(data, model) {
  ## Remember to make CheckDMI and check_dmi consistent
  message1 <- "Model list is for multiple subjects\nNo s column is found in
 data model instance"
  message2 <- "List of models must be same length as number of subjects"
  message3 <- "data must be a data frame"

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
  message4 <- paste("data must have columns named:", paste(fnams, collapse=" "),
    "R", "RT")

  if ( !is.data.frame(data) ) stop(message3)

  # check data
  if ( !all(c(fnams, "R", "RT") %in% names(data)) ) stop(message4)
  for ( i in fnams ) {
    factori <- factors[[i]]
    if ( !all(sort(factori) == sort(levels(data[,i]))) )
      stop(paste("Factor", i, "must have levels:", paste(factori, collapse=" ")))
  }

  if ( !all(sort(resps) == sort(levels(data[, "R"]))) )
    stop(paste("data$R must have levels:", paste(resps, collapse=" ")))
  if ( !is.numeric(data$RT) ) stop("data$RT must be of type numeric")

  list(issm = subject_models, model=modeli)
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


nadf <- function(n, facs, levs) {
  # create a data-frame container
  data <- data.frame(lapply(facs, rep, n))

  for (i in names(levs)) data[[i]] <- factor(as.character(data[[i]]), levels=levs[[i]])
  data$R <- NA
  data$RT <- NA
  return(data)
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
  level_array <- MakeLevelArray(FR)
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

