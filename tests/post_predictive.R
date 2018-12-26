predict_one <- function(samples, npost, rand, factors, seed) {
  model <- attributes(samples$data)$model
  facs <- names(attr(model, "factors"))
  class(samples$data) <- c("data.frame", "list")

  if (!is.null(factors)) {
    if (any(is.na(factors))) factors <- facs
    if (!all(factors %in% facs)) 
      stop(paste("Factors argument must contain one or more of:", 
                 paste(facs, collapse=",")))
  }
  
  resp <- names(attr(model, "responses"))
  ns   <- table(samples$data[,facs], dnn = facs)
  npar   <- samples$n.pars
  nchain <- samples$n.chains
  nmc    <- samples$nmc
  ntsample <- nchain * nmc
  pnames   <- samples$p.names
  # str(samples$theta) ## nchain x npar x nmc
  # str(thetas)        ## (nchainx nmc) x npar
  thetas <- matrix(aperm(samples$theta, c(3,1,2)), ncol = npar)
  colnames(thetas) <- pnames
  
  if (is.na(n.post)) {
    use <- 1:ntsample 
  } else {
    if (rand) {
      use <- sample(1:ntsample, n.post, replace = F) 
    } else {
      use <- round(seq(1, ntsample, length.out = n.post))
    }
  }
  
  npost  <- length(use)
  posts   <- thetas[use, ]
  nttrial <- sum(ns) ## number of total trials
  out <- data.frame(matrix(nrow = npost*nttrial, ncol = dim(samples$data)[2]))

  # for (i in 1:npost) {
  #   tmp <- ggdmc:::simulate_one(model, n = ns, ps = posts[i,], seed = NULL)
  #   out[(1+(i-1)*nttrial):(i*nttrial), names(tmp)] <- tmp
  #   if ( (i %% report) == 0) cat(".")
  # }
  
  ## should replace with parallel
  v <- lapply(1:npost, function(i) {
    ggdmc:::simulate_one(model, n = ns, ps = posts[i,], seed = NULL)
  })
  out <- data.table::rbindlist(v)
  names(out) <- names(samples$data)
  return(out)
}


tmp0 <- predict_one(fit, 100, TRUE, NA, NULL)
head(tmp0)
tail(tmp0)
str(tmp0)  

predict_one <- function(samples, n.post=100, probs=c(1:99)/100, rand=TRUE,
                       bw="nrd0",report=10,save.simulation=FALSE,factors=NA,
                       save.simulation.as.attribute=FALSE,ignore.R2=FALSE,censor=c(NA,NA),
                       gglist=FALSE, probs.gglist=c(0.1, 0.5, 0.9),CI.gglist=c(0.025, 0.975))
{
  
  samples <- fit
  n.post <- 100
  rand <- TRUE
  factors <- NA
  res <- prepare_one(samples, n.post, rand, factors)
  sim   <- res[[1]]
  posts <- res[[2]]
  ns    <- res[[3]]
  nttrial <- sum(ns)
  report=10
  

  # probs=c(1:99)/100
  
  # bw="nrd0",
  # save.simulation=FALSE,
  # save.simulation.as.attribute=FALSE,ignore.R2=FALSE,censor=c(NA,NA),
  # gglist=FALSE, 
  # probs.gglist=c(0.1, 0.5, 0.9),
  # CI.gglist=c(0.025, 0.975)

  # make list of posterior preditive density, quantiles and response p(robability)
  # NB: quantiles only calcualted for 2 or more RTs
  
  sim <- prepare_df(samples, n.post, rand, factors)
  cat(paste0("\nSimulating (\'.\'=", report, "): "))

  str(sim)
  str(posts)
  dplyr::tbl_df(sim)
  dplyr::tbl_df(posts)
  
  posts[i, ]
  ns

  head(ggdmc:::simulate_one)
  dplyr::tbl_df(tmp)
  dplyr::tbl_df(sim)
  tail(sim)
  i <- 2
  
  for (i in 1:n.post) {
    tmp <- ggdmc:::simulate_one(model, n = ns, ps = posts[i,], seed = NULL)
    sim[(1+(i-1)*nttrial):(i*nttrial), names(tmp)] <- tmp
    if ( (i %% report) == 0) cat(".")
  }

  reps <- rep(1:n.post,each=dim(samples$data)[1])
  
  if (!is.na(censor[1])) fast <- sim[,"RT"] < censor[1] else fast <- rep(FALSE,dim(sim)[1])
  if (!is.na(censor[2])) slow <- sim[,"RT"] > censor[2] else slow <- rep(FALSE,dim(sim)[1])
  ok <- !fast & !slow 
  sim <- sim[ok,]
  reps <- reps[ok]
  
  if ( save.simulation ) {
    sim <- cbind(reps,sim)
    attr(sim,"data") <- samples$data
    sim
  } else {
    sim.dqp <- get.dqp(sim,facs=factors,probs,n.post,ns=ns,bw=bw)
    dat.dqp <- get.dqp(sim=samples$data,factors,probs,bw=bw)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    out <- c(sim.dqp,dat.dqp)
    dpqs <- vector(mode="list",length=length(n.post))
    for (i in 1:n.post) {
      simi <- sim[reps==i,]
      dpqs[[i]] <- get.dqp(simi,factors,probs,1)
    }
    attr(out,"dpqs") <- dpqs
    if (save.simulation.as.attribute) 
      attr(out,"sim") <- cbind(reps,sim)
    if (gglist) attr(out, "gglist") <- 
      get.fitgglist.dmc(sim=cbind(reps,sim),data=samples$data,factors=factors, noR=FALSE, 
                        quantiles.to.get= probs.gglist, CI = CI.gglist)
    out
  }
}



predict_hier <- function(object, npost = 100, rand = TRUE, factors = NA, 
                         xlim = NA, seed = NULL) {
  hyper <- attr(object, "hyper")
  nhpar <- hyper$n.pars
  pp.prior <- hyper$pp.prior
  
  phi1 <- matrix(aperm(hyper$phi[[1]], c(3, 1, 2)), ncol = nhpar)
  phi2 <- matrix(aperm(hyper$phi[[2]], c(3, 1, 2)), ncol = nhpar)
  str(phi1)
  str(phi2)
  rprior(pp.prior[[1]])
  rprior(pp.prior[[2]])
  print(pp.prior[[1]])
  
  model <- attributes(object$data)$model
  facs <- names(attr(model, "factors"))
  class(object$data) <- c("data.frame", "list")
  
  if (!is.null(factors)) {
    if (any(is.na(factors))) factors <- facs
    if (!all(factors %in% facs)) 
      stop(paste("Factors argument must contain one or more of:", 
                 paste(facs, collapse=",")))
  }
  
  resp <- names(attr(model, "responses"))
  ns   <- table(object$data[,facs], dnn = facs)
  npar   <- object$n.pars
  nchain <- object$n.chains
  nmc    <- object$nmc
  ntsample <- nchain * nmc
  pnames   <- object$p.names
  # str(samples$theta) ## nchain x npar x nmc
  # str(thetas)        ## (nchainx nmc) x npar
  thetas <- matrix(aperm(object$theta, c(3,1,2)), ncol = npar)
  colnames(thetas) <- pnames
  
  if (is.na(npost)) {
    use <- 1:ntsample 
  } else {
    if (rand) {
      use <- sample(1:ntsample, npost, replace = F) 
    } else {
      use <- round(seq(1, ntsample, length.out = npost))
    }
  }
  
  npost  <- length(use)
  posts   <- thetas[use, ]
  nttrial <- sum(ns) ## number of total trials
  out <- data.frame(matrix(nrow = npost*nttrial, ncol = dim(object$data)[2]))
  
  # for (i in 1:npost) {
  #   tmp <- ggdmc:::simulate_one(model, n = ns, ps = posts[i,], seed = NULL)
  #   out[(1+(i-1)*nttrial):(i*nttrial), names(tmp)] <- tmp
  #   if ( (i %% report) == 0) cat(".")
  # }
  
  ## should replace with parallel
  v <- lapply(1:npost, function(i) {
    ggdmc:::simulate_one(model, n = ns, ps = posts[i,], seed = seed)
  })
  out <- data.table::rbindlist(v)
  # names(out) <- names(object$data)
  reps <- rep(1:npost, each = nttrial)
  out <- cbind(reps, out)
  
  if (!any(is.na(xlim))) {
    out <- out[RT > xlim[1] & RT < xlim[2]] 
  }
  
  attr(out, "data") <- object$data
  return(out)
}

predict_hier0 <- function(samples,n.post=100,probs=c(1:99)/100,bw="nrd0",
         save.simulation=FALSE,factors=NA,save.subject.posts=FALSE,cores=1,ignore.R2=FALSE,
         probs.gglist =c(0.1, 0.5, 0.9), CI.gglist =c(0.025, 0.975),censor=c(NA,NA))
  # apply post.predict to each subject
{
  os <- get.os()
  if ( cores>1 & length(samples)>1 & os=="windows") { 
    cat("No progress indication in multi-core\n")
    require(snowfall,quietly=TRUE)
    require(rlecuyer,quietly=TRUE)
    sfInit(parallel=TRUE, cpus=cores, type="SOCK")
    sfClusterSetupRNG()
    sfLibrary(msm)
    sfLibrary(rtdists)
    sfExportAll()
    out <- sfLapply(samples,post.predict.dmc,n.post=n.post,probs=probs,bw=bw,
                    factors=factors,save.simulation=save.simulation,gglist=TRUE,
                    save.simulation.as.attribute=TRUE,ignore.R2=ignore.R2,censor=censor)
    sfStop() 
  } else if (cores>1) {
    cat("No progress indication in multi-core\n")
    require(parallel, quietly=TRUE)
    out <- mclapply(samples,post.predict.dmc,n.post=n.post,probs=probs,bw=bw,
                    factors=factors,save.simulation=save.simulation,gglist=TRUE,
                    save.simulation.as.attribute=TRUE,ignore.R2=ignore.R2,censor=censor)
  } else {
    out <- lapply(samples,post.predict.dmc,n.post=n.post,probs=probs,bw=bw,
                  factors=factors,save.simulation=save.simulation,gglist=TRUE,
                  save.simulation.as.attribute=TRUE,ignore.R2=ignore.R2,censor=censor)
  }
  
  if ( !save.simulation ) { # Get averages
    sim <- do.call(rbind,lapply(out,function(x){attr(x,"sim")}))
    if ( (any(names(samples[[1]]$data)=="R2")) && !ignore.R2 ) for (i in 1:length(samples)) {
      samples[[i]]$data$R <- 
        paste(as.character(samples[[i]]$data$R),as.character(samples[[i]]$data$R2),sep="") 
      samples[[i]]$data$R[samples[[i]]$data$R2=="DK"] <- "DK"
      samples[[i]]$data$R <- factor(samples[[i]]$data$R)
    }
    dat <- do.call(rbind,lapply(samples,function(x){x$data}))
    facs <- names(attr(attributes(samples[[1]]$data)$model,"factors"))
    if (!is.null(factors)) {
      if (any(is.na(factors))) factors <- facs
      if (!all(factors %in% facs)) 
        stop(paste("Factors argument must contain one or more of:",paste(facs,collapse=",")))
    }
    sim.dqp <- get.dqp(sim[,-1],factors,probs,n.post=1,bw=bw)
    dat.dqp <- get.dqp(sim=dat,factors,probs,bw=bw)
    names(dat.dqp) <- paste("data",names(dat.dqp),sep=".")
    av <- c(sim.dqp,dat.dqp)
    dpqs <- vector(mode="list",length=length(n.post))
    for (i in 1:n.post) {
      simi <- sim[sim$reps==i,-1]
      dpqs[[i]] <- get.dqp(simi,factors,probs,1)
    }
    attr(av,"dpqs") <- dpqs
    # Strip global sim attribute and dpqs for each participant
    out <- lapply(out,function(x){
      attr(x,"sim") <- NULL
      if (!save.subject.posts) attr(x,"dpqs") <- NULL
      x
    })
    attr(av, "gglist") <- get.fitgglist.dmc(sim,dat,factors=factors, noR=FALSE, 
                                            quantiles.to.get = probs.gglist, CI= CI.gglist)
    attr(out,"av") <- av
  }
  out
}