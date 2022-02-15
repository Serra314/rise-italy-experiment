library(ggplot2)
library(viridis)
library(inlabru)
library(INLA)
library(dplyr)
library(data.table)
library(metR)
library(matrixStats)
library(parallel)
library(mvtnorm)
library(MASS)
library(raster)
library(rgeos)
library(mnormt)
library(foreach)
library(gnorm)
library(sf)

## for all the functions below.

## th.p - theta parameters (vector)
## th0.p - theta0 parameters (vector)
## Sigma - 2x2 covariance matrix (matrix)
## Chol.M - Cholensky decomposition of hypothetical Sigma matrix (matrix)
## ------------------------------ ##
## xx - single x-coordinate (scalar)
## yy - single y-coordinate (scalar)
## tt - single time point (scalar)
## ------------------------------ ##
## xs - multiple x-coordinates (vector)
## ys - multiple y-coordinates (vector)
## ts - multiple times (vector)
## ------------------------------ ##
## xh - observed x-coordinates (vector)
## yh - observed y-coordinates (vector)
## th - observed times (vector)
## mh - observed magnitudes (vector)
## ------------------------------ ##
## M0 - magnitude of completeness (scalar)
## T1 - lower time interval extreme (scalar)
## T2 - upper time interval extreme (scalar)
## bdy - Polygon representing the study region (SpatialPolygon)

# create spatial polygon given coordinates
create_poly <- function(coord.matrix){
  x_coords <- c(coord.matrix[,1], coord.matrix[1,1])
  y_coords <- c(coord.matrix[,2], coord.matrix[1,2])
  poly1 <- sp::Polygon(cbind(x_coords,y_coords))
  bdy <- sp::Polygons(list(poly1), ID = "A")
  bdy <- sp::SpatialPolygons(list(bdy))
  as(bdy, 'SpatialPolygonsDataFrame')
}


# magnitude triggering function (vector mh)
g.m <- function(th.p, mh, M0){
  #expm1( ( expm1(th.p[3]) + 1 ) * (mh - M0) ) + 1
  exp(exp(th.p[3])*(mh - M0))
}

# time triggering function (scalar tt and vector th)
g.t <- function(th.p, tt, th){
  output <- rep(0, length(th))
  idx.0 <- th >= tt
  
  gamma.h <- (tt - th[!idx.0])/exp(th.p[4])
  output[!idx.0] <- expm1( (-1-exp(th.p[5]))*log(gamma.h + 1) ) + 1
  output
}


## FUNCTIONS TO GET SIGMA FROM THETA BASED ON BARTLETT'S DECOMPOSITION
# function to get the values of the diagonal of A
#   mandatory input : th --> value of parameter theta 
#
#   optional input : df.p --> degrees of freedom (default set to 1)
# 
#   output : scalar, value of a, should be a Chi-square with df.p degrees of freedom 
# 
#   notes : it assumes th comes from a normal 0, 1 distribution
get_a <- function(th, df.p = 1) {
  th_pos <- th >= 0
  res <- th
  res[!th_pos] <- qchisq(pnorm(th[!th_pos], log.p = TRUE), df = df.p, log.p = TRUE)
  res[th_pos] <- qchisq(pnorm(th[th_pos], log.p = TRUE, lower.tail = FALSE), df = df.p, log.p = TRUE,
                        lower.tail = FALSE) 
  res
}

# function to create Sigma from values of theta applying the Bartlett's decomposition
# function to get the values of the diagonal of A
#   mandatory input : theta6/7/8 --> value of parameter theta 
#                   : Chol.M --> Cholensky decomposition of a covariance matrix Sigma' 
# 
#   output : 2x2 Covariance matrix 
# 
#   notes : it assumes theta6/7/8 comes from a normal 0, 1 distribution.
#         : Chol.M is the Cholensky decomposition of Sigma' where Sigma' has to be seen as a prior 
#           covariance matrix

Sigma.from.theta <- function(theta6, theta7, theta8, Chol.M){
  # create A
  A = matrix(c(get_a(theta6, df.p = 2), 0,
               theta8, get_a(theta7, df.p = 1)), byrow = TRUE, ncol = 2)
  # calculate Sigma
  temp = Chol.M %*% A
  Sigma =  t(temp) %*% temp
  Sigma              
}


## alternative
g.s <- function(th.p, xx, yy, xh, yh, Sigma = NULL, Chol.M = NULL){
  if(length(xh) != length(yh)){
    stop('xh and yh has to be of the same length')
  }
  if(is.null(Sigma)){
    if(is.null(Chol.M)){
      stop('Please provide either Sigma or Chol.M')  
    }
    else{
      Sigma <- Sigma.from.theta(th.p[6], th.p[7], th.p[8], Chol.M)
    }
  }
  #mean.h <- cbind(xh, yh)
  #sapply(1:nrow(mean.h), \(idx) dmvnorm(c(xx,yy), mean = mean.h[idx,], sigma = Sigma))
   mean.dif <- cbind(xx -xh, yy - yh)
   det.s <- det(Sigma)
   S.inv <- pd.solve(Sigma)
   (1/(2*pi))*((det(Sigma))^(-1/2))*exp((-1/2)*diag( mean.dif %*% S.inv %*% t(mean.dif) ) )
   
}



# triggering function.
g.x <- function(th.p, tt, xx, yy, th, xh, yh, mh, M0, Sigma = NULL, Chol.M = NULL){
   log.output <- th.p[2] + log1p(g.m(th.p, mh, M0)-1) + log1p(g.t(th.p, tt, th)-1) + log1p(g.s(th.p, xx, yy, xh, yh, Sigma, Chol.M)-1)
   expm1(log.output) + 1
}

## conditional intensity
lambda <- function(th.p, tt, xx, yy, th, xh, yh, mh, M0, Sigma = NULL, Chol.M = NULL){
  expm1(th.p[1]) + 1 + sum(g.x(th.p, tt, xx, yy, th, xh, yh, mh, M0, Sigma, Chol.M))
}


lambda_deb <- function(th.p, tt, xx, yy, th, xh, yh, mh, M0, Sigma = NULL, Chol.M = NULL){
  list(background = exp(th.p[1]), 
       K = exp(th.p[2]),
       gm = g.m(th.p, mh, M0),
       gt = g.t(th.p, tt, th),
       gs = g.s(th.p, xx, yy, xh, yh, Sigma = Sigma))
}



## integrated intensity for single point
logLambda.h.vec <- function(th.p, th, xh, yh, mh, M0, T1, T2, bdy, Sigma = NULL, Chol.M = NULL){
  if(is.null(Sigma)){
    if(is.null(Chol.M)){
      stop('Please provide either Sigma or Chol.M')  
    }
    else{
      Sigma <- Sigma.from.theta(th.p[6], th.p[7], th.p[8], Chol.M)
    }
  }
  #print(th)
  if(any(th > T2)){
    return('Error: th > T2')
  }
  
  # integral in space
  lower_ <- c(bdy@bbox[1,1], bdy@bbox[2,1])
  upper_ <- c(bdy@bbox[1,2], bdy@bbox[2,2])
  p.loc <- cbind(xh, yh)
  I.s <- sapply(1:nrow(p.loc), \(idx)
                abs(pmvnorm(lower = lower_, upper = upper_, 
                        mean = p.loc[idx,], sigma = Sigma,
                        keepAttr = FALSE)))
  
  # integral in time
  Tl <- sapply(th, \(x) max(x, T1))
  gamma.l <- (Tl - th)/exp(th.p[4])
  gamma.u <- (T2 - th)/exp(th.p[4])
  w.l <- (gamma.l + 1)^(-exp(th.p[5]))  
  w.u <- (gamma.u + 1)^(-exp(th.p[5]))
  
  # output
  th.p[2] + exp(th.p[3])*(mh - M0) + th.p[4] - th.p[5] + log1p(w.l - 1) + log1p(-w.u/w.l) + log1p(I.s - 1)
}

ETAS.loglik <- function(th.p, th, xh, yh, mh, M0, T1, T2, bdy, Sigma = NULL, Chol.M = NULL, 
                         prec.prior = rep(1,5)){
  
  sd.prior <- sqrt(1/prec.prior)
  logIntegral.trig <- logLambda.h.vec(th.p, th, xh, yh, mh, M0, T1, T2, bdy, Sigma,
                                      Chol.M)
  Integral. <- exp(th.p[1])*(T2 - T1) + exp(logSumExp(logIntegral.trig))
  log.lambdas <- unlist(mclapply(1:length(th), \(ii) 
                        log1p(lambda(th.p, th[ii], xh[ii], yh[ii],
                                     th, xh, yh, mh, M0, Sigma, Chol.M) - 1), 
                        mc.cores = 5))
  
  - Integral. + sum(log.lambdas) +
    sum(dnorm(th.p[1:5], sd = sd.prior, log = TRUE))
}


## integrated time-triggering function
It <- function(th.p, th, T2){
  gamma.u <- (T2 - th)/exp(th.p[4])
  exp(th.p[4] - th.p[5])*(1 - (gamma.u + 1)^(-exp(th.p[5])))
}
## inverse of integrated time-triggering function
Inv.It <- function(th.p, omega, th){
  th + exp(th.p[4])*( (1 - omega*exp(-th.p[4] + th.p[5]))^(-1/exp(th.p[5])) - 1)
}

## sampling times
sample.time <- function(th.p, n.ev, th, T2){
  if(n.ev == 0){
    df <- data.frame(ts = 1, x = 1, y = 1, mags = 1, gen = 0)
    return(df[-1,])
    }
  bound.l <- 0 #It(th.p, th, T)
  bound.u <- It(th.p, th, T2)
  unif.s <- runif(n.ev, min = bound.l, max = bound.u)
  t.sample <- Inv.It(th.p, unif.s, th)
  t.sample
}

# sampling locations
sample.loc <- function(th.p, xh, yh, n.ev, bdy, Sigma = NULL, Chol.M = NULL, crsobj = NULL){
  
  if(is.null(Sigma)){
    if(is.null(Chol.M)){
      stop('Please provide either Sigma or Chol.M')  
    }
    else{
      Sigma <- Sigma.from.theta(th.p[6], th.p[7], th.p[8], Chol.M)
    }
  }
  
  # initialize empty SpatialPoints
  samp.points <- SpatialPoints(data.frame(x = 0, y = 0))
  samp.points <- samp.points[-1,]
  if(!is.null(crsobj)){
    proj4string(samp.points) <- crsobj
  }
  # initialize number of placed events
  num <- 0
  
  # until we placed all the events
  while (num < n.ev){
    # sample from the 2D Gaussian kernel without boundaries
    pts.matrix <- rmvnorm(n.ev, mean = c(xh, yh), sigma = Sigma)
    # transform the sample in SpatialPoints and project
    pts <- data.frame(x = pts.matrix[,1], y = pts.matrix[,2])
    coordinates(pts) <- ~ x + y
    if(!is.null(crsobj)){
      proj4string(pts) <- crsobj
      pts <- spTransform(pts, crsobj)
    }
    # discard the ones outside bdy
    pts <- crop(pts, bdy)
    # merge sampled points
    if(length(pts) > 0){
      samp.points <- rbind(samp.points, pts)
      num <- length(samp.points)
    }
    else{
      #print('no retained')
    }
  }
  
  ## if we retained a number of points > n.ev, we select n.ev events at random
  kp <- sample(seq(1, length(samp.points), by=1), n.ev, replace=FALSE)
  samp.points <- samp.points[kp,]
  samp.points
}


# sampling triggered events
sample.triggered <- function(th.p, beta.p, th, xh, yh, n.ev, M0, T1, T2, bdy, 
                             Sigma = NULL, Chol.M = NULL, crsobj = NULL){
  # if the number of events to be placed is zero returns an empty data.frame
  if(n.ev == 0){
    samp.points <- data.frame(x = 1, y = 1, ts = 1, mags = 1)
    samp.points <- samp.points[-1,]
    return(samp.points)
  }
  else{
    
    # sample times
    samp.ts <- sample.time(th.p, n.ev, th, T2)
    # sample magnitudes
    samp.mags <- rexp(n.ev, rate = beta.p) + M0
    # sample locations
    samp.locs <- sample.loc(th.p, xh, yh, n.ev, bdy, Sigma, Chol.M, crsobj)
    
    # build output dataset
    samp.points <- data.frame(ts = samp.ts, mags = samp.mags, x = samp.locs@coords[,1],
                              y = samp.locs@coords[,2])
    # return only the ones with time different from NA (the one with NA are outside the interval T1, T2)
    # even though it should not happen given how we built sample.omori
    return(samp.points[!is.na(samp.points$ts),])
  }
}


sample.generation <- function(th.p, beta.p, Ht, M0, T1, T2, bdy,
                              Sigma = NULL, Chol.M = NULL, crsobj = NULL, ncore = 2){
  # number of parents
  n.parent <- nrow(Ht)
  # calculate the aftershock rate for each parent in history
  trig.rates <- exp(logLambda.h.vec(th.p, Ht$ts, Ht$x, Ht$y, Ht$mags, M0, T1, T2, bdy,
                                 Sigma, Chol.M))
  # extract number of aftershock for each parent
  n.ev.v <- sapply(1:n.parent, function(x) rpois(1, trig.rates[x]))
  
  # if no aftershock has to be generated returns empty data.frame
  if(sum(n.ev.v) == 0){
    app <- data.frame(x = 1, y = 1, ts = 1, mags = 1)
    app <- app[-1,]
    return(app)
  }
  
  # identify parent with number of aftershocks > 0 
  idx.p <- which(n.ev.v > 0)

  #print(sample.triggered(theta.v, beta.p, Sigma, Chol.M, n.ev.v[idx.p[1]], Ht[idx.p[1],], T1, T2, M0, bdy, crsobj))
  # sample (in parallel) the aftershocks for each parent 
  sample.list <- lapply(idx.p, function(idx) 
    sample.triggered(th.p, beta.p, Ht$ts[idx], Ht$x[idx], Ht$y[idx], n.ev.v[idx], M0, T1, T2, bdy, 
                     Sigma, Chol.M, crsobj))#,
    #mc.cores = ncore)
  
  # bind the data.frame in the list and return
  sample.pts <- bind_rows(sample.list) 
  sample.pts
}


sample.ETAS <- function(th.p, beta.p, M0, T1, T2, bdy, Sigma = NULL, Chol.M = NULL, crsobj = NULL, Ht = NULL, ncore = 2, 
                        bk.field.list = NULL, debugg = TRUE){
  # if the upper extreme greater than lower
  if(T2 < T1){
    stop('Error - right-end of time interval greater than left-end')
  }
  
  if(debugg){
    n.bkg <- rpois(1, exp(th.p[1])*(T2 - T1))
  }
  else{
    if(is.null(crsobj)){
      n.bkg <- rpois(1, exp(th.p[1])*(T2 - T1)*(gArea(bdy)))
    } else{
      # background number of events (area calculated in squared km)
      n.bkg <- rpois(1, exp(th.p[1])*(T2 - T1)*(area(bdy)/1000000))
    }
  }
  
  # if no background events are generated initialize an empty data.frame
  if(n.bkg == 0){
    bkg.df <- data.frame(x = 1, y = 1, ts = 1, mags = 1)
    bkg.df <- bkg.df[-1,]
  }
  else{
    # sample bkg events
    # if no bk.field.list element is passed it assumes uniform background rate
    if(is.null(bk.field.list)){
      bkg.locs <- spsample(bdy, n.bkg, 'random', iter = 10)
      if(!is.null(crsobj)){
        proj4string(bkg.locs) <- crsobj
        bkg.locs <- spTransform(bkg.locs, crsobj)
      }
      bkg.df <- data.frame(x = bkg.locs@coords[,1], 
                           y = bkg.locs@coords[,2], 
                           ts = runif(n.bkg, T1, T2), 
                           mags = rexp(n.bkg, beta.p) + M0, 
                           gen = 1)
    }
    # otherwise it samples using the information provided
    else{
      bkg.locs <- point_sampler(bk.field.list$loglambda, 
                                bdy, bk.field.list$mesh, n.bkg, crsobj)
      
      bkg.df <- data.frame(x = bkg.locs@coords[,1], 
                           y = bkg.locs@coords[,2], 
                           ts = runif(n.bkg, T1, T2), 
                           mags = rexp(n.bkg, beta.p) + M0, 
                           gen = 1)
    }
  }
  
  # if known events are provided
  if(!is.null(Ht)){
    #### TODO : the function has to generate all the points.
    # sample a generation from the known events
    gen.from.past <- sample.generation(th.p, beta.p, Ht, M0, T1, T2,  bdy, Sigma, Chol.M, crsobj, ncore)
    # if at least an aftershock is produced
    if(nrow(gen.from.past) > 0){
      # set generation
      gen.from.past$gen = 0
      # Merge first generation and background events
      Gen.list <- list(rbind(gen.from.past, bkg.df))  
    }
    else{
      Gen.list <- list(bkg.df)
    }
    
  }
  else{
    Gen.list <- list(bkg.df)
  }
  # stop if we have no background events and no events generated from known observations
  if(nrow(Gen.list[[1]]) == 0){
    #print(exp(theta.v[1])*(T2 - T1)*(area(bdy)/1000000))
    #stop('No events generated - increase theta1')
    return(Gen.list)
  }
  
  # initialize flag and gen counter
  flag = TRUE
  gen = 1
  # this goes until the condition inside the loop is met
  while(flag){
    # set parents
    parents <- Gen.list[[gen]]
    #print(c(T1,T2))
    #print(range(parents$ts))
    # generate aftershocks
    triggered <- sample.generation(th.p, beta.p, parents, 
                                   M0, T1, T2, bdy, Sigma, Chol.M, crsobj, ncore)
    #print(nrow(triggered))
    # stop the loop if there are no more aftershocks
    if(nrow(triggered) == 0){
      flag = FALSE}
    else{
      # set generations
      triggered$gen = gen + 1
      # store new generation
      Gen.list[[gen + 1]] = triggered
      # update generation counter
      gen = gen + 1
    }
  }
  Gen.list
}


### functions for INLA
ETAS.fit_temporal <- function(sample.s, N.breaks.min, max.length, M0, T1, T2, bdy, 
                          Sigma = NULL, Chol.M = NULL,
                          prior.mean = rep(0,5), prior.prec = rep(1,5),
                          bru.opt = list(bru_verbose = 3,
                                         bru_max_iter = 50)){
  
  # create different poisson counts models
  # first for background
  df.0 <- data.frame(x = 0.5, y = 0.5, ts = 0.5, mags = 2.5, counts = 0, exposures = 1)
  # this is the expression of log(Lambda0)
  form.0 <- counts ~ th1 + log(T2 - T1)
  
  # first likelihood
  lik.0 <- inlabru::like(formula = form.0,
                         data = df.0,
                         family = 'poisson',
                         options = list(E = df.0$exposures))
  
  
  ### create time bins 
  
  df.j <- foreach(idx = 1:nrow(sample.s), .combine = rbind) %do% {
    
    tt <- sample.s$ts[idx]
    N.breaks <- max(N.breaks.min, ceiling((T2 - tt)/max.length))
    
    kk <- seq_len(N.breaks) - 1
    Time.bins.matrix <- tt + cbind(kk, kk + 1) / N.breaks * (T2 - tt)
    
    data.frame(x = rep(sample.s$x[idx], each = nrow(Time.bins.matrix)),
               y = rep(sample.s$y[idx], each = nrow(Time.bins.matrix)),
               ts = rep(sample.s$ts[idx], each = nrow(Time.bins.matrix)),
               mags = rep(sample.s$mags[idx], each = nrow(Time.bins.matrix)),
               counts = 0,
               exposures = 1,
               bin.start = Time.bins.matrix[,1],
               bin.end = Time.bins.matrix[,2])
    
  }
  
  
  logLambda.h.inla <- function(th1, th2, th3, th4, th5, ts, x, y, mags, M0, T1.v, T2.v, bdy, 
                               Sigma, Chol.M){
    th.p <- c(th1[1], th2[1], th3[1], th4[1], th5[1])
    unlist(mclapply(1:length(ts), \(idx)
                    logLambda.h.vec(th.p, ts[idx], x[idx], y[idx], mags[idx], M0, T1.v[idx], T2.v[idx], bdy, 
                                    Sigma, Chol.M),
                    mc.cores = 5))
  }
  # creating formula for past events contributions to integrated lambda
  form.j.part <- counts ~ logLambda.h.inla(th1, 
                                           th2, 
                                           th3, 
                                           th4, 
                                           th5,
                                           ts, x, y, mags, M0, bin.start, bin.end, bdy, 
                                           Sigma, Chol.M)
  # second for triggered part of the integral
  lik.j.part <- inlabru::like(formula = form.j.part,
                              data = df.j,
                              family = 'poisson',
                              options = list(E = df.j$exposures)) 
  
  # third is for the sum of the log intensities
  df.s <- data.frame(x = sample.s$x, y = sample.s$y, ts = sample.s$ts, mags = sample.s$mags, 
                     counts = 1, exposures = 0)
  
  loglambda.inla <- function(th1, th2, th3, th4, th5, tt, xx, yy, 
                             th, xh, yh, mh, M0, Sigma = NULL, Chol.M = NULL){
    th.p <- c(th1[1], th2[1], th3[1], th4[1], th5[1])
    
    # print('loglambda')
    # print(th.p)
    # print(Sigma)
    unlist(mclapply(1:length(tt), \(idx) 
                    log(lambda(th.p, tt[idx], xx[idx], yy[idx], 
                               th, xh, yh, mh, M0, Sigma, Chol.M)),
                    mc.cores = 5))
  }
  # creating formula for summation part
  form.s.part <- counts ~ loglambda.inla(th1, 
                                         th2, 
                                         th3, 
                                         th4, 
                                         th5, ts, x, y,
                                         sample.s$ts, sample.s$x, 
                                         sample.s$y, sample.s$mags,
                                         M0, Sigma, Chol.M)
  lik.s.part <- inlabru::like(formula = form.s.part,
                              data = df.s,
                              family = 'poisson',
                              options = list(E = df.s$exposures)) 
  cmp.part <- counts ~ -1 + 
    th1(1, model = 'linear', mean.linear = prior.mean[1], prec.linear = prior.prec[1]) + 
    th2(1, model = 'linear', mean.linear = prior.mean[2], prec.linear = prior.prec[2]) +
    th3(1, model = 'linear', mean.linear = prior.mean[3], prec.linear = prior.prec[3]) +
    th4(1, model = 'linear', mean.linear = prior.mean[4], prec.linear = prior.prec[4]) +
    th5(1, model = 'linear', mean.linear = prior.mean[5], prec.linear = prior.prec[5]) 
  
  bru(lik.0, lik.j.part, lik.s.part, components = cmp.part, 
      options = bru.opt)
}


# special case with isotropic spatial triggering function
ETAS.fit <- function(sample.s, N.breaks.min, max.length, M0, T1, T2, bdy, 
                           prior.mean = rep(0,6), prior.prec = rep(1,6),
                           bru.opt = list(bru_verbose = 3,
                                          bru_max_iter = 50)){
  

  # create different poisson counts models
  # first for background
  df.0 <- data.frame(x = 0.5, y = 0.5, ts = 0.5, mags = 2.5, counts = 0, exposures = 1)
  # this is the expression of log(Lambda0)
  form.0 <- counts ~ th1 + log(T2 - T1)
  
  # first likelihood
  lik.0 <- inlabru::like(formula = form.0,
                         data = df.0,
                         family = 'poisson',
                         options = list(E = df.0$exposures))
  
  
  ### create time bins 
  
  df.j <- foreach(idx = 1:nrow(sample.s), .combine = rbind) %do% {
    
    tt <- sample.s$ts[idx]
    N.breaks <- max(N.breaks.min, ceiling((T2 - tt)/max.length))
    
    kk <- seq_len(N.breaks) - 1
    Time.bins.matrix <- tt + cbind(kk, kk + 1) / N.breaks * (T2 - tt)
    
    data.frame(x = rep(sample.s$x[idx], each = nrow(Time.bins.matrix)),
               y = rep(sample.s$y[idx], each = nrow(Time.bins.matrix)),
               ts = rep(sample.s$ts[idx], each = nrow(Time.bins.matrix)),
               mags = rep(sample.s$mags[idx], each = nrow(Time.bins.matrix)),
               counts = 0,
               exposures = 1,
               bin.start = Time.bins.matrix[,1],
               bin.end = Time.bins.matrix[,2])
    
  }
  
  
  logLambda.h.inla <- function(th1, th2, th3, th4, th5, th6, 
                               ts, x, y, mags, M0, T1.v, T2.v, bdy){
    th.p <- c(th1[1], th2[1], th3[1], th4[1], th5[1])
    Sigma_ <- matrix(c(exp(th6[1]), 0, 0, exp(th6[1])), byrow = TRUE, ncol = 2)
    # print('logLambda')
    # print(th.p)
    # print(Sigma)
    unlist(mclapply(1:length(ts), \(idx)
           logLambda.h.vec(th.p, ts[idx], x[idx], y[idx], mags[idx], M0, T1.v[idx], T2.v[idx], bdy, 
                           Sigma = Sigma_, Chol.M = NULL),
           mc.cores = 5))
  }
  # creating formula for past events contributions to integrated lambda
  form.j.part <- counts ~ logLambda.h.inla(th1, th2, th3, th4, th5, th6,
                                           ts, x, y, mags, M0, bin.start, bin.end, bdy)
  # second for triggered part of the integral
  lik.j.part <- inlabru::like(formula = form.j.part,
                              data = df.j,
                              family = 'poisson',
                              options = list(E = df.j$exposures)) 
  
  # third is for the sum of the log intensities
  df.s <- data.frame(x = sample.s$x, y = sample.s$y, ts = sample.s$ts, mags = sample.s$mags, 
                     counts = 1, exposures = 0)
  
  loglambda.inla <- function(th1, th2, th3, th4, th5, th6, tt, xx, yy, 
                             th, xh, yh, mh, M0){
    th.p <- c(th1[1], th2[1], th3[1], th4[1], th5[1])
    Sigma_ <- matrix(c(exp(th6[1]), 0, 0, exp(th6[1])), byrow = TRUE, ncol = 2)
    # print('loglambda')
    # print(th.p)
    # print(Sigma)
    unlist(mclapply(1:length(tt), \(idx) 
           log(lambda(th.p, tt[idx], xx[idx], yy[idx], 
                      th, xh, yh, mh, M0, Sigma = Sigma_, Chol.M = NULL)),
           mc.cores = 5))
  }
  # creating formula for summation part
  form.s.part <- counts ~ loglambda.inla(th1, th2, th3, th4, th5, th6, 
                                         ts, x, y,
                                         sample.s$ts, sample.s$x, 
                                         sample.s$y, sample.s$mags, M0)
  lik.s.part <- inlabru::like(formula = form.s.part,
                              data = df.s,
                              family = 'poisson',
                              options = list(E = df.s$exposures)) 
  cmp.part <- counts ~ -1 + 
    th1(1, model = 'linear', mean.linear = prior.mean[1], prec.linear = prior.prec[1]) + 
    th2(1, model = 'linear', mean.linear = prior.mean[2], prec.linear = prior.prec[2]) +
    th3(1, model = 'linear', mean.linear = prior.mean[3], prec.linear = prior.prec[3]) +
    th4(1, model = 'linear', mean.linear = prior.mean[4], prec.linear = prior.prec[4]) +
    th5(1, model = 'linear', mean.linear = prior.mean[5], prec.linear = prior.prec[5]) +
    th6(1, model = 'linear', mean.linear = prior.mean[6], prec.linear = prior.prec[6])
  
  bru(lik.0, lik.j.part, lik.s.part, components = cmp.part, 
      options = bru.opt)
}


CalcG.dx <- function(x_, ..., eps = .Machine$double.eps^(1/3)) {
  (qgnorm(pnorm(x_ + eps, log.p = TRUE), ..., log.p = TRUE) - 
     qgnorm(pnorm(x_ - eps, log.p = TRUE), ..., log.p = TRUE)) / (2*eps)
}

CalcG.sx <- function(x_, ..., eps = .Machine$double.eps^(1/3)) {
  (qgnorm(pnorm(x_ + eps, log.p = TRUE, lower.tail = TRUE), ..., log.p = TRUE, lower.tail = TRUE) - 
     qgnorm(pnorm(x_ - eps, log.p = TRUE, lower.tail = TRUE), ..., log.p = TRUE, lower.tail = TRUE)) / (2*eps)
}


bru_ft <- function (qfun, x, ..., tail.split. = 0) 
{
  if (is.null(tail.split.)) {
    upper <- x >= 0
  }
  else {
    upper <- x >= tail.split.
  }
  res <- numeric(length(x))
  if (sum(upper) > 0) {
    to.app <- x[upper] > 7
    if(sum(to.app) > 0){
      f0 <- qfun(pnorm(7, log.p = TRUE), ..., log.p = TRUE) 
      f0.der <- CalcG.dx(7, ...) 
      print(f0.der)
      res[upper][to.app] <- f0 + (x[upper][to.app] - 7)*f0.der
    }
    if(sum(!to.app) > 0){
      res[upper][!to.app] <-  qfun(pnorm(x[upper][!to.app], 
                                         log.p = TRUE), ..., log.p = TRUE)
    }
  }
  # this works fine for N(0,1)
  if (sum(!upper) > 0) {
    to.app <- x[!upper] < -7
    if(sum(to.app) > 0){
      f0 <- qfun(pnorm(-7, log.p = TRUE, lower.tail = TRUE), ..., lower.tail = TRUE, log.p = TRUE) 
      f0.der <- CalcG.sx(-7, ...) 
      res[!upper][to.app] <- f0 + (x[!upper][to.app] + 7)*f0.der
    }
    if(sum(!to.app) > 0){
      res[!upper][!to.app] <-  qfun(pnorm(x[!upper][!to.app], 
                                          log.p = TRUE, lower.tail = TRUE), ..., log.p = TRUE, lower.tail = TRUE)
    }
  }
  res
}


th_tranf <- function(th., mu., alpha., beta.){
  bru_ft(qgnorm, th., mu., 
         alpha., beta.)
}



# special case with isotropic spatial triggering function
ETAS.fit_gnprior <- function(sample.s, N.breaks.min, max.length, M0, T1, T2, bdy, 
                             prior.mu = rep(0,6), 
                             prior.alpha = rep(1,6),
                             prior.beta = rep(1,6),
                             bru.opt = list(bru_verbose = 3,
                             bru_max_iter = 50)){
  

  # create different poisson counts models
  # first for background
  df.0 <- data.frame(x = 0.5, y = 0.5, ts = 0.5, mags = 2.5, counts = 0, exposures = 1)
  # this is the expression of log(Lambda0)
  form.0 <- counts ~ th_tranf(th1, prior.mu[1], prior.alpha[1], prior.beta[1]) + log(T2 - T1)
  
  # first likelihood
  lik.0 <- inlabru::like(formula = form.0,
                         data = df.0,
                         family = 'poisson',
                         options = list(E = df.0$exposures))
  
  
  ### create time bins 
  
  df.j <- foreach(idx = 1:nrow(sample.s), .combine = rbind) %do% {
    
    tt <- sample.s$ts[idx]
    N.breaks <- max(N.breaks.min, ceiling((T2 - tt)/max.length))
    
    kk <- seq_len(N.breaks) - 1
    Time.bins.matrix <- tt + cbind(kk, kk + 1) / N.breaks * (T2 - tt)
    
    data.frame(x = rep(sample.s$x[idx], each = nrow(Time.bins.matrix)),
               y = rep(sample.s$y[idx], each = nrow(Time.bins.matrix)),
               ts = rep(sample.s$ts[idx], each = nrow(Time.bins.matrix)),
               mags = rep(sample.s$mags[idx], each = nrow(Time.bins.matrix)),
               counts = 0,
               exposures = 1,
               bin.start = Time.bins.matrix[,1],
               bin.end = Time.bins.matrix[,2])
    
  }
  
  
  logLambda.h.inla <- function(th1, th2, th3, th4, th5, th6, 
                               ts, x, y, mags, M0, T1.v, T2.v, bdy){
    th.p <- c(th1[1], th2[1], th3[1], th4[1], th5[1])
    Sigma_ <- matrix(c(exp(th6[1]), 0, 0, exp(th6[1])), byrow = TRUE, ncol = 2)
    # print('logLambda')
    #print(th.p)
    # print(th6[1])
    # print(Sigma_)
    unlist(mclapply(1:length(ts), \(idx)
                    logLambda.h.vec(th.p, ts[idx], x[idx], y[idx], mags[idx], M0, T1.v[idx], T2.v[idx], bdy, 
                                    Sigma = Sigma_, Chol.M = NULL),
                    mc.cores = 5))
  }
  # creating formula for past events contributions to integrated lambda
  form.j.part <- counts ~ logLambda.h.inla(th_tranf(th1, prior.mu[1], prior.alpha[1], prior.beta[1]), 
                                           th_tranf(th2, prior.mu[2], prior.alpha[2], prior.beta[2]),
                                           th_tranf(th3, prior.mu[3], prior.alpha[3], prior.beta[3]), 
                                           th_tranf(th4, prior.mu[4], prior.alpha[4], prior.beta[4]), 
                                           th_tranf(th5, prior.mu[5], prior.alpha[5], prior.beta[5]), 
                                           th_tranf(th6, prior.mu[6], prior.alpha[6], prior.beta[6]),
                                           ts, x, y, mags, M0, bin.start, bin.end, bdy)
  # second for triggered part of the integral
  lik.j.part <- inlabru::like(formula = form.j.part,
                              data = df.j,
                              family = 'poisson',
                              options = list(E = df.j$exposures)) 
  
  # third is for the sum of the log intensities
  df.s <- data.frame(x = sample.s$x, y = sample.s$y, ts = sample.s$ts, mags = sample.s$mags, 
                     counts = 1, exposures = 0)
  
  loglambda.inla <- function(th1, th2, th3, th4, th5, th6, tt, xx, yy, 
                             th, xh, yh, mh, M0){
    th.p <- c(th1[1], th2[1], th3[1], th4[1], th5[1])
    Sigma_ <- matrix(c(exp(th6[1]), 0, 0, exp(th6[1])), byrow = TRUE, ncol = 2)
    # print('loglambda')
    # print(th.p)
    # print(Sigma)
    unlist(mclapply(1:length(tt), \(idx) 
                    log(lambda(th.p, tt[idx], xx[idx], yy[idx], 
                               th, xh, yh, mh, M0, Sigma = Sigma_, Chol.M = NULL)),
                    mc.cores = 5))
  }
  # creating formula for summation part
  form.s.part <- counts ~ loglambda.inla(th_tranf(th1, prior.mu[1], prior.alpha[1], prior.beta[1]), 
                                         th_tranf(th2, prior.mu[2], prior.alpha[2], prior.beta[2]),
                                         th_tranf(th3, prior.mu[3], prior.alpha[3], prior.beta[3]), 
                                         th_tranf(th4, prior.mu[4], prior.alpha[4], prior.beta[4]), 
                                         th_tranf(th5, prior.mu[5], prior.alpha[5], prior.beta[5]), 
                                         th_tranf(th6, prior.mu[6], prior.alpha[6], prior.beta[6]), 
                                         ts, x, y,
                                         sample.s$ts, sample.s$x, 
                                         sample.s$y, sample.s$mags, M0)
  lik.s.part <- inlabru::like(formula = form.s.part,
                              data = df.s,
                              family = 'poisson',
                              options = list(E = df.s$exposures)) 
  cmp.part <- counts ~ -1 + 
    th1(1, model = 'linear', mean.linear = 0, prec.linear = 1) + 
    th2(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th3(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th4(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th5(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th6(1, model = 'linear', mean.linear = 0, prec.linear = 1)
  
  bru(lik.0, lik.j.part, lik.s.part, components = cmp.part, 
      options = bru.opt)
}
