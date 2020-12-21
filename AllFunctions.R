## Packages ##
  library(lars)
  library(tempdisagg)

## Storing results for different scenarios of AR parameter rho and dimension p ##
  # Functions involved: SparseTD.simulation()
  coef3 <- rhofinal3 <- lambdafinal3 <- series3 <- cl.coef3 <- cl.series3 <- cl.rho3 <- trueseries3 <- list()
  coef2 <- rhofinal2 <- lambdafinal2 <- series2 <- cl.coef2 <- cl.series2 <- cl.rho2 <- trueseries2 <- list()
  rfcoef2 <- rfrhofinal2 <- rflambdafinal2 <- rfseries2 <- list()
  rfcoef3 <- rfrhofinal3 <- rflambdafinal3 <- rfseries3 <- list()
  pvector <- c(30,90,150) # dimensions
  rhovector <- c(0.2,0.5,0.8) # rho parameter
  N <- 1 # number of experiments
  n <- 100 # low frequency sample size
  k <- 10 # number of non-zero betas
  sig <- 1 # variance of disturbance
  
  
  start.time = Sys.time()
  for(rhoref in 1:3){
    for(pref in 1:3) {
      sim1 <- SparseTD.simulation(N,n,pvector[pref],k,rhovector[rhoref],sig)
      
      coef2[[pref]] <- sim1$coefficients
      rhofinal2[[pref]] <- sim1$rho
      series2[[pref]] <- sim1$series
      lambdafinal2[[pref]] <- sim1$lambda
      
      rfcoef2[[pref]] <- sim1$rfcoefficients
      rfrhofinal2[[pref]] <- sim1$rfrho
      rfseries2[[pref]] <- sim1$rfseries
      rflambdafinal2[[pref]] <- sim1$rflambda
      
      cl.coef2[[pref]] <- sim1$cl.coef
      cl.series2[[pref]] <- sim1$cl.series
      cl.rho2[[pref]] <- sim1$cl.rho
      trueseries2[[pref]] <- sim1$true.series
    }
    coef3[[rhoref]] <- coef2
    rhofinal3[[rhoref]] <- rhofinal2
    series3[[rhoref]] <- series2
    lambdafinal3[[rhoref]] <- lambdafinal2
    
    rfcoef3[[rhoref]] <- rfcoef2
    rfrhofinal3[[rhoref]] <- rfrhofinal2
    rfseries3[[rhoref]] <- rfseries2
    rflambdafinal3[[rhoref]] <- rflambdafinal2
    
    cl.coef3[[rhoref]] <- cl.coef2
    cl.series3[[rhoref]] <- cl.series2
    cl.rho3[[rhoref]] <- cl.rho2
    trueseries3[[rhoref]] <- trueseries2
  }
  end.time <- Sys.time()
  
###############################################################################################################
  #### FUNCTIONS #####
  
## Main simulation loop ##
  # Functions involved: Y.sim(), SparseTD.estimates(), chowlin(), agg.mat()
  SparseTD.simulation <- function(N,n,p,k,truerho,truevariance) {
    
    m=4*n; 
    phi = rep(0,p) # CHANGE THIS TO rep(1,p) FOR RANDOM WALK DESIGN!!
    sig = rep(1,p) 
    beta <- c(rep(5,k), rep(0,p-k))
    rho.final <- lambda.final <- cl.rho <- rep(NA,N) 
    bhat.final <- series.final <- cl.coef <- cl.series <- trueseries <- list()
    rfrho.final <- rflambda.final <- rep(NA,N)
    rfbhat.final <- rfseries.final <- list()
    
    for(trial in 1:N) {
      set.seed(trial)
      data <- Y.sim(m, beta, phi, sig, truerho, truevariance)
      X <- data$X
      z.true <- data$Y
      C <- agg.mat(n,4) # CHANGE THIS TO 3 FOR QUARTERLY TO MONTHLY
      y <- C%*%z.true
      
      optimiser <- SparseTD.estimates(X,y)
      
      rho.final[trial] <- optimiser$rho.estimate
      bhat.final[[trial]] <- optimiser$beta.estimate
      series.final[[trial]] <- optimiser$series.estimate
      lambda.final[trial] <- optimiser$lambda.estimate
      
      rfrho.final[trial] <- optimiser$rfrho.estimate
      rfbhat.final[[trial]] <- optimiser$rfbeta.estimate
      rfseries.final[[trial]] <- optimiser$rfseries.estimate
      rflambda.final[trial] <- optimiser$rflambda.estimate
      
      if(p < n) {
        chowlinestimates <- chowlin(X,y)
        cl.coef[[trial]] <- chowlinestimates$bhat
        cl.rho[trial] <- chowlinestimates$rho
        cl.series[[trial]] <- chowlinestimates$yhat
      }
      trueseries[[trial]] <- z.true # this is so we can compare z estimate to this truth
    }
    optimiser.list <- list("rho" = rho.final, "coefficients" = bhat.final, "series" = series.final, "lambda" = lambda.final,
                           "cl.rho" = cl.rho, "cl.coef" = cl.coef, "cl.series" = cl.series, "true.series" = trueseries,
                           "rfrho" = rfrho.final, "rfcoefficients" = rfbhat.final, "rfseries" = rfseries.final,
                           "rflambda" = rflambda.final)
    return(optimiser.list)
  }
 
## Estimation loop ## 
  # Functions used: AR.cov(), SParseTD.LARS.RF(), SParseTD.LARS(), agg.mat()
  SparseTD.estimates <- function(X,Y) {
    # 1.
    grid <- seq(0.01,0.99,by=0.01) # rho grid search
    ngrid <- length(grid)
    m <- nrow(X)
    bestrho <- rep(NA, ngrid)
    rfbestrho <- rep(NA, ngrid)
    
    # 2. 
    for(map in 1:ngrid) {
      V <- AR.cov(m, grid[map], 1) # define toeplitz form covariance for a rho in grid
      std.lars <- SparseTD.LARS(X,Y,V) # run algorithm using this covariance without re-fit
      std.lars.rf <- SparseTD.LARS.RF(X,Y,V) # # run algorithm using this covariance with re-fit
      bestrho[map] <- std.lars$bic.score # store the BIC for this rho - no re-fit
      rfbestrho[map] <- std.lars.rf$bic.score # store the BIC for this rho - re-fit
    }
    
    # 3. 
    rho.estimate <- grid[which.min(bestrho)] # rho minimising BIC  found
    vcov <- AR.cov(m,rho.estimate,1) # find the corresponding beta esstimate for this best rho
    lasso <- SparseTD.LARS(X,Y,vcov)
    beta.estimate <- lasso$beta.hat # store beta estimate
    lambda.estimate <- lasso$lambda.hat # store corresponding lambda 
    
    rfrho.estimate <- grid[which.min(rfbestrho)]
    rfvcov <- AR.cov(m,rfrho.estimate,1)
    rflasso <- SparseTD.LARS.RF(X,Y,rfvcov)
    rfbeta.estimate <- rflasso$beta.hat
    rflambda.estimate <- rflasso$lambda.hat
    
    # 4. 
    n <- length(Y) # find z estimate using this beta 
    C <- agg.mat(n,4) # CHANGE TO 3 FOR QUARTERLY TO MONTHLY 
    
    S <- C%*%tcrossprod(vcov,C)
    invS <- solve(S)
    series.estimate <- X%*%beta.estimate + tcrossprod(vcov,C)%*%(invS%*%(Y - C%*%(X%*%beta.estimate)))
    
    rfS <- C%*%tcrossprod(rfvcov,C)
    rfinvS <- solve(rfS)
    rfseries.estimate <- X%*%rfbeta.estimate + tcrossprod(rfvcov,C)%*%(rfinvS%*%(Y - C%*%(X%*%rfbeta.estimate)))
    
    output <- list("rho.estimate" = rho.estimate, "beta.estimate" = beta.estimate, 
                   "lambda.estimate" = lambda.estimate, "series.estimate" = series.estimate,
                   "rfrho.estimate" = rfrho.estimate, "rfbeta.estimate" = rfbeta.estimate,
                    "rfseries.estimate" = rfseries.estimate, "rflambda.estimate" = rflambda.estimate)
    return(output)
  }

## Algorithm for spTD WITH RE-FIT ## 
  # Functions used: agg.mat(), k.index(), tdbicFunc(), refit()
  SparseTD.LARS.RF <- function(X,Y,cov) {
    
    # 1.
    n <- length(Y)
    p <- ncol(X)
    C <- agg.mat(n,4)  # CHANGE TO 3 FOR QUARTERLY TO MONTHLY 
    Xl <- C%*%X
    S <- C%*%tcrossprod(cov,C)
    rotmat <- t(solve(chol(S)))
    A <- rotmat%*%Xl
    B <- rotmat%*%Y
    
    # 2. 
    lars <- lars(A,B,intercept=F)
    betamat <- lars$beta
    nbetamat = nrow(betamat)
    
    
    # 3. 
    npath <- k.index(betamat, n) # if p>n then it works out where floor(n/2) position is
    rfBIC <- rep(NA, npath)
    rfBIC[1] <- tdbicFunc(A,B,S,rep(0,p))
    rfbeta <- list()
    rfbeta[[1]] <- rep(0,p)
    for(path in 2:npath) {
      rfbeta[[path]] <- refit(A,B,betamat[path,],p) # first re-fit each path of LARS algorithm 
      rfBIC[path] <- tdbicFunc(A,B,S,rfbeta[[path]]) # second find the BIC score of this re-fitted beta estimate
    }
    rfbic.score <- min(rfBIC)
    idx.rfbic <- which.min(rfBIC)
    rfbeta.hat <- rfbeta[[idx.rfbic]]
    rflambda.hat <- lars$lambda[idx.rfbic]
    
    output <- list("bic.score" = rfbic.score, "beta.hat" = rfbeta.hat, "lambda.hat" = rflambda.hat)
    return(output)
  }
  
## Algorithm for spTD WITHOUT RE-FIT ##
  # Functions used: agg.mat(). k.index(), tdbicFunc()
  SparseTD.LARS <- function(X,Y,cov) {
    
    # 1.
    n <- length(Y)
    p <- ncol(X)
    C <- agg.mat(n,4) # CHANGE TO 3 FOR QUARTERLY TO MONTHLY 
    Xl <- C%*%X
    S <- C%*%tcrossprod(cov,C)
    rotmat <- t(solve(chol(S)))
    A <- rotmat%*%Xl
    B <- rotmat%*%Y
    
    # 2. 
    lars <- lars(A,B,intercept=F)
    betamat <- lars$beta
    nbetamat = nrow(betamat)
    
    
    # 3. 
    npath <- k.index(betamat, n)
    BIC <- rep(NA, npath)
    for(path in 1:npath) {
      BIC[path] <- tdbicFunc(A,B,S,betamat[path,])
    }
    bic.score <- min(BIC)
    idx.bic <- which.min(BIC)
    beta.hat <- betamat[idx.bic,]
    
    if(npath != idx.bic) {
      lambda.hat <- lars$lambda[idx.bic]
    }
    else {
      lambda.hat = 0
    }
    
    output <- list("bic.score" = bic.score, "lambda.hat" = lambda.hat, "beta.hat" = beta.hat)
    return(output)
  }
  
## Y.sim() - simulates data 
  Y.sim <- function(nq,beta,phi,sigma,phi.err,sigma.err) {
    p = length(phi)
    X <- matrix(NA, nrow = nq, ncol = p)
    for(isim in 1:p) {
      X[,isim] <- AR.sim(n=nq,phi=phi[isim],sigma=sigma[isim]) # CHANGE TO AR2.sim FOR RANDOM WALK!!
    }
    err <- AR.sim(n=nq,phi=phi.err,sigma=sigma.err)
    Y <- X%*%beta + err
    lt <- list("Y"=Y,"X"=X,"err"=err)
    return(lt)
  }

## AR.sim()/AR2.sim() - AR simulator
  AR.sim <- function(n, phi, sigma) {
    
    x <- rep(0,n)
    eps <- rnorm(n,0,sigma)
    x[1] <- rnorm(1,0,sigma/(1-phi^2))
    
    for(t in 2:n) {
      x[t] <- phi*x[t-1] + eps[t]
    }
    return(x)
  }
  AR2.sim <- function(n, phi, sigma) { # ONLY USE THIS FOR RANDOM WALK CASE
    
    x <- rep(0,n)
    eps <- rnorm(n,0,sigma)
    x[1] <- 0
    
    for(t in 2:n) {
      x[t] <- phi*x[t-1] + eps[t]
    }
    return(x)
  }
  
## chowlin() - does chow-lin temporal disaggregation 
  chowlin <- function(X,Y) {
    tsX <- ts(X, deltat = 1/4) # CHANGE TO 1/12 FOR MONTHLY
    tsY <- ts(Y) # CHANGE TO 1/4 FOR QUARTERLY
    cl <- td(tsY ~ 0 + tsX, to = "quarterly", truncated.rho = 0) # CHANGE TO "MONTHLY"
    bhat_cl <- cl$coefficients
    estrho <- cl$rho
    fittedval <- cl$values
    lst <- list("bhat" = bhat_cl, "yhat" = fittedval, "rho" = estrho)
    return(lst)
  }
  
## agg.mat() - aggregation matrix (M is low freq sample size and s is 4 for quarterly)
  agg.mat <- function(M,s) { 
    
    idmat <- diag(1,nrow=M)
    agvec <- rep(1,s)
    agmat <- t(kronecker(idmat,agvec))
    return(agmat)
    
  }
  
## AR.cov() - makes toeplitz matrix for AR(1) disturbances. Always set sigma to be 1!
  TheoACF <- function(phi,sigma,k) {
    return(sigma/(1-phi^2)*phi^abs(k))
  }
  AR.cov <- function(nq, phi, sigma) {
    return(toeplitz(sapply(0:(nq-1),TheoACF,phi=phi,sigma=sigma)))
  }
  
## k.index() - finds the floor(n/2) position 
  nonzeros <- function(vec) {
    sum(vec != 0)
  }
  k.index <- function(matrix, n) {
    count <- apply(matrix, 1, nonzeros)
    if(max(count) > n/2) {
      kindex <- min(which(count > n/2))
    }
    else {
      kindex <- nrow(matrix)
    }
  }
  
## tdbicFunc() - de-biased BIC formula 
  tdbicFunc <- function(X,Y,covariance,beta) {
    years <- length(Y)
    shat <- sum(beta != 0)
    u <- Y-X%*%beta
    log.lik <- -(years-shat)/2 - years/2*log(2*pi/(years-shat)*crossprod(u)) - log(det(covariance))/2
    BIC <- -2*(log.lik) + log(years)*shat
    return(BIC)
  }
  
## refit() - refits back into OLS 
  refit <- function(X,Y,bhat,p) {
    active <- which(bhat != 0 )
    X_new <- X[,active]
    lmfit <- lm(Y ~ 0 + X_new)
    bhat_lm <- lmfit$coefficients
    bhat_lm_full <- rep(0,p)
    bhat_lm_full[active] <- bhat_lm
    
    return(bhat_lm_full)
  }
  
  
  
  