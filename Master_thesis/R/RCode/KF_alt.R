Kalman.filter <- function(Y,X,F,Q,A,H,R,xi.10,P.10){
  T<-dim(Y)[1]
  n<-dim(Y)[2]
  r<-dim(xi.10)[1]
  
  matrix.xi.tt <- matrix(NA,T,r)
  matrix.xi.tp1t <- matrix(NA,T,r)
  matrix.P.tt <- matrix(NA,T,r^2)
  matrix.P.tp1t <- matrix(NA,T,r^2)
  
  log <- 0
  
  for (t in 1:T) {
    indic.available.prev <- which(!is.na(Y[t,]))
    
    Y_na <- matrix(Y[t,])[indic.available.prev,]
    H_na <- H[,indic.available.prev]
    R_na <- R[indic.available.prev,indic.available.prev]
    A_na <- t(matrix(A[,indic.available.prev]))
    
    if (t==1) {
      P.ttm1 <- P.10
      xi.ttm1 <- xi.10
    }
    
    else {
      P.ttm1 <- matrix(matrix.P.tp1t[t-1,],r,r)
      xi.ttm1 <- matrix(matrix.xi.tp1t[t-1,],r,1)
    }
    
    gain_t <- P.ttm1 %*% H %*% 
      solve(t(H) %*% P.ttm1 %*% H + R)
    
    mu_t <- t(A) %*% matrix(X[t,]) + t(H) %*% xi.ttm1
    sigma_t <- (t(H) %*% P.ttm1 %*% H + R)[indic.available.prev,indic.available.prev]
    
    prev_err <- (matrix(Y[t,]) - mu_t)[indic.available.prev,]
    
    xi.tt <- xi.ttm1 + gain_t[,indic.available.prev] %*% prev_err
    P.tt <- P.ttm1 - gain_t %*% t(H) %*% P.ttm1
    
    xi.tp1t <- F %*% xi.tt 
    P.tp1t <- F %*% P.tt %*% t(F) + Q
    
    matrix.xi.tt[t,] <- c(xi.tt)
    matrix.xi.tp1t[t,] <- c(xi.tp1t)
    matrix.P.tt[t,] <- c(P.tt)
    matrix.P.tp1t[t,] <- c(P.tp1t)
    
    log <- log -length(prev_err)/2*log(2*pi) - .5*log(det(sigma_t)) - 
      .5*t(prev_err) %*% solve(sigma_t) %*% prev_err
  }
  
  return(list(matrix.xi.tt,matrix.P.tt,matrix.xi.tp1t,matrix.P.tp1t,c(log)))
}

Kalman.smoother <- function(Y,X,F,Q,A,H,R,xi.10,P.10) {
  
  results_KF <- Kalman.filter(Y,X,F,Q,A,H,R,xi.10,P.10)
  
  T<-dim(Y)[1]
  n<-dim(Y)[2]
  r<-dim(xi.10)[1]
  
  matrix.xi.tt <- results_KF[[1]]
  matrix.P.tt <- results_KF[[2]]
  matrix.xi.tp1t <- results_KF[[3]]
  matrix.P.tp1t <- results_KF[[4]]
  
  matrix.xi.tT <- matrix(NA,T,r)
  matrix.P.tT <- matrix(NA,T,r^2)
  
  matrix.xi.tT[T,]<-matrix.xi.tt[T,]
  matrix.P.tT[T,]<-matrix.P.tt[T,]
  
  for (t in (T-1):1) {
    P.tt <- matrix(matrix.P.tt[t,],r,r)
    P.tp1t <- matrix(matrix.P.tp1t[t,],r,r)
    P.tp1T <- matrix(matrix.P.tT[t+1,],r,r)
    
    xi.tt <- matrix(matrix.xi.tt[t,],r,1)
    xi.tp1t <- matrix(matrix.xi.tp1t[t,],r,1)
    xi.tp1T <- matrix(matrix.xi.tT[t+1,],r,1)
    
    J_t <- P.tt %*% t(F) %*% solve(P.tp1t)
    xi.tT <- xi.tt + J_t %*% (xi.tp1T - xi.tp1t)
    P.tT <- P.tt + J_t %*% (P.tp1T - P.tp1t) %*% t(J_t)
    
    matrix.xi.tT[t,]<- c(xi.tT)
    matrix.P.tT[t,]<- c(P.tT)
  }
  
  return(list(matrix.xi.tT,matrix.P.tT))
}