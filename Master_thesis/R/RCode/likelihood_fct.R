minus.log.likelihood.4.optimization <- function(vector.param,
                                                Y,X,xi.10,P.10,obs){
  T <- dim(Y)[1]
  n <- dim(Y)[2]
  nb <- n/obs
  r <- dim(xi.10)[1]
  
  vector_0 <- matrix(0,nb,1)
  vector_1 <- matrix(1,nb,1)
  matrix_0 <- matrix(0,nb,nb)
  I <- diag(nb)
  
  # Build model:
  A<-matrix(0,1,n)
  
  {
    lambda_pi <- matrix(c((max_lambda)/(1+exp(-vector.param[1:(nb-1)])),1),nb,1)
    
    # 1st line
    tH <- vector_1
    tH <- cbind(tH,I)
    tH <- cbind(tH,vector_0)
    tH <- cbind(tH,matrix_0)
    tH <- cbind(tH,lambda_pi)
    tH <- cbind(tH,I)
    tH <- cbind(tH,I)
    tH <- cbind(tH,matrix_0)
    tH <- cbind(tH,matrix_0)
    
    # 2nd line
    tH_2 <- vector_1
    tH_2  <- cbind(tH_2,I)
    tH_2  <- cbind(tH_2,vector_1)
    tH_2  <- cbind(tH_2,I)
    tH_2  <- cbind(tH_2,lambda_pi)
    tH_2  <- cbind(tH_2,I)
    tH_2  <- cbind(tH_2,matrix_0)
    tH_2  <- cbind(tH_2,I)
    tH_2  <- cbind(tH_2,matrix_0)
    tH <- rbind(tH,tH_2)
    
    # 3rd line
    tH_3 <- vector_0
    tH_3  <- cbind(tH_3,matrix_0)
    tH_3  <- cbind(tH_3,vector_0)
    tH_3  <- cbind(tH_3,matrix_0)
    tH_3  <- cbind(tH_3,lambda_pi)
    tH_3  <- cbind(tH_3,I)
    tH_3  <- cbind(tH_3,matrix_0)
    tH_3  <- cbind(tH_3,matrix_0)
    tH_3  <- cbind(tH_3,I)
    tH <- rbind(tH,tH_3)
    rm(tH_2,tH_3)
    
    H <- t(tH)
    
  } # H
  
  { # Phi matrix
    phi_vector <- max_phi/(1+exp(-vector.param[nb:(nb+(obs*2+obs)-1)])) 
    # phi <- c("phi_RR","phi_RRL","phi_Rpi","phi_RLR","phi_RLRL","phi_RLpi","phi_piR","phi_piRR","phi_pipi")
    
    Phi <- cbind(diag(x=phi_vector[1],nb),diag(x=phi_vector[2],nb),diag(x=phi_vector[3],nb))
    Phi <- rbind(Phi,cbind(diag(x=phi_vector[4],nb),diag(x=phi_vector[5],nb),diag(x=phi_vector[6],nb)))
    Phi <- rbind(Phi,cbind(diag(x=phi_vector[7],nb),diag(x=phi_vector[8],nb),diag(x=phi_vector[9],nb)))
    
    F <- diag(r)
    F[(r-n+1):(r),(r-n+1):(r)]<- Phi
    
    
  } # F
  
  { 
    sigma_rw <- min_sigma_rw + (max_sigma_rw-min_sigma_rw)/(1+exp(-vector.param[nb+(obs*2+obs)]))
    # max_sigma_r/(1+exp(-vector.param[nb+10]))
    sigma_rCA <- min_sigma_r + (max_sigma_r-min_sigma_r)/(1+exp(-vector.param[nb+(obs*2+obs)+1]))
    sigma_rFR <- min_sigma_r + (max_sigma_r-min_sigma_r)/(1+exp(-vector.param[nb+(obs*2+obs)+2]))
    sigma_rDE <- min_sigma_r + (max_sigma_r-min_sigma_r)/(1+exp(-vector.param[nb+(obs*2+obs)+3]))
    sigma_rIT <- min_sigma_r + (max_sigma_r-min_sigma_r)/(1+exp(-vector.param[nb+(obs*2+obs)+4]))
    sigma_rJP <- min_sigma_r + (max_sigma_r-min_sigma_r)/(1+exp(-vector.param[nb+(obs*2+obs)+5]))
    sigma_rCH <- min_sigma_r + (max_sigma_r-min_sigma_r)/(1+exp(-vector.param[nb+(obs*2+obs)+6]))
    sigma_rGB <- min_sigma_r + (max_sigma_r-min_sigma_r)/(1+exp(-vector.param[nb+(obs*2+obs)+7]))
    sigma_rUS <- min_sigma_r + (max_sigma_r-min_sigma_r)/(1+exp(-vector.param[nb+(obs*2+obs)+8]))
    
    sigma_tsw <- min_sigma_tsw+ (max_sigma_tsw-min_sigma_tsw)/(1+exp(-vector.param[nb+(obs*2+obs)+9]))
    sigma_tsCA <- min_sigma_ts + (max_sigma_ts-min_sigma_ts)/(1+exp(-vector.param[nb+(obs*2+obs)+10]))
    sigma_tsFR <- min_sigma_ts + (max_sigma_ts-min_sigma_ts)/(1+exp(-vector.param[nb+(obs*2+obs)+11]))
    sigma_tsDE <- min_sigma_ts + (max_sigma_ts-min_sigma_ts)/(1+exp(-vector.param[nb+(obs*2+obs)+12]))
    sigma_tsIT <- min_sigma_ts + (max_sigma_ts-min_sigma_ts)/(1+exp(-vector.param[nb+(obs*2+obs)+13]))
    sigma_tsJP <- min_sigma_ts + (max_sigma_ts-min_sigma_ts)/(1+exp(-vector.param[nb+(obs*2+obs)+14]))
    sigma_tsCH <- min_sigma_ts + (max_sigma_ts-min_sigma_ts)/(1+exp(-vector.param[nb+(obs*2+obs)+15]))
    sigma_tsGB <- min_sigma_ts + (max_sigma_ts-min_sigma_ts)/(1+exp(-vector.param[nb+(obs*2+obs)+16]))
    sigma_tsUS <- min_sigma_ts + (max_sigma_ts-min_sigma_ts)/(1+exp(-vector.param[nb+(obs*2+obs)+17]))
    
    sigma_piw <- min_sigma_piw+ (max_sigma_piw-min_sigma_piw)/(1+exp(-vector.param[nb+(obs*2+obs)+18]))
    sigma_piCA <- min_sigma_pi + (max_sigma_pi-min_sigma_pi)/(1+exp(-vector.param[nb+(obs*2+obs)+19]))
    sigma_piFR <- min_sigma_pi + (max_sigma_pi-min_sigma_pi)/(1+exp(-vector.param[nb+(obs*2+obs)+20]))
    sigma_piDE <- min_sigma_pi + (max_sigma_pi-min_sigma_pi)/(1+exp(-vector.param[nb+(obs*2+obs)+21]))
    sigma_piIT <- min_sigma_pi + (max_sigma_pi-min_sigma_pi)/(1+exp(-vector.param[nb+(obs*2+obs)+22]))
    sigma_piJP <- min_sigma_pi + (max_sigma_pi-min_sigma_pi)/(1+exp(-vector.param[nb+(obs*2+obs)+23]))
    sigma_piCH <- min_sigma_pi + (max_sigma_pi-min_sigma_pi)/(1+exp(-vector.param[nb+(obs*2+obs)+24]))
    sigma_piGB <- min_sigma_pi + (max_sigma_pi-min_sigma_pi)/(1+exp(-vector.param[nb+(obs*2+obs)+25]))
    sigma_piUS <- min_sigma_pi + (max_sigma_pi-min_sigma_pi)/(1+exp(-vector.param[nb+(obs*2+obs)+26]))
    sigma_pi_vec <- matrix(c(sigma_piCA, sigma_piFR, sigma_piDE,sigma_piIT,sigma_piJP,sigma_piCH,sigma_piGB,sigma_piUS))
    sigma_r_vec <- matrix(c(sigma_rCA, sigma_rFR, sigma_rDE,sigma_rIT,sigma_rJP,sigma_rCH,sigma_rGB,sigma_rUS))
    sigma_ts_vec <- matrix(c(sigma_tsCA, sigma_tsFR, sigma_tsDE,sigma_tsIT,sigma_tsJP,sigma_tsCH,sigma_tsGB,sigma_tsUS))
    
    sigma_pi <- diag(nb)
    sigma_ts <- diag(nb)
    sigma_r <- diag(nb)
    sigma_Rtilde <- diag(nb)
    sigma_RLtilde <- diag(nb)
    sigma_pitilde <- diag(nb)
    
    sigma_RtildeCA <- vector.param[nb+(obs*2+obs)+27]
    sigma_RtildeFR <- vector.param[nb+(obs*2+obs)+28]
    sigma_RtildeDE <- vector.param[nb+(obs*2+obs)+29]
    sigma_RtildeIT <- vector.param[nb+(obs*2+obs)+30]
    sigma_RtildeJP <- vector.param[nb+(obs*2+obs)+31]
    sigma_RtildeCH <- vector.param[nb+(obs*2+obs)+32]
    sigma_RtildeGB <- vector.param[nb+(obs*2+obs)+33]
    sigma_RtildeUS <- vector.param[nb+(obs*2+obs)+34]
    sigma_Rtilde_vec <- matrix(c(sigma_RtildeCA, sigma_RtildeFR, sigma_RtildeDE,sigma_RtildeIT,sigma_RtildeJP,sigma_RtildeCH,sigma_RtildeGB,sigma_RtildeUS))
    
    sigma_RLtildeCA <- vector.param[nb+(obs*2+obs)+35]
    sigma_RLtildeFR <- vector.param[nb+(obs*2+obs)+36]
    sigma_RLtildeDE <- vector.param[nb+(obs*2+obs)+37]
    sigma_RLtildeIT <- vector.param[nb+(obs*2+obs)+38]
    sigma_RLtildeJP <- vector.param[nb+(obs*2+obs)+39]
    sigma_RLtildeCH <- vector.param[nb+(obs*2+obs)+40]
    sigma_RLtildeGB <- vector.param[nb+(obs*2+obs)+41]
    sigma_RLtildeUS <- vector.param[nb+(obs*2+obs)+42]
    sigma_RLtilde_vec <- matrix(c(sigma_RLtildeCA, sigma_RLtildeFR, sigma_RLtildeDE,sigma_RLtildeIT,sigma_RLtildeJP,sigma_RLtildeCH,sigma_RLtildeGB,sigma_RLtildeUS))
    
    sigma_pitildeCA <- vector.param[nb+(obs*2+obs)+43]
    sigma_pitildeFR <- vector.param[nb+(obs*2+obs)+44]
    sigma_pitildeDE <- vector.param[nb+(obs*2+obs)+45]
    sigma_pitildeIT <- vector.param[nb+(obs*2+obs)+46]
    sigma_pitildeJP <- vector.param[nb+(obs*2+obs)+47]
    sigma_pitildeCH <- vector.param[nb+(obs*2+obs)+48]
    sigma_pitildeGB <- vector.param[nb+(obs*2+obs)+49]
    sigma_pitildeUS <- vector.param[nb+(obs*2+obs)+50]
    sigma_pitilde_vec <- matrix(c(sigma_pitildeCA, sigma_pitildeFR, sigma_pitildeDE,sigma_pitildeIT,sigma_pitildeJP,sigma_pitildeCH,sigma_pitildeGB,sigma_pitildeUS))
    
    for (i in 1:nb){
      sigma_pi[i,i] <- sigma_pi_vec[i]
      sigma_r[i,i] <- sigma_r_vec[i]
      sigma_ts[i,i] <- sigma_ts_vec[i]
      sigma_Rtilde[i,i] <- sigma_Rtilde_vec[i]
      sigma_RLtilde[i,i] <- sigma_RLtilde_vec[i]
      sigma_pitilde[i,i] <- sigma_pitilde_vec[i]
    }
    
    Q.1.2 <- diag(r)
    Q.1.2[1,1]<-sigma_rw
    Q.1.2[(2:(1+nb)),((2:(1+nb)))]<-sigma_r
    Q.1.2[1+nb+1,1+nb+1]<-sigma_tsw
    Q.1.2[(1+nb+2):(1+nb+1+nb),(1+nb+2):(1+nb+1+nb)]<-sigma_ts
    Q.1.2[1+nb+1+nb+1,1+nb+1+nb+1]<-sigma_piw
    Q.1.2[(1+nb+1+nb+2):(1+nb+1+nb+1+nb),(1+nb+1+nb+2):(1+nb+1+nb+1+nb)]<-sigma_pi
    Q.1.2[(1+nb+1+nb+1+nb+1):(1+nb+1+nb+1+nb+nb),
          (1+nb+1+nb+1+nb+1):(1+nb+1+nb+1+nb+nb)]<-sigma_Rtilde
    Q.1.2[(1+nb+1+nb+1+nb+nb+1):(1+nb+1+nb+1+nb+nb+nb),
          (1+nb+1+nb+1+nb+nb+1):(1+nb+1+nb+1+nb+nb+nb)]<-sigma_RLtilde
    Q.1.2[(1+nb+1+nb+1+nb+nb+nb+1):(1+nb+1+nb+1+nb+nb+nb+nb),
          (1+nb+1+nb+1+nb+nb+nb+1):(1+nb+1+nb+1+nb+nb+nb+nb)]<-sigma_pitilde
    
    Q <- Q.1.2%*%t(Q.1.2)
  } # Q
  
  {
    R.1.2 <- 0.0001 * diag(n)
    R <- R.1.2%*%t(R.1.2)
    X <- matrix(1,T,1)
  } # R, X
  
  { 
    Q_cyc <- Q[(1+nb+1+nb+1+nb+1):(1+nb+1+nb+1+nb+nb+nb+nb),(1+nb+1+nb+1+nb+1):(1+nb+1+nb+1+nb+nb+nb+nb)]
    P10_cyc <- matrix(solve(diag(n*n)-kronecker(Phi,Phi))%*%vec(Q_cyc),2*nb+nb,2*nb+nb)
    P.10[(1+nb+1+nb+1+nb+1):(1+nb+1+nb+1+nb+nb+nb+nb),
         (1+nb+1+nb+1+nb+1):(1+nb+1+nb+1+nb+nb+nb+nb)]<- P10_cyc
  } # P.10
  
  # Call KF:
  
  logL <- Kalman.filter(Y,X,F,Q,A,H,R,xi.10,P.10)[[5]]
  
  return(-logL)
}