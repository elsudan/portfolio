optim<-0

setwd("~/Desktop/Master_thesis/R/RCode")
options(scipen=999)

source("./KF_alt.R")
source("./likelihood_fct.R")

library(readxl)
library(matrixcalc)
library(optimx)


##### Get the data
{ 
  data <- read_xlsx("../Data/data_inf_var_clean_switz.xlsx",sheet=1)
  year <- data$year
  data$year <- NULL
  obs=3
  T<-dim(data)[1]
  n<-dim(data)[2]
  nb <- n/obs
  r <- 6*nb + obs
  A <- matrix(0,1,n)
  
  vector_0 <- matrix(0,nb,1)
  vector_1 <- matrix(1,nb,1)
  matrix_0 <- matrix(0,nb,nb)
  I <- diag(nb)
  
  R.1.2 <- 0 * diag(n) # s.t. R = (R.1.2)%*%t(R.1.2)
  R <- (R.1.2)%*%t(R.1.2)
  X=matrix(1,T,1)
} # parameters

Y<-data.matrix(data)

##### Initial conditions
{ 
  init_cond_values <- matrix(c(2,0,2.5,0,0,0,0,0,0, # Expected values
                             1,1,2,0.5,0.5,1,0.5,0.5,1), # Standard deviations
                             2,(2*obs)+obs,byrow=TRUE)
  xi.10 <- matrix(init_cond_values[1,obs+1],r,1)
  xi.10[1] <- init_cond_values[1,1]
  xi.10[1+nb+1] <- init_cond_values[1,2]
  xi.10[1+nb+1+nb+1] <- init_cond_values[1,3]
  
  P.10 <- init_cond_values[2,obs+1]^2*diag(r)
  P.10[1,1]<-init_cond_values[2,1]^2
  P.10[(1+nb+1),(1+nb+1)]<-init_cond_values[2,2]^2
  P.10[(1+nb+1+nb+1),(1+nb+1+nb+1)]<-init_cond_values[2,3]^2
  for(i in (1+nb+1+nb+2):(1+nb+1+nb+1+nb)){
    P.10[i,i]<-init_cond_values[2,2*obs]^2
  }
  for(i in (1+nb+1+nb+1+nb+1):(1+nb+1+nb+1+nb+nb)){
    P.10[i,i]<-init_cond_values[2,(2*obs)+1]^2
  }
  for(i in (1+nb+1+nb+1+nb+nb+1):(1+nb+1+nb+1+nb+nb+nb)){
    P.10[i,i]<-init_cond_values[2,(2*obs)+2]^2
  }
  for(i in (1+nb+1+nb+1+nb+nb+nb+1):(1+nb+1+nb+1+nb+nb+nb+nb)){
    P.10[i,i]<-init_cond_values[2,(2*obs)+3]^2
  }
}

##### Restrictions
{
  max_lambda <- 1.6
  max_sigma_rw <- sqrt(1/50)
  max_sigma_r <- sqrt(1/25)
  min_sigma_rw <- .09
  min_sigma_r <- .09
  max_sigma_tsw <- sqrt(1/50)
  max_sigma_ts <- sqrt(1/25)
  min_sigma_tsw <- .09
  min_sigma_ts <- .09
  max_sigma_piw <- sqrt(1/10)
  max_sigma_pi <- sqrt(1/4)
  min_sigma_piw <- 0
  min_sigma_pi <- 0
  max_phi <- .99
  
  restriction <- c(max_lambda,1,
                  max_sigma_rw,max_sigma_r,
                  max_sigma_tsw,max_sigma_ts,
                  max_sigma_piw,max_sigma_pi,
                  max_phi)
  
  sigmarw <- -log((max_sigma_rw-min_sigma_rw)/(0.1-min_sigma_rw)-1)
  sigmar <- -log((max_sigma_r-min_sigma_r)/(0.1-min_sigma_r)-1)
  sigmatsw <- -log((max_sigma_tsw-min_sigma_tsw)/(0.1-min_sigma_tsw)-1)
  sigmats <- -log((max_sigma_ts-min_sigma_ts)/(0.1-min_sigma_ts)-1)
  sigmapiw <- -log((max_sigma_piw-min_sigma_piw)/(sqrt(1/50)-min_sigma_piw)-1)
  sigmapi <- -log((max_sigma_pi-min_sigma_pi)/(sqrt(1/50)-min_sigma_pi)-1)
  
  phi <- -log((max_phi/0.2)-1)
  phi_ownlag <- -log((max_phi/0.001)-1)
}

##### Priors
param.0 <- c(matrix(-log(max_lambda-1),1,nb-1), # lambda
             phi_ownlag,phi,phi,phi,phi_ownlag,phi,phi,phi,phi_ownlag, # Phi
             sigmarw,matrix(sigmar,nb),
             sigmatsw,matrix(sigmats,nb),
             sigmapiw,matrix(sigmapi,nb),
             matrix(2,nb),matrix(2,nb),matrix(4,nb)) # Sigma

rm(sigmar,sigmarw,sigmats,sigmatsw,sigmapi,sigmapiw,phi,phi_ownlag)

log <- minus.log.likelihood.4.optimization(param.0,Y,X,xi.10,P.10,obs)

##### Optimization
if (optim==1){
  for(i in 1:10){
    for(METHOD in c("nlminb","Nelder-Mead")){
      res.optim <- optimx(param.0,
                          minus.log.likelihood.4.optimization,
                          Y=Y,X=X,xi.10=xi.10,P.10=P.10,obs=obs,
                          method = METHOD,
                          hessian=FALSE,
                          control = list(trace=TRUE,
                                         maxit=ifelse(METHOD=="Nelder-Mead",
                                                      2000,100),
                                         kkt=FALSE)
      )
      
      param.est <- c(as.matrix(res.optim)[1:length(param.0)])
      param.0 <- param.est
    }
  }
  save(param.est,file="./param_est.RData")
} else {
  load(file="./param_est.RData")
}

##### Get the estimated parameters
{
  {
    lambda_pi.est <- matrix(c((max_lambda)/(1+exp(-param.est[1:(nb-1)])),1),nb,1)
    row.names(lambda_pi.est)<-c("lambda CA","lambda FR","lambda DE","lambda IT","lambda JP",
                                "lambda CH","lambda GB","lambda US")
    # 1st line
    tH.est <- vector_1
    tH.est <- cbind(tH.est,I)
    tH.est <- cbind(tH.est,vector_0)
    tH.est <- cbind(tH.est,matrix_0)
    tH.est <- cbind(tH.est,lambda_pi.est)
    tH.est <- cbind(tH.est,I)
    tH.est <- cbind(tH.est,I)
    tH.est <- cbind(tH.est,matrix_0)
    tH.est <- cbind(tH.est,matrix_0)
    
    # 2nd line
    tH.est_2 <- vector_1
    tH.est_2  <- cbind(tH.est_2,I)
    tH.est_2  <- cbind(tH.est_2,vector_1)
    tH.est_2  <- cbind(tH.est_2,I)
    tH.est_2  <- cbind(tH.est_2,lambda_pi.est)
    tH.est_2  <- cbind(tH.est_2,I)
    tH.est_2  <- cbind(tH.est_2,matrix_0)
    tH.est_2  <- cbind(tH.est_2,I)
    tH.est_2  <- cbind(tH.est_2,matrix_0)
    tH.est <- rbind(tH.est,tH.est_2)
    
    # 3rd line
    tH.est_3 <- vector_0
    tH.est_3  <- cbind(tH.est_3,matrix_0)
    tH.est_3  <- cbind(tH.est_3,vector_0)
    tH.est_3  <- cbind(tH.est_3,matrix_0)
    tH.est_3  <- cbind(tH.est_3,lambda_pi.est)
    tH.est_3  <- cbind(tH.est_3,I)
    tH.est_3  <- cbind(tH.est_3,matrix_0)
    tH.est_3  <- cbind(tH.est_3,matrix_0)
    tH.est_3  <- cbind(tH.est_3,I)
    tH.est <- rbind(tH.est,tH.est_3)
    rm(tH.est_2,tH.est_3)
    
    H.est <- t(tH.est)
    
  } # H.est
  
  { phi_vector.est <- matrix(c(max_phi/(1+exp(-param.est[nb:(nb+(obs*2+obs)-1)]))))
    # phi <- c("phi_RR","phi_RRL","phi_Rpi","phi_RLR","phi_RLRL","phi_RLpi","phi_piR","phi_piRR","phi_pipi")
    
    Phi.est <- cbind(diag(x=phi_vector.est[1],nb),diag(x=phi_vector.est[2],nb),diag(x=phi_vector.est[3],nb))
    Phi.est <- rbind(Phi.est,cbind(diag(x=phi_vector.est[4],nb),diag(x=phi_vector.est[5],nb),diag(x=phi_vector.est[6],nb)))
    Phi.est <- rbind(Phi.est,cbind(diag(x=phi_vector.est[7],nb),diag(x=phi_vector.est[8],nb),diag(x=phi_vector.est[9],nb)))
    
    F.est <- diag(r)
    F.est[(r-n+1):(r),(r-n+1):(r)]<- Phi.est
    
  } # F.est
  
  { 
    sigma_rw <- min_sigma_rw + (max_sigma_rw-min_sigma_rw)/(1+exp(-param.est[nb+(obs*2+obs)]))
    # max_sigma_r/(1+exp(-param.est[nb+10]))
    sigma_rCA <- min_sigma_r + (max_sigma_r-min_sigma_r)/(1+exp(-param.est[nb+(obs*2+obs)+1]))
    sigma_rFR <- min_sigma_r + (max_sigma_r-min_sigma_r)/(1+exp(-param.est[nb+(obs*2+obs)+2]))
    sigma_rDE <- min_sigma_r + (max_sigma_r-min_sigma_r)/(1+exp(-param.est[nb+(obs*2+obs)+3]))
    sigma_rIT <- min_sigma_r + (max_sigma_r-min_sigma_r)/(1+exp(-param.est[nb+(obs*2+obs)+4]))
    sigma_rJP <- min_sigma_r + (max_sigma_r-min_sigma_r)/(1+exp(-param.est[nb+(obs*2+obs)+5]))
    sigma_rCH <- min_sigma_r + (max_sigma_r-min_sigma_r)/(1+exp(-param.est[nb+(obs*2+obs)+6]))
    sigma_rGB <- min_sigma_r + (max_sigma_r-min_sigma_r)/(1+exp(-param.est[nb+(obs*2+obs)+7]))
    sigma_rUS <- min_sigma_r + (max_sigma_r-min_sigma_r)/(1+exp(-param.est[nb+(obs*2+obs)+8]))
    
    sigma_tsw <- min_sigma_tsw+ (max_sigma_tsw-min_sigma_tsw)/(1+exp(-param.est[nb+(obs*2+obs)+9]))
    sigma_tsCA <- min_sigma_ts + (max_sigma_ts-min_sigma_ts)/(1+exp(-param.est[nb+(obs*2+obs)+10]))
    sigma_tsFR <- min_sigma_ts + (max_sigma_ts-min_sigma_ts)/(1+exp(-param.est[nb+(obs*2+obs)+11]))
    sigma_tsDE <- min_sigma_ts + (max_sigma_ts-min_sigma_ts)/(1+exp(-param.est[nb+(obs*2+obs)+12]))
    sigma_tsIT <- min_sigma_ts + (max_sigma_ts-min_sigma_ts)/(1+exp(-param.est[nb+(obs*2+obs)+13]))
    sigma_tsJP <- min_sigma_ts + (max_sigma_ts-min_sigma_ts)/(1+exp(-param.est[nb+(obs*2+obs)+14]))
    sigma_tsCH <- min_sigma_ts + (max_sigma_ts-min_sigma_ts)/(1+exp(-param.est[nb+(obs*2+obs)+15]))
    sigma_tsGB <- min_sigma_ts + (max_sigma_ts-min_sigma_ts)/(1+exp(-param.est[nb+(obs*2+obs)+16]))
    sigma_tsUS <- min_sigma_ts + (max_sigma_ts-min_sigma_ts)/(1+exp(-param.est[nb+(obs*2+obs)+17]))
    
    sigma_piw <- min_sigma_piw+ (max_sigma_piw-min_sigma_piw)/(1+exp(-param.est[nb+(obs*2+obs)+18]))
    sigma_piCA <- min_sigma_pi + (max_sigma_pi-min_sigma_pi)/(1+exp(-param.est[nb+(obs*2+obs)+19]))
    sigma_piFR <- min_sigma_pi + (max_sigma_pi-min_sigma_pi)/(1+exp(-param.est[nb+(obs*2+obs)+20]))
    sigma_piDE <- min_sigma_pi + (max_sigma_pi-min_sigma_pi)/(1+exp(-param.est[nb+(obs*2+obs)+21]))
    sigma_piIT <- min_sigma_pi + (max_sigma_pi-min_sigma_pi)/(1+exp(-param.est[nb+(obs*2+obs)+22]))
    sigma_piJP <- min_sigma_pi + (max_sigma_pi-min_sigma_pi)/(1+exp(-param.est[nb+(obs*2+obs)+23]))
    sigma_piCH <- min_sigma_pi + (max_sigma_pi-min_sigma_pi)/(1+exp(-param.est[nb+(obs*2+obs)+24]))
    sigma_piGB <- min_sigma_pi + (max_sigma_pi-min_sigma_pi)/(1+exp(-param.est[nb+(obs*2+obs)+25]))
    sigma_piUS <- min_sigma_pi + (max_sigma_pi-min_sigma_pi)/(1+exp(-param.est[nb+(obs*2+obs)+26]))
    sigma_pi_vec <- matrix(c(sigma_piCA, sigma_piFR, sigma_piDE,sigma_piIT,sigma_piJP,sigma_piCH,sigma_piGB,sigma_piUS))
    sigma_r_vec <- matrix(c(sigma_rCA, sigma_rFR, sigma_rDE,sigma_rIT,sigma_rJP,sigma_rCH,sigma_rGB,sigma_rUS))
    sigma_ts_vec <- matrix(c(sigma_tsCA, sigma_tsFR, sigma_tsDE,sigma_tsIT,sigma_tsJP,sigma_tsCH,sigma_tsGB,sigma_tsUS))
    sigma_vector.est<-c(sigma_rw,sigma_tsw,sigma_piw,sigma_r_vec,sigma_ts_vec,sigma_pi_vec)
    
    sigma_pi <- diag(nb)
    sigma_ts <- diag(nb)
    sigma_r <- diag(nb)
    sigma_Rtilde <- diag(nb)
    sigma_RLtilde <- diag(nb)
    sigma_pitilde <- diag(nb)
    
    sigma_RtildeCA <- param.est[nb+(obs*2+obs)+27]
    sigma_RtildeFR <- param.est[nb+(obs*2+obs)+28]
    sigma_RtildeDE <- param.est[nb+(obs*2+obs)+29]
    sigma_RtildeIT <- param.est[nb+(obs*2+obs)+30]
    sigma_RtildeJP <- param.est[nb+(obs*2+obs)+31]
    sigma_RtildeCH <- param.est[nb+(obs*2+obs)+32]
    sigma_RtildeGB <- param.est[nb+(obs*2+obs)+33]
    sigma_RtildeUS <- param.est[nb+(obs*2+obs)+34]
    sigma_Rtilde_vec <- matrix(c(sigma_RtildeCA, sigma_RtildeFR, sigma_RtildeDE,sigma_RtildeIT,sigma_RtildeJP,sigma_RtildeCH,sigma_RtildeGB,sigma_RtildeUS))
    
    sigma_RLtildeCA <- param.est[nb+(obs*2+obs)+35]
    sigma_RLtildeFR <- param.est[nb+(obs*2+obs)+36]
    sigma_RLtildeDE <- param.est[nb+(obs*2+obs)+37]
    sigma_RLtildeIT <- param.est[nb+(obs*2+obs)+38]
    sigma_RLtildeJP <- param.est[nb+(obs*2+obs)+39]
    sigma_RLtildeCH <- param.est[nb+(obs*2+obs)+40]
    sigma_RLtildeGB <- param.est[nb+(obs*2+obs)+41]
    sigma_RLtildeUS <- param.est[nb+(obs*2+obs)+42]
    sigma_RLtilde_vec <- matrix(c(sigma_RLtildeCA, sigma_RLtildeFR, sigma_RLtildeDE,sigma_RLtildeIT,sigma_RLtildeJP,sigma_RLtildeCH,sigma_RLtildeGB,sigma_RLtildeUS))
    
    sigma_pitildeCA <- param.est[nb+(obs*2+obs)+43]
    sigma_pitildeFR <- param.est[nb+(obs*2+obs)+44]
    sigma_pitildeDE <- param.est[nb+(obs*2+obs)+45]
    sigma_pitildeIT <- param.est[nb+(obs*2+obs)+46]
    sigma_pitildeJP <- param.est[nb+(obs*2+obs)+47]
    sigma_pitildeCH <- param.est[nb+(obs*2+obs)+48]
    sigma_pitildeGB <- param.est[nb+(obs*2+obs)+49]
    sigma_pitildeUS <- param.est[nb+(obs*2+obs)+50]
    sigma_pitilde_vec <- matrix(c(sigma_pitildeCA, sigma_pitildeFR, sigma_pitildeDE,sigma_pitildeIT,sigma_pitildeJP,sigma_pitildeCH,sigma_pitildeGB,sigma_pitildeUS))
    
    for (i in 1:nb){
      sigma_pi[i,i] <- sigma_pi_vec[i]
      sigma_r[i,i] <- sigma_r_vec[i]
      sigma_ts[i,i] <- sigma_ts_vec[i]
      sigma_Rtilde[i,i] <- sigma_Rtilde_vec[i]
      sigma_RLtilde[i,i] <- sigma_RLtilde_vec[i]
      sigma_pitilde[i,i] <- sigma_pitilde_vec[i]
    }
    
    Q.1.2.est <- diag(r)
    Q.1.2.est[1,1]<-sigma_rw
    Q.1.2.est[(2:(1+nb)),((2:(1+nb)))]<-sigma_r
    Q.1.2.est[1+nb+1,1+nb+1]<-sigma_tsw
    Q.1.2.est[(1+nb+2):(1+nb+1+nb),(1+nb+2):(1+nb+1+nb)]<-sigma_ts
    Q.1.2.est[1+nb+1+nb+1,1+nb+1+nb+1]<-sigma_piw
    Q.1.2.est[(1+nb+1+nb+2):(1+nb+1+nb+1+nb),(1+nb+1+nb+2):(1+nb+1+nb+1+nb)]<-sigma_pi
    Q.1.2.est[(1+nb+1+nb+1+nb+1):(1+nb+1+nb+1+nb+nb),
          (1+nb+1+nb+1+nb+1):(1+nb+1+nb+1+nb+nb)]<-sigma_Rtilde
    Q.1.2.est[(1+nb+1+nb+1+nb+nb+1):(1+nb+1+nb+1+nb+nb+nb),
          (1+nb+1+nb+1+nb+nb+1):(1+nb+1+nb+1+nb+nb+nb)]<-sigma_RLtilde
    Q.1.2.est[(1+nb+1+nb+1+nb+nb+nb+1):(1+nb+1+nb+1+nb+nb+nb+nb),
          (1+nb+1+nb+1+nb+nb+nb+1):(1+nb+1+nb+1+nb+nb+nb+nb)]<-sigma_pitilde
    
    Q.est <- Q.1.2.est%*%t(Q.1.2.est)
  } # Q.est
  
  {
    R.est <- R
    A.est <- A
  } # R, A
  
  { 
    Q_cyc <- Q.est[(1+nb+1+nb+1+nb+1):(1+nb+1+nb+1+nb+nb+nb+nb),(1+nb+1+nb+1+nb+1):(1+nb+1+nb+1+nb+nb+nb+nb)]
    P10_cyc <- matrix(solve(diag(n*n)-kronecker(Phi.est,Phi.est))%*%vec(Q_cyc),2*nb+nb,2*nb+nb)
    P.10[(1+nb+1+nb+1+nb+1):(1+nb+1+nb+1+nb+nb+nb+nb),
         (1+nb+1+nb+1+nb+1):(1+nb+1+nb+1+nb+nb+nb+nb)]<- P10_cyc
  } # P.10
  
}

##### Apply Kalman Filter with estimated parameters
{
  xi.est <- Kalman.filter(Y,X,F.est,Q.est,A.est,H.est,R.est,xi.10,P.10)[[1]]
  P.est <- Kalman.filter(Y,X,F.est,Q.est,A.est,H.est,R.est,xi.10,P.10)[[2]]
}

##### Smoothed estimations
{
  xi.est.smoothed <- Kalman.smoother(Y,X,F.est,Q.est,A.est,H.est,R.est,xi.10,P.10)[[1]]
  xi.est.smoothed <- rbind(xi.est.smoothed,matrix(NA,2020-tail(year,n=1),r))
  P.est.smoothed <- Kalman.smoother(Y,X,F.est,Q.est,A.est,H.est,R.est,xi.10,P.10)[[2]]
  P.est.smoothed <- rbind(P.est.smoothed,matrix(NA,2020-tail(year,n=1),r*r))
  
  colnames(xi.est.smoothed) <- c("rw","ri_Canada","ri_France","ri_Germany",
                               "ri_Italy","ri_Japan","ri_Switzerland","ri_UK","ri_USA",
                               "tsw","tsi_Canada","tsi_France","tsi_Germany",
                               "tsi_Italy","tsi_Japan","tsi_Switzerland","tsi_UK","tsi_USA",
                               "pw","pii_Canada","pii_France","pii_Germany",
                               "pii_Italy","pii_Japan","pii_Switzerland","pii_UK","pii_USA",
                               "Rtilde_Canada","Rtilde_France","Rtilde_Germany",
                               "Rtilde_Italy","Rtilde_Japan","Rtilde_Switzerland","Rtilde_UK","Rtilde_USA",
                               "RLtilde_Canada","RLtilde_France","RLtilde_Germany",
                               "RLtilde_Italy","RLtilde_Japan","RLtilde_Switzerland","RLtilde_UK","RLtilde_USA",
                               "piTilde_Canada","piTilde_France","piTilde_Germany",
                               "piTilde_Italy","piTilde_Switzerland","piTilde_Japan","piTilde_UK","piTilde_USA"
                               )
}

#### Create Del Negro variables
{
  # r_world + r_i
  {
    r_world <- data.frame(year=1870:2020)
    r_world$rw<-xi.est.smoothed[,1]
    r_i <- data.frame(year=1870:2020)
    countries=c("CA","FR","GE","IT","JP","CH","UK","US")
    for (j in 1:nb){
      r_i[,j+1]<-xi.est.smoothed[,1+j]+r_world[,2]
      colnames(r_i)[j+1]<-paste("r",countries[j],sep="_")
    }
    r_world<-merge(r_world,data.frame(year=1870:2020),by="year",all=TRUE)
    r_i<-merge(r_i,data.frame(year=1870:2020),by="year",all=TRUE)
  }
  
  # pi_world + pi_i
  {
    pi_world <- data.frame(year=1870:2020)
    pi_world$piw <- xi.est.smoothed[,1+nb+1+nb+1]
    pi_i <- data.frame(year=1870:2020)
    for (j in 1:nb){
      pi_i[,j+1]<-xi.est.smoothed[,1+nb+1+nb+1+j]+lambda_pi.est[j]*pi_world[,2]
      colnames(pi_i)[j+1]<-paste("pi",countries[j],sep="_")
    }
    pi_world<-merge(pi_world,data.frame(year=1870:2020),by="year",all=TRUE)
    pi_i<-merge(pi_i,data.frame(year=1870:2020),by="year",all=TRUE)
  }
  
  # ts_world + ts_i
  {
    ts_world <- data.frame(year=1870:2020)
    ts_world$tsw<-xi.est.smoothed[,1+nb+1]
    ts_i <- data.frame(year=1870:2020)
    countries=c("CA","FR","GE","IT","JP","CH","UK","US")
    for (j in 1:nb){
      ts_i[,j+1]<-xi.est.smoothed[,1+nb+1+j]+ts_world[,2]
      colnames(ts_i)[j+1]<-paste("ts",countries[j],sep="_")
    }
    ts_world<-merge(ts_world,data.frame(year=1870:2020),by="year",all=TRUE)
    ts_i<-merge(ts_i,data.frame(year=1870:2020),by="year",all=TRUE)
  }
}


### Calculations on r_world and r_i
{
  countries_DN <- c("CAN","FRA","DEU","ITA","JPN","CHE","GBR","USA")
  
  print(r_world[which(r_world$year==2008),2]-r_world[which(r_world$year==2016),2])
  for (i in 1:nb){
    print(cbind(countries_DN[i],
                xi.est.smoothed[which(r_world$year==1990),i+1]-xi.est.smoothed[which(r_world$year==2016),i+1],
                r_i[which(r_world$year==1990),i+1]-r_i[which(r_world$year==2016),i+1],
                (r_i[which(r_world$year==1990),i+1]-r_i[which(r_world$year==2016),i+1])-(xi.est.smoothed[which(r_world$year==1981),i+1]-xi.est.smoothed[which(r_world$year==2016),i+1])))
  }
}

#### Save the data
{
  save(init_cond_values,nb,obs,file="../RData/Switz/results_KF/init_conditions.RData")

  save(X,F.est,Q.est,A.est,H.est,R.est,xi.10,P.10,nb,obs,file="../RData/Switz/results_KF/matrices_est.RData")
  
  save(lambda_pi.est,phi_vector.est,sigma_vector.est,nb,obs,file="../RData/Switz/results_KF/estimated_param.RData")

  save(lambda_pi.est,xi.est.smoothed,P.est.smoothed,nb,obs,file="../RData/Switz/results_KF/estimated_var.RData")

  save(r_world,ts_world,pi_world,r_i,ts_i,pi_i,nb,obs,file="../RData/Switz/results_KF/estimated_DN_var.RData")
  
  save(restriction,file="../RData/Switz/results_KF/restrictions.RData")
}
