setwd("~/Desktop/Thesis/R/regressions_tests")
options(scipen=999)

## Libraries
{
  library(readxl)
  library(matrixcalc)
  library(optimx)
  library(tidyr)
  library(sandwich)
  library(data.table)
  library(plm)
}

#### Retrieve data
{
  countries_DN <- c("CAN","FRA","DEU","ITA","JPN","CHE","GBR","USA")
  colors=c("orangered1","gold","royalblue3","seagreen3","plum4","hotpink","turquoise","peachpuff4")
  
  ### Get the previously estimated variables
  load(file="../RData/Switz/results_KF/estimated_DN_var.RData")  
  load(file="../RData/Switz/results_KF/estimated_var.RData")
  
  ### Loads all the variables stored in the folder "exo", stores names in the list "variables"
  variables <- list()
  for (item in list.files(path="../RData/Switz/exo")){
    if (item!="global"){
      var <-load(file=paste("../RData/Switz/exo/",item,sep=""))
      variables[[var]]<-load(file=paste("../RData/Switz/exo/",item,sep=""))
    }
  }
  
  ### Adding expected growth to the list of variables
  var <-load("../RData/Switz/expected_growth/exp_growth.RData")
  variables[[var[2]]]<-var[2]
}

### Creating long form data frame for panel regression
{
  df <- data.frame(year=matrix(1870:2020,nb*dim(r_world)[1],1))
  df$country <- rbind(matrix(countries_DN[1],dim(r_world)[1],1),
                      matrix(countries_DN[2],dim(r_world)[1],1),
                      matrix(countries_DN[3],dim(r_world)[1],1),
                      matrix(countries_DN[4],dim(r_world)[1],1),
                      matrix(countries_DN[5],dim(r_world)[1],1),
                      matrix(countries_DN[6],dim(r_world)[1],1),
                      matrix(countries_DN[7],dim(r_world)[1],1),
                      matrix(countries_DN[8],dim(r_world)[1],1))
  for (variable in c(variables)){
    new_df<-data.frame()
    for (country in countries_DN){
      new<-data.frame(year=1870:2020)
      new$new<-get(variable)[,which(colnames(get(variable))==
                             paste(country,substr(colnames(get(variable))[2],
                                                  4,nchar(colnames(get(variable))[2])),sep=""))]
      new$country<-country
      new_df<-rbind(new_df,new)
    }
    colnames(new_df)[2]<-variable
    df<-merge(df,new_df,by=c("year","country"))
  }
  r_i_long <- data.frame(year=matrix(1870:2020,nb*dim(r_world)[1],1))
  r_i_long$r_i <- rbind(matrix((xi.est.smoothed[,2]),dim(r_world)[1],1),
                        matrix((xi.est.smoothed[,3]),dim(r_world)[1],1),
                        matrix((xi.est.smoothed[,4]),dim(r_world)[1],1),
                        matrix((xi.est.smoothed[,5]),dim(r_world)[1],1),
                        matrix((xi.est.smoothed[,6]),dim(r_world)[1],1),
                        matrix((xi.est.smoothed[,7]),dim(r_world)[1],1),
                        matrix((xi.est.smoothed[,8]),dim(r_world)[1],1),
                        matrix((xi.est.smoothed[,9]),dim(r_world)[1],1))
  r_i_long$country <- rbind(matrix(countries_DN[1],dim(r_world)[1],1),
                        matrix(countries_DN[2],dim(r_world)[1],1),
                        matrix(countries_DN[3],dim(r_world)[1],1),
                        matrix(countries_DN[4],dim(r_world)[1],1),
                        matrix(countries_DN[5],dim(r_world)[1],1),
                        matrix(countries_DN[6],dim(r_world)[1],1),
                        matrix(countries_DN[7],dim(r_world)[1],1),
                        matrix(countries_DN[8],dim(r_world)[1],1))
  r_i_long<-data.frame(r_i_long)
  r_i_long$year<-c(r_i_long$year)
  r_i_long$r_i<-c(r_i_long$r_i)
  r_i_long$country<-c(r_i_long$country)
  
  df<-merge(df,r_i_long,by=c("year","country"))
}

### Keeping only the long data frame "df" and the list of variable names "variables"
rm(list=ls()[! ls() %in% c("variables","df")])

### Doing univariate panel regression with country fixed effect and storing the results in table
{
  all_reg_coeff_pval<-data.table(variable=character(),coeff=numeric(),pval=numeric())
  for (variable in variables){
    formula <- paste("r_i~",variable,sep="")
    plm<-plm(formula=formula,data=df,
                effect="individual",
                index="country")
    coeff <- summary(plm)$coefficients[1,1]
    r2 <- plm$r.squared[1]
    
    var_name <- variable
    HAC<-vcovNW(plm)[1]
    tstat<- coeff/sqrt(HAC)
    dF<-plm$df.residual
    pval<- 2*pt(-abs(tstat),df=dF)
    
    all_reg_coeff_pval<-rbind(all_reg_coeff_pval,data.table(variable=var_name,coeff=coeff,pval=pval))
    }
}

### Choosing the variables to use in the multivariate regression based on indices in "all_reg_coeff_pval" table
var_indices <- c(2,5,8,9,10,12,15,16,21)

### Panel regressions on chosen variables in "var_indices"
{
  ## Doing the univariate regressions and storing results
  all_reg_coeff_pval<-data.table(variable=character(),coeff=numeric(),pval=numeric())
  for (i in var_indices){
    formula <- paste("r_i~",names(variables)[i],sep="")
    plm<-plm(formula=formula,data=df,
             effect="individual",
             index="country")
    coeff <- summary(plm)$coefficients[1,1]
    r2 <- summary(plm)$r.squared[1]
    
    var_name <- names(variables)[i]
    HAC<-vcovNW(plm)[1]
    tstat<- coeff/sqrt(HAC)
    dF<-plm$df.residual
    pval<- 2*pt(-abs(tstat),df=dF)
    
    all_reg_coeff_pval<-rbind(all_reg_coeff_pval,data.table(variable=var_name,coeff=coeff,pval=pval))
    if (exists("var_names")==FALSE){
      var_names<-matrix(var_name)
      coeff_list<-matrix(coeff)
      pval_list<-matrix(round(pval,3))
      r2_list<-matrix(r2)
      nbobs_list<-matrix(nobs(plm))
    } else {
      var_names<-cbind(var_names,matrix(var_name))
      coeff_list<-cbind(coeff_list,matrix(coeff))
      pval_list<-cbind(pval_list,matrix(round(pval,3)))
      r2_list<-cbind(r2_list,matrix(r2))
      nbobs_list<-cbind(nbobs_list,matrix(nobs(plm)))
    }
  }
  
  ## Doing the multivariate regression on all variables in "var_indices"
  # Creating the formula
  formula <- "r_i~"
  for (i in var_indices){
    if (formula=="r_i~"){
      formula <- paste(formula,names(variables)[i],sep="")
      vec_name <- matrix(names(variables)[i])
    } else {
      formula <- paste(formula,names(variables)[i],sep="+")
      vec_name[length(vec_name)+1]<-matrix(names(variables)[i])
    }
  }
  plm<-plm(formula=formula,data=df,
           effect="individual",
           index="country")
  allcoeff <- summary(plm)$coefficients[,1]
  
  HAC<-diag(vcovNW(plm))
  tstat<- allcoeff/sqrt(HAC)
  dF<-plm$df.residual
  allpval<- 2*pt(-abs(tstat),df=dF)
  allnbobs<-nobs(plm)
  allr2<-summary(plm)$r.squared[2]
  
  # Tables with results
  robust_reg_result<-data.table(variable=matrix(vec_name),pval=matrix(allpval),coeff=matrix(allcoeff))
  
  # Savings results
  save(var_names,coeff_list,pval_list,r2_list,nbobs_list,allcoeff,allpval,allnbobs,allr2,
       file="../RData/Switz/tables/panel/reg_new.RData")
  
  rm(list=ls()[! ls() %in% c("variables","df","all_reg_coeff_pval")])
}
