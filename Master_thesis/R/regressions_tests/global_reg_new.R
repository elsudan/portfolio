setwd("~/Desktop/Thesis/R/regressions_tests")
library(zoo)
library(sandwich)
library(data.table)


### Get & prepare data
{
  # Contains the estimates of r_world, r_i, p_world, etc..
  load(file="../RData/Switz/results_KF/estimated_DN_var.RData")
  # Loads lists which contain data frames of global variables.
  # The name of the list element is the name of the variable, all data frames are of length T=151
  # Global_var_list2 and global_var_list3 only differ by the weight to create the global variables
  load(file="../RData/Switz/exo/global/global_var.RData")
  # If using another list of variables, the name of the list should be changed to global_var_list
  global_var_list<-global_var_list2
  
  # Loads the US-only measures
  new_var <- load(file="../RData/Switz/exo/global/US_educattain_change.RData") # Load var and store name under "new_var"
  global_var_list[length(global_var_list)+1]<-list(get(new_var)) # Append new var to the list
  names(global_var_list)[length(global_var_list)]<-new_var # Change the name of the list element
  
  new_var<-load(file="../RData/Switz/exo/global/US_educattain.RData")
  global_var_list[length(global_var_list)+1]<-list(get(new_var))
  names(global_var_list)[length(global_var_list)]<-new_var
  
  new_var<-load(file="../RData/Switz/exo/global/US_gfd2gdp.RData")
  global_var_list[length(global_var_list)+1]<-list(get(new_var))
  names(global_var_list)[length(global_var_list)]<-new_var
  
  new_var<-load(file="../RData/Switz/cy_model/US_spread.RData")
  global_var_list[length(global_var_list)+1]<-list(get(new_var))
  names(global_var_list)[length(global_var_list)]<-new_var
  
  # Loads measures of catch-up growth :
  new_var<-load(file="../RData/Switz/exo/global/catchup.RData")
  global_var_list[length(global_var_list)+1]<-list(get(new_var))
  names(global_var_list)[length(global_var_list)]<-new_var
  
  new_var<-load(file="../RData/Switz/exo/global/catchup2.RData")
  global_var_list[length(global_var_list)+1]<-list(get(new_var))
  names(global_var_list)[length(global_var_list)]<-new_var
}

### Keep only the list of variables and the estimation of r_world
rm(list=ls()[! ls() %in% c("global_var_list","r_world")])

### Creates a table where we will store results of univariate regression
all_reg_coeff_pval<-data.table(variable=character(),coeff=numeric(),sign=character(),pval=numeric())

### OLS on all the elements of the list of variables
# Univariate (including constant) regressions with HAC SE, storing results in the table all_reg_coeff_pval
{
  for (i in 1:length(global_var_list)){
    var_name <- names(global_var_list)[[i]]
    variable <- global_var_list[[i]]
    ols <- lm(r_world$rw~variable[,2],na.action=na.exclude)
    coeff <- ols$coefficients[2]
    HAC<-vcovHAC(ols)[2,2]
    tstat<- coeff/sqrt(HAC)
    df<-ols$df.residual
    
    pval<- 2*pt(-abs(tstat),df=df)
    r2 <- summary(ols)$r.squared
    if (pval<0.01){
      sign<-"***"
    } else if (pval<0.05){
      sign<-"**"
    } else if (pval<0.1){
      sign<-"*"
    } else {
      sign<-" "
    }
    
    all_reg_coeff_pval<-rbind(all_reg_coeff_pval,data.table(variable=var_name,coeff=coeff,
                                                            sign=sign,pval=pval))
  }
}

### Remove unnecessary elements
rm(list=ls()[! ls() %in% c("global_var_list","r_world","all_reg_coeff_pval")])

### Choose variables by looking to the all_reg_coeff_pval table and according to p-values
{
  # Store the indices in increasing order
  var_indices <- c(3,5,6,7,9,12,13,15,18,20,22,23,24)
  # Create a new list keeping only the variables of choices
  global_kept_var<-global_var_list[var_indices]
  # Saving the list
  save(global_kept_var,file="../RData/Switz/global_kept/globals.RData")
}

### Doing the univariate OLS again and storing the results
{
  all_reg_coeff_pval<-data.table(variable=character(),coeff=numeric(),pval=numeric())
  for (i in var_indices){
    var_name <- names(global_var_list)[[i]]
    variable <- global_var_list[[i]]
    ols <- lm(r_world$rw~variable[,2],na.action=na.exclude)
    coeff <- ols$coefficients[2]
    HAC<-vcovHAC(ols)[2,2]
    tstat<- coeff/sqrt(HAC)
    df<-ols$df.residual
    
    pval<- 2*pt(-abs(tstat),df=df)
    r2 <- summary(ols)$r.squared
    
    all_reg_coeff_pval<-rbind(all_reg_coeff_pval,data.table(variable=var_name,coeff=coeff,pval=pval))
    if (exists("var_names")==FALSE){
      var_names<-matrix(var_name)
      coeff_list<-matrix(coeff)
      pval_list<-matrix(round(pval,3))
      r2_list<-matrix(r2)
      nbobs_list<-matrix(nobs(ols))
    } else {
      var_names<-cbind(var_names,matrix(var_name))
      coeff_list<-cbind(coeff_list,matrix(coeff))
      pval_list<-cbind(pval_list,matrix(round(pval,3)))
      r2_list<-cbind(r2_list,matrix(r2))
      nbobs_list<-cbind(nbobs_list,matrix(nobs(ols)))
    }
  }
}

### Doing multivariate regression on all variables selected in global_kept_var
{
  # Creates the formula
  formula <- "r_world$rw~"
  formula <- paste(formula,"global_var_list$",names(global_var_list)[var_indices[1]],"[,2]",sep="")
  for (i in var_indices[2:length(var_indices)]){
    formula <- paste(formula,"global_var_list$",sep="+")
    formula <- paste(formula,names(global_var_list)[i],"[,2]",sep="")
  }
  # Doing OLS and storing robust estimations + adjusted r2
  ols <- lm(formula=formula,na.action=na.exclude)
  allcoeff <- ols$coefficients[2:(length(var_indices)+1)]
  HAC<-diag(vcovHAC(ols)[2:(length(var_indices)+1),2:(length(var_indices)+1)])
  tstat<- allcoeff/sqrt(HAC)
  df<-ols$df.residual
  allpval<- 2*pt(-abs(tstat),df=df)
  allnbobs<-nobs(ols)
  allr2<-summary(ols)$adj.r.squared
  summary(ols)
  # Shows the coefficients and robust pvalues of the regression in a table
  robust_reg_result <- data.table(variable=matrix(var_names),pval=matrix(allpval),coeff=matrix(allcoeff))
  
  # Storing the results of the univariates regressions + the multivariate
  save(var_names,coeff_list,pval_list,r2_list,nbobs_list,allcoeff,allpval,allnbobs,allr2,
       file="../RData/Switz/tables/global/reg.RData")
}

rm(list=ls()[! ls() %in% c("global_var_list","r_world","all_reg_coeff_pval","robust_reg_result")])