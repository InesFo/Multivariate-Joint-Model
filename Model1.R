# Model 1-------------------------------------------------------------------
# YS_ij = DY * (b0 + b1*time_ij + U0_i + ZY_ij) +
#         DS * (b2 + U0_i + ZS_ij)



# libraries ---------------------------------------------------------------
library(nlme)
library(dplyr)
library(beepr)

start.time <- Sys.time()




# Definition of parameters ------------------------------------------------

{
  # simulation values v0	tauY	tauS
  simValues<-rbind(c(1,0.1,0.1),
                   c(1,0.25,0.1))
  colnames(simValues)<-c("niu0","tauY","tauS")
  
  # number subjects
  n = 200
  
  # number of observations per subject
  # all subjects are observed same number of times
  m = 5
  mi <- rep(m,n)
  N <- sum(mi) 
  
  # betas b0, b1, b2
  beta <- c(15, 1.2, 2)
  
  # number of simulations per each k set
  nsim = 1000
}




# for each k = 1, ..., simvalues
for (k in 1:dim(simValues)[1]) {
  
  # niu0 = sd(U0)  
  niu0 = simValues[k,1]
  
  # tauY = sd(Zy)
  tauY <- simValues[k,2]
  
  # tauS = sd(Zs)
  tauS <- simValues[k,3]
  
  
  # lists to save results
  lmelist = list()
  datalist = list()
  
  # extensive list of seeds
  seeds<-1:(nsim*100)
  
  
  #Create data frame with the parameters
  value <- c(beta[1],beta[2],beta[3],niu0,tauY,tauS,"NA","NA")
  parameter <- c("beta.DY","beta.DY:time","beta.DS","sd(U0)","sd(Y)","sd(S)","seed","nsubj")
  df.par<-as.data.frame(cbind(parameter,value))
  
  
  # Function to create data set with simulated data ----------------------------------------------
  
  createVardf <- function(npart,mobsvec){
    
    # dataframe: Y11, Y12, ... , Ynm, ..., S1, ..., Sn 
    N.tot = sum(mi)
    
    # participant id
    id <- c(rep(1:npart, mobsvec),1:npart)
    # Y measurement time
    time <- c(sequence(mobsvec,from=0),rep(0,npart))
    # indicator variable (DY = 1 for rows with Y and DY = 0 for rows with S)
    DY <- c(rep(1,N.tot),rep(0,npart))
    # indicator variable (DS = 1 for rows with S and DS = 0 for rows with Y)
    DS <- c(rep(0,N.tot),rep(1,npart))
    # factor variable (variable = Y for rows with Y and variable = S for rows with S)
    variable <- c(rep("Y",N.tot),rep("S",npart))
    data.frame(id,time,DY,DS,variable)
  }
  
  # starting index for successful models
  i.suc = 0
  # starting index of total simulations
  i.tot = 0
  
  
  # fit several models and save results on lmelist[[]] ----------------------
  
  
  # while the number of successfully saved models is less than the number of simulations: iterate
  while (i.suc < nsim ) {
    
    # increment index of total simulations
    i.tot=i.tot+1
    
    # set seed
    set.seed(seeds[i.tot])
    
    # Get Z and U distributions -----------------------------------------------
    Zy <- rnorm(N, 0, tauY)
    Zs <- rnorm(n, 0, tauS)
    U <- rnorm(n, 0, niu0)
    
    # Generate data -----------------------------------------------------------
    dfvar<- createVardf(n,mi)
    
    # Add columns to dfvar with parameters to compute Y_ij
    dfvar$b0 <- c(rep(beta[1],N),rep(0,n))
    dfvar$b1 <- c(rep(beta[2],N),rep(0,n))
    dfvar$b2 <- c(rep(0,N),rep(beta[3],n))
    dfvar$Z <- c(Zy,Zs)
    dfvar$U0 <- c(rep(U, mi),U)
    # multivariate data frame in the long format
    mvdf.l <- dfvar %>% mutate(value = DY*(b0 + b1*time + U0 + Z) + DS*(b2 + U0 + Z)) 
    
    
    # Select Y measurements before failure time F (note: F = exp(S)) ------------------------------
    
    # vector with failure times (length = n)
    Fvalues<-exp(mvdf.l[mvdf.l$variable=="S","value"])
    # add vector with failure times to data frame mvdf.l
    mvdf.l$Fvalues<-c(rep(Fvalues, times=mi),Fvalues)
    # select rows in which Y observed time<Fvalues
    mvdf.l<-mvdf.l[mvdf.l$time<mvdf.l$Fvalues,]
    
    
    
    # Fit lme model -----------------------------------------------------------
    m.YS1<-tryCatch(lme(value ~ 0 +  DY + DY:time + DS, data =mvdf.l,
                        random = ~ 1 | id, 
                        weights = varIdent(form = ~1 |variable)),
                    warning = function(w) {"warning"},
                    error = function(e) {"error"})
    
    # if successfully fitted model
    if (class(m.YS1)=="lme") {
      
      # increment index for successful models
      i.suc = i.suc + 1
      
      # attribute seed used to generate the dataset to the simulated dataset
      attr(mvdf.l,"seed")<-seeds[i.tot]
      # save the simulated dataset in list datalist[]
      datalist[[(i.suc)]] <- mvdf.l
      
      # attribute seed used to generate the dataset to the fitted model
      attr(m.YS1,"seed")<-seeds[i.tot]
      # attribute number of subjects (nsubj) to the fitted model
      attr(m.YS1, "nsubj") <- length(unique(mvdf.l[mvdf.l$variable=="Y","id"]))
      # save results of fitted model in list lmelist[]
      lmelist[[i.suc]] <- m.YS1
      
      # print in the console
      print(paste0("fitted model ",k,"-", i.suc))  }
    
  }
  
  # When i.suc = nsim 
  # Save results of each model in df.par dataframe -----------------------------
  for (i in 1:nsim) {
    
    #Create dataframes with estimates of parameters
    
    # b0,b1,b2
    betas.res<-fixed.effects(lmelist[[i]])
    betas.df<-data.frame(parameter=c("beta.DY","beta.DS","beta.DY:time"),value=betas.res)
    
    #niu0
    niu0.res<-VarCorr(lmelist[[i]])[1,2]
    niu0.df<-data.frame(parameter="sd(U0)", value = niu0.res)
    
    #tauY, tauS
    tauY.res<-as.double(VarCorr(lmelist[[i]])[2,2])
    tauS.res<-as.double(coef(lmelist[[i]]$modelStruct$varStruct, unconstrained=FALSE, allCoef=TRUE)[2])*tauY.res
    tau.df<-data.frame(parameter=c("sd(Y)","sd(S)"), value=c(tauY.res,tauS.res))
    
    #seed
    seed.res <- data.frame(parameter="seed",value=attr(lmelist[[i]],"seed"))
    
    #nsubj
    nsubj.res <- data.frame(parameter="nsubj",value=attr(lmelist[[i]],"nsubj"))
    
    
    # join rows with estimated parameters
    df.res<-rbind(betas.df,niu0.df,tau.df,seed.res,nsubj.res)
    
    # reorder rows such that the order is the same as in the data frame with the parameters (df.par)
    df.res <- df.res %>% slice(match(df.par$parameter, parameter))
    # round all values to 4 decimal places
    df.res$value[1:(nrow(df.res)-2)] <- round(as.double(df.res$value[1:(nrow(df.res)-2)]),6)
    #change to decimal
    df.res$value<- format(as.double(df.res$value), scientific = FALSE) 
    
    # new column with values from lmelist[[i]]
    new <- df.res$value
    # append new column to df.par
    df.par[ , ncol(df.par) + 1] <- new
    # rename column
    colnames(df.par)[ncol(df.par)] <- paste0("Trial", i)
    # rename rows
    rownames(df.par)<-1:nrow(df.par)
    
    # df.par has the form (e.g. for nsim = 3):
    #
    # parameter      value    Trial1    Trial2    Trial3    
    #        beta.DY    15 15.460356 14.556277 15.190695 
    #   beta.DY:time   1.2  1.212939  1.227745  1.212804  
    #        beta.DS     2  2.465556  1.632999  2.193105  
    #         sd(U0)     1  1.277459  1.226924  0.719341  
    #          sd(Y)  0.25  0.223334  0.287319  0.213496  
    #          sd(S)   0.1  0.118267  0.131146  0.105031 
    #           seed    NA  1.000000  2.000000  3.000000  
    #          nsubj    NA 200.00000 200.00000 200.00000 
    # 
  }
  
  # Check mean, SD, min, max, CV of estimated parameters in the nsim simulations ---------------------------------------------------------
  
  # copy the df.par dataframe
  df.par.final<- df.par
  # change columns values to double
  df.par.final[3:ncol(df.par.final)] <- lapply(df.par.final[3:ncol(df.par.final)], as.double)
  
  # compute mean of estimated parameters in the nsim simulations
  df.par.final<-transform(df.par.final, Mean=apply(df.par.final[,3:ncol(df.par.final)],1, mean, na.rm = TRUE))
  # compute sd of estimated parameters in the nsim simulations
  df.par.final<-transform(df.par.final, SD  =apply(df.par.final[,3:(ncol(df.par.final)-1)],1, sd, na.rm = TRUE))
  # compute min of estimated parameters in the nsim simulations
  df.par.final<-transform(df.par.final, Min =apply(df.par.final[,3:(ncol(df.par.final)-2)],1, min, na.rm = TRUE))
  # compute max of estimated parameters in the nsim simulations
  df.par.final<-transform(df.par.final, Max =apply(df.par.final[,3:(ncol(df.par.final)-3)],1, max, na.rm = TRUE))
  # compute coefficient of variation of estimated parameters in the nsim simulations
  df.par.final$CV <- round(df.par.final$SD/df.par.final$Mean,4)
  
  # show all columns except the last with rge statistics regarding the estimated parameters in the nsim simulations
  df.par.final[,c(1,2,(ncol(df.par.final)-4):ncol(df.par.final))]
  
  # the result is like:
  #
  #     parameter value       Mean         SD       Min       Max     CV
  #       beta.DY    15 14.9965422 0.42379212 14.499698 15.754384 0.0283
  #  beta.DY:time   1.2  1.2040700 0.02217429  1.160382  1.227745 0.0184
  #       beta.DS     2  2.0121927 0.40440815  1.476410  2.689868 0.2010
  #        sd(U0)     1  0.9049757 0.26049116  0.391741  1.277459 0.2878
  #         sd(Y)  0.25  0.2415466 0.02706554  0.198809  0.287319 0.1121
  #         sd(S)   0.1  0.1068460 0.04301845  0.027400  0.194384 0.4026
  #          seed    NA  5.5000000 3.02765035  1.000000 10.000000 0.5505
  #         nsubj    NA 200.000000 0.00000000 200.00000 200.00000 0.0000
  
  
  beep(2)
  
  # Save results ------------------------------------------------------------
  write.csv(df.par.final, file = paste0("./df.par.final-niu_",niu0,"tauY_",tauY,"tauS_",tauS,".csv"))
  saveRDS(df.par.final,   file = paste0("./df.par.final-niu_",niu0,"tauY_",tauY,"tauS_",tauS,".rds"))
  saveRDS(lmelist, file = paste0("./lmelist-niu_",niu0,"tauY_",tauY,"tauS_",tauS,".rds")) 
  
  # print in the console
  print(paste("saved fitted model", k))
  
}

# set end time
end.time <- Sys.time()
# set total duration
time.taken <- round(end.time - start.time,2)
# print total duration
time.taken


