load("NimbleInputs8.RData")
library(nimble)
library(nimbleEcology)

###Nimble Code for full model (M39)
MA39<-nimbleCode({
  
###detection 'fixed' effects
  for (e in 1:2){ 
    mu_a[e]~dnorm(0, sd=1) 
    sig_a[e]~T(dnorm(0, sd = 1), 0, )
      for (i in 1:nspec){
        a[e,i]~dnorm(mu_a[e], sd=sig_a[e])
      }
    }
  
  sig_det~T(dnorm(0, sd = 1), 0, ) 

###occupancy/transition effects
  for (i in 1:30){
    mu_b[i]~dnorm(0, sd=0.75)
    sig_b[i]~T(dnorm(0, sd = 0.75), 0, )
    for (s in 1:nspec){
      b[i, s]~dnorm(mu_b[i], sd=sig_b[i])
    }
  }

###penalization for species loadings...
for (e in 1:2){  
  theta[1, e]~dgamma(10, 1)
  theta[2, e]~dgamma(10, 1)
  theta[3, e]~dgamma(10, 1)
  
  for (k in 1:3){
    tau[k, e]<-prod(theta[1:k, e])
    for (i in 1:nspec){
      lambda[i, k, e]~dt(0, tau=tau[k, e], 3)
    }  
  }  
}
  
  for (s in 1:nsites){
    for (i in 1:nspec){
      psi[i, s] <- iprobit(inprod(b[1:8, i], X1[s,1:8])+inprod(lambda[i, 1:3, 1], eta[s, 1:3, 1])) 
      phi[i, s] <- iprobit(inprod(b[9:19, i], X2[s,1:11])+inprod(lambda[i, 1:3, 2], eta[s, 1:3, 2]))
      gamma[i, s] <- iprobit(inprod(b[20:30, i], X2[s,1:11])+inprod(lambda[i, 1:3, 2], eta[s, 1:3, 2]))
      ###If we remove lambda and eta from the model, it reverts back to model 12.
      ###If we keep lambda and eta as is, but remove the last three predictors in X2 and re-index b--
      ###something like inprod(b[9:16, i], X2[s,1:8]) and inprod(b[17:24, i], X2[s,1:8])--
      ###we revert back to model 37.
      ###If we remove the the third index for lambda and the associated loop on L29, 
      ###we revert back to model 30. 

      y[s,1:2,1:5,i]~dDynOcc_ssm(init=psi[i, s], probPersist=phi[i, s],
                                  probColonize = gamma[i, s], p=p[s,1:2, 1:5,i],
                                  start=Starts2[s,1:2], end=Ends[s,1:2]) 

     ####one could also derive the point-wise log-likelihoods within the model as....
     ###loglike[s, i]<-dDynOcc_ssm(y[s,1:2,1:5,i],init=psi[i, s], probPersist=phi[i, s],
     #                             probColonize = gamma[i, s], p=p[s,1:2, 1:5,i],
     #                             start=Starts2[s,1:2], end=Ends[s,1:2], log=1) 
     ### ...note, this takes up a lot of storage and involves a redundant calculation.
     ###It is likely easier to derive a log-likelihood matrix after the fact.
    }
  }

  


  for (s in 1:nsites){
    for (e in 1:2){
      ###latent site factors...
      eta[s, 1, e]~dnorm(0, 1)
      eta[s, 2, e]~dnorm(0, 1)
      eta[s, 3, e]~dnorm(0, 1)
      for (j in 1:5){
        ###visit-specific random effect on detection
        eps[s, e, j]~dnorm(0, sd=sig_det)
        for (i in 1:nspec){
          logit(p[s, e, j, i])<-eps[s, e, j]+a[e, i] 
        }
      }
    }
  }
  
  
}) ###model end.

Constants<-list(Ends=ifelse(N_Surveys>5, 5, N_Surveys),
                nspec=185, nsites=320, Starts2=matrix(1, 320, 2), 
                X1=cbind(rep(1, 320),
                         as.numeric(scale(Site_Data2$MeanHistoricPPT)),
                         as.numeric(scale(Site_Data2$MeanHistoricTMin)),
                         as.numeric(scale(Site_Data2$HistoricTMinFPC2)),
                         as.numeric(scale(Site_Data2$Prop_Ag_Hist)),
                         as.numeric(scale(Site_Data2$PropUrbanHist)),
                         as.numeric(scale(Site_Data2$MeanHistoricTMin))^2,
                         as.numeric(scale(Site_Data2$HistoricTMinFPC2))^2),
                X2=cbind(rep(1, 320),
                         as.numeric(scale(Site_Data2$MeanModernPPT)),
                         as.numeric(scale(Site_Data2$MeanModernTMin)),
                         as.numeric(scale(Site_Data2$ModernTMinFPC2)),
                         as.numeric(scale(Site_Data2$PropAgMod)),
                         as.numeric(scale(Site_Data2$PropUrbanMod)),
                         as.numeric(scale(Site_Data2$MeanModernTMin))^2,
                         as.numeric(scale(Site_Data2$ModernTMinFPC2))^2,
                         as.numeric(scale(Site_Data2$MeanModernPPT-Site_Data2$MeanHistoricPPT)),
                         as.numeric(scale(Site_Data2$MeanModernTMin-Site_Data2$MeanHistoricTMin)),
                         as.numeric(scale(Site_Data2$ModernTMinFPC2-Site_Data2$HistoricTMinFPC2))))

Data<-list(y=y[, , 1:5, ])

##non-intuitive and probably not neccessary, but providing initial values for y seems to avoid unhappy warnings.
yIn<-array(NA, dim=dim(y))
yIn[is.na(y)]<-0


Inits <- list(y=yIn[,,1:5,], mu_b=rnorm(30, 0, 0.1), sig_b=rep(.25, 30), 
              a=matrix(rnorm(370, 0, .1), 2, 185), mu_a=rep(0, 2), sig_a=rep(.25, 2),
              b=matrix(rnorm(5550, 0, .1), 30, 185), 
eps=array(rnorm(3200, 0, .1), dim=c(320, 2, 5)),
              eta=array(0, dim=c(185, 3, 2)), lambda=array(0, dim=c(320, 3, 2)),
              theta=matrix(10, 3, 2))


Mod39 <- nimbleModel(code = MA39, name = 'M39', constants = Constants,
                         data=Data, inits=Inits, calculate=FALSE)

Mod39Conf <- configureMCMC(Mod39,
                               monitors = c("mu_a", "mu_b", "sig_a", "sig_b", "a", "b", "sig_det", "eps", "eta", "lambda"), UseConjugacy=FALSE) 

###Some hyper-sd's appear close to zero and are easier to sample on the log scale.
Mod39Conf$removeSamplers(c("sig_b[5]", "sig_b[13]", "sig_b[24]"))
Mod39Conf$addSampler(target ='sig_b[5]', type = 'RW', control=list(log=TRUE))
Mod39Conf$addSampler(target ='sig_b[13]', type='RW', control=list(log=TRUE))
Mod39Conf$addSampler(target = 'sig_b[24]', type= 'RW', control=list(log=TRUE))

Rmcmc<-buildMCMC(Mod39Conf)
compMCMC <- compileNimble(Rmcmc, Mod39)


###In practice, given HPC contstraints, we were able to run chains for about 70,000 iterations per job,
###with one chain per node. We stored the final samples, and used these as starting values for subsequent runs.
###In an ideal world, doing all of this at once (maybe in parallel), a useful posterior sample for this model
###might involve 6 chains with a burn-in of ~150000 iterations, niter=200000+, and thin=20 or so.

samps<-runMCMC(mcmc = compMCMC$Rmcmc,
               niter=70000, nburnin=50000, thin=10, 
               nchains=1) 





###Model 36--this is the most complete "autologistic model".


MA1_P1<-nimbleCode({
  for (e in 1:2){ ###detection 'fixed' effects
    mu_a[e]~dnorm(0, sd=1) 
    sig_a[e]~T(dnorm(0, sd = 1), 0, )
    for (i in 1:nspec){
      a[e,i]~dnorm(mu_a[e], sd=sig_a[e])
    }
  }
  
  sig_det~T(dnorm(0, sd = 1), 0, ) 
  
  for (i in 1:19){
    mu_b[i]~dnorm(0, sd=0.75)
    sig_b[i]~T(dnorm(0, sd = 0.75), 0, )
    for (s in 1:nspec){
      b[i, s]~dnorm(mu_b[i], sd=sig_b[i])
    }
  }
  
  mu_delta~dnorm(0, sd=0.75)
  sig_delta~T(dnorm(0, sd = 0.75), 0, )
  
  for (i in 1:nspec){
    delta[i]~dnorm(mu_delta, sd=sig_delta) 
    ###in text, we call this phi, but that will be used below. This is the additive effect of having occupied a site in the previous time.
  }
  
  for (e in 1:2){
    theta[1, e]~dgamma(10, 1)
    theta[2, e]~dgamma(10, 1)
    theta[3, e]~dgamma(10, 1)
    
    for (k in 1:3){
      tau[k, e]<-prod(theta[1:k, e])
      for (i in 1:nspec){
        lambda[i, k, e]~dt(0, tau=tau[k, e], 3)
      }
    }
  }
  
  for (s in 1:nsites){
    for (i in 1:nspec){
      psi[i, s] <- iprobit(inprod(b[1:8, i], X1[s,1:8])+inprod(lambda[i, 1:3, 1], eta[s, 1:3, 1]))
      phi[i, s] <- iprobit(inprod(b[9:19, i], X2[s,1:11])+delta[i]+inprod(lambda[i, 1:3, 2], eta[s, 1:3, 2]))
      gamma[i, s] <-iprobit(inprod(b[9:19, i], X2[s, 1:11])+inprod(lambda[i, 1:3, 2], eta[s, 1:3, 2]))
      y[s,1:2,1:5,i]~dDynOcc_ssm(init=psi[i, s], probPersist=phi[i, s],
                                 probColonize = gamma[i, s], p=p[s,1:2, 1:5,i],
                                 start=Starts2[s,1:2], end=Ends[s,1:2])

     ###Again, can drop lambda and eta and revert back to model 10, can hold lambda constant and drop back
     ###to model 27, can get rid of the last three predictors in X2 and b[17:19,] and drop back to model 34,
     ###can remove delta and drop back to model 31, etc.
    }
  }
  
  for (s in 1:nsites){
    for (e in 1:2){
      eta[s, 1, e]~dnorm(0, 1)
      eta[s, 2, e]~dnorm(0, 1)
      eta[s, 3, e]~dnorm(0, 1)
      for (j in 1:5){
        eps[s, e, j]~dnorm(0, sd=sig_det)
        for (i in 1:nspec){
          logit(p[s, e, j, i])<-eps[s, e, j]+a[e, i] 
        }
      }
    }
  }
  
 
  
})



###Data and constants are the same as for M39.
###here are some initial values that are guaranteed to start sampling.
Inits <- list(y=yIn[,,1:5,], mu_b=rep(0, 19), sig_b=rep(.25, 19), 
              a=matrix(0, 2, 185), mu_a=rep(0, 2), sig_a=rep(.25, 2),
              b=matrix(0, 19, 185), eps=array(0, dim=c(320, 2, 5)),
              mu_delta=0, sig_delta=.25, delta=rep(0, 185), sig_det=0.25,
              eta=array(0, dim=c(185, 3, 2)), lambda=array(0, dim=c(320, 3, 2)),
              theta=matrix(10, 3, 2))





###Model 33. This is the most complete "Time-varying" model without dependence.
MA33<-nimbleCode({
  for (e in 1:2){ ###detection 'fixed' effects
    mu_a[e]~dnorm(0, sd=1) 
    sig_a[e]~T(dnorm(0, sd = 1), 0, )
    for (i in 1:nspec){
      a[e,i]~dnorm(mu_a[e], sd=sig_a[e])
    }
  }
  
  sig_det~T(dnorm(0, sd = 1), 0, ) 
  
  for (i in 1:19){
    mu_b[i]~dnorm(0, sd=.75)
    sig_b[i]~T(dnorm(0, sd = .75), 0, )
    for (s in 1:nspec){
      b[i, s]~dnorm(mu_b[i], sd=sig_b[i])
    }
  }

for (e in 1:2){  
  theta[1,e]~dgamma(10, 1)
  theta[2,e]~dgamma(10, 1)
  theta[3,e]~dgamma(10, 1)
  
  for (k in 1:3){
    tau[k,e]<-prod(theta[1:k, e])
    for (i in 1:nspec){
      lambda[i, k, e]~dt(0, tau=tau[k,e], 3)
    }  
  }  
}  
  
  for (s in 1:nsites){
    for (i in 1:nspec){
      psi[i, s, 1] <- iprobit(inprod(b[1:8, i], X1[s,1:8])+inprod(lambda[i, 1:3, 1], eta[s, 1:3, 1]))
      psi[i, s, 2] <- iprobit(inprod(b[9:19, i], X2[s,1:11])+inprod(lambda[i, 1:3, 2], eta[s, 1:3, 2]))
      y[s,1:2,1:5,i]~dDynOcc_ssm(init=psi[i, s, 1], probPersist=psi[i, s, 2],
                                  probColonize = psi[i, s, 2], p=p[s,1:2, 1:5,i],
                                  start=Starts2[s,1:2], end=Ends[s,1:2]) 
    }
  }
  
  for (s in 1:nsites){
    for (e in 1:2){
      eta[s, 1, e]~dnorm(0, 1)
      eta[s, 2, e]~dnorm(0, 1)
      eta[s, 3, e]~dnorm(0, 1)
      for (j in 1:5){
        eps[s, e, j]~dnorm(0, sd=sig_det)
        for (i in 1:nspec){
          logit(p[s, e, j, i])<-eps[s, e, j]+a[e, i] 
        }
      }
    }
  }
  
 
})


Constants<-list(Ends=ifelse(N_Surveys>5, 5, N_Surveys),
                nspec=185, nsites=320, Starts2=matrix(1, 320, 2),
                X1=cbind(rep(1, 320),
                         as.numeric(scale(Site_Data2$MeanHistoricPPT)),
                         as.numeric(scale(Site_Data2$MeanHistoricTMin)),
                         as.numeric(scale(Site_Data2$HistoricTMinFPC2)),
                         as.numeric(scale(Site_Data2$Prop_Ag_Hist)),
                         as.numeric(scale(Site_Data2$PropUrbanHist)),
                         as.numeric(scale(Site_Data2$MeanHistoricTMin))^2,
                         as.numeric(scale(Site_Data2$HistoricTMinFPC2))^2),
                X2=cbind(rep(1, 320),
                         as.numeric(scale(Site_Data2$MeanModernPPT)),
                         as.numeric(scale(Site_Data2$MeanModernTMin)),
                         as.numeric(scale(Site_Data2$ModernTMinFPC2)),
                         as.numeric(scale(Site_Data2$PropAgMod)),
                         as.numeric(scale(Site_Data2$PropUrbanMod)),
                         as.numeric(scale(Site_Data2$MeanModernTMin))^2,
                         as.numeric(scale(Site_Data2$ModernTMinFPC2))^2,
                         as.numeric(scale(Site_Data2$MeanModernPPT-Site_Data2$MeanHistoricPPT)),
                         as.numeric(scale(Site_Data2$MeanModernTMin-Site_Data2$MeanHistoricTMin)),
                         as.numeric(scale(Site_Data2$ModernTMinFPC2-Site_Data2$HistoricTMinFPC2))))

###Note, when comparing era-specific coefficients for this model structure (e.g., did coefficients change or not)
###The predictors are on different scales for different time-periods. Need to rescale the coefficients after 
###fitting to make the comparison valid. 




Data<-list(y=y[,,1:5,])

yIn<-array(NA, dim=dim(y))
yIn[is.na(y)]<-0

Inits <- list(y=yIn[,,1:5,], mu_b=rnorm(19, 0, .1), sig_b=rep(.25, 19), 
              a=matrix(rnorm(370, 0, .1), 2, 185), mu_a=rnorm(2, 0, .1), sig_a=rep(.25, 2),
              b=matrix(rnorm(3515, 0, .1), 19, 185), sig_det=0.5, 
              eps=array(rnorm(3200, 0, .1), dim=c(320, 2, 5)),
              eta=array(0, dim=c(185, 3, 2)), lambda=array(0, dim=c(320, 3, 2)),
              theta=matrix(10, 3, 2))

#,mu_a=rep(0, 2), 
#sig_a=rep(.5, 2),

Mod34 <- nimbleModel(code = MA34, name = 'MA34', constants = Constants,
                         data=Data, inits=Inits, calculate=FALSE)

Mod34Conf <- configureMCMC(Mod34,
                               monitors = c("mu_a", "mu_b", "sig_a", "sig_b", "a", "b", "sig_det", "eps",
                                            "theta", "lambda", "eta"), UseConjugacy=FALSE) 

Mod34Conf$removeSamplers(c("sig_b[5]", "sig_b[13]"))
Mod34Conf$addSampler(target="sig_b[5]", type="RW", control=list(log=TRUE))
Mod34Conf$addSampler(target="sig_b[13]", type="RW", control=list(log=TRUE))

Rmcmc<-buildMCMC(Mod34Conf)
compMCMC <- compileNimble(Rmcmc, Mod34)
samps<-runMCMC(mcmc = compMCMC$Rmcmc,
               niter=70000, nburnin=60000, thin=5, 
               nchains=1)
