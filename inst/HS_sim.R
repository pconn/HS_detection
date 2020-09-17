#test model using simulated data

HS_data = read.csv('file:///C:/Users/paul.conn/git/HS_detection/HarborSeal_DisturbanceInAleutians_20190702_SKH.csv')

Poly_df = data.frame(ID = c("OB06","OB05","OB04","OB07","OB08","OA02","OA00","OB13","OB12","PA04","PA03","PA07","PA11","QA02","QA01","PA14","PA12","SA19","SA20","SA03","SA00","UA06","UA03","UA00","TA22","TA21","XD00","XD01","XD04","XD05","XD06","XD07","XD09","XD10"))
Poly_df$strata = c(rep("E",17),rep("C",9),rep("W",8))
n_poly = nrow(Poly_df)
Poly_df$mean_count = rep(0,n_poly)
Poly_df$priority = rep("NA",n_poly)  
for(ipoly in 1:n_poly){
  Cur_dat = HS_data[as.character(HS_data$polyid)==Poly_df$ID[ipoly],]
  Poly_df$mean_count[ipoly]=round(mean(Cur_dat$num_seals_photos+Cur_dat$num_disturbed))
  Poly_df$priority[ipoly]=Cur_dat[1,"priority"]  #1- undefined, 2=high, 3=low, 4=medium
} 
Poly_df$mean_count[18:20]=mean(Poly_df$mean_count[21:26])  #these weren't surveyed since 2015
Poly_df$mean_count[34]=1  #change mean count = 0 to mean count = 1

library(mefa)
Design_df = rep(Poly_df,each=4)
Design_df$visit = rep(c(1,1,2,2),n_poly)
Design_df$rep = rep(c(1,2),n_poly)
Weather = sample(c(0,1),n_poly*2,prob=c(2,1),replace=TRUE)
Design_df$weather = rep(Weather,each=2)
Thermal = sample(c(0,1),n_poly*2,replace=TRUE)
Design_df$thermal = as.vector(rbind(Thermal,1-Thermal))
Design_df$trial = rep(c(0,1),2*n_poly)
Design_df$alt = rep(800,n_poly*4)
Which_thermal = which(Design_df$thermal==1)
Design_df$alt[Which_thermal]=runif(length(Which_thermal),800,1500)
Design_df$alt = Design_df$alt/1000

exp.count <- function(Design_df){
  beta0 = -0.1
  beta1 = .05 
  alpha0 = -2.9
  alpha1 = 5
  alpha2 = -2   
  weather_eff = -.2
  sigma_v = 0.1
  sigma_t = 0.02
  X= as.matrix(cbind(log(Design_df$mean_count+0.01),(Design_df$visit-1),(Design_df$visit-1)*Design_df$alt,Design_df$thermal,Design_df$thermal*Design_df$alt,Design_df$thermal*Design_df$alt^2,Design_df$weather*Design_df$thermal))
  Par = matrix(c(1,beta0,beta1,alpha0,alpha1,alpha2,weather_eff),ncol=1)
  Lin_pred = X%*%Par + rnorm(nrow(X),0,sigma_t)
  Epsilon_v = rnorm(nrow(X)/2,0,sigma_v)
  Lin_pred = Lin_pred + rep(Epsilon_v,each=2) 
  Exp_count = exp(Lin_pred)
  Out=list(Exp_count = Exp_count,X=X)
  return(Out)
}
Out = exp.count(Design_df)

crap = Out$Exp_count
#lm(log(crap)~0+X[,1]+X[,2]+X[,3]+X[,4]+X[,5]+X[,6]+X[,7])

#jags model
library(R2jags)
jags.HS <- function(){
  for(i in 1:n_count){
    Count[i]~dpois(Exp_count[i])
    #Count[i]~dnorm(Exp_count[i],1.0)
    Exp_count[i]<-exp(X[i,1]*intercept + X[i,2]*beta0 + X[i,3]*beta1 + X[i,4]*alpha0+X[i,5]*alpha1+X[i,6]*alpha2+X[i,7]*weather_eff+Eps_visit[i]+Eps_count[i])
    Eps_count[i]~dnorm(0,tau_t)
  }
  
  Eps_visit = X_vis %*% Eps_visit_2
  for(i in 1:n_vis){
    Eps_visit_2[i]~dnorm(0,tau_v)
  }
  
  #priors
  intercept~dnorm(0,.01)
  beta0~dnorm(0,.01)
  beta1~dnorm(0,.01)
  alpha0~dnorm(0,.01)
  alpha1~dnorm(0,.01)
  alpha2~dnorm(0,.01)
  weather_eff~dnorm(0,.01)
  tau_t ~ dgamma(1.0,0.1)
  tau_v ~ dgamma(1.0,0.1)
}

n_count = 4*n_poly
X_vis = matrix(0,n_count,n_count/2)
counter=1
for(i in 1:n_count){
  X_vis[i,counter]=1
  if(i%%2==0)counter = counter+1
}
data <- list(X=Out$X,Count=as.integer(as.vector((Out$Exp_count))),n_count=as.integer(n_poly*4),X_vis=X_vis,n_vis=2*n_poly)

init.values <- function(){
  list(intercept=runif(1,.8,1.2),beta0=.1*rnorm(1),beta1=.1*rnorm(1),alpha0=.1*rnorm(1),alpha1=.1*rnorm(1),
       alpha2=.1*rnorm(1),weather_eff=.1*rnorm(1),tau_v=runif(1,50,150),tau_t=runif(1,50,150))
}

# MCMC settings
set.seed(1482312)  
ada=1000  
iter=13000      # 500000
nchains=3  
nburnin=3000      #  100000
nthin=1  #100  
inits <- init.values()  

## Parameters to be monitored (= to estimate)
params =  c("intercept","beta0","beta1","alpha0","alpha1","alpha2","weather_eff","tau_v","tau_t")
jags_save = c(params) #,"Exp_count","Eps_visit")
jags_fit = jags(data=data,
                inits=init.values,
                params,
                n.iter=iter,
                model.file=jags.HS,
                DIC=TRUE,
                parameters.to.save=jags_save,
                n.chains=nchains,
                n.burnin=nburnin,
                n.thin=nthin,
                working.directory=getwd()) 

library(coda)
library(mcmcplots)
samps <- as.mcmc(jags_fit)

pdf("MCMC.pdf")
  plot(samps)
dev.off()

mcmcplot(samps)