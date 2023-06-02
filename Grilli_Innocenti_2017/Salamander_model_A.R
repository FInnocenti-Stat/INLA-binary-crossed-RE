###Example Salamander data, model A  Karim and Zeger####
rm(list=ls(all=T))
## load the INLA library
library (INLA)
library(foreign)
library(lme4)

## load the salamander data
salam <- read.dta("salam.dta")   # Dataset from Stata Manual
salam
head(salam)
attach(salam)
wmf=wsm*wsf
wmf                             # Interaction tra wsm e wsf

data<-as.data.frame(cbind(y=y,wsm=wsm, wsf=wsf, wmf=wmf,fe=female, m=male, group=group ))
#######################################  LAPLACE  #############################

modLp<-glmer(y~ 1+ wsm+wsf+wmf + (1|fe) + (1|m) , data, family=binomial(link=logit),nAGQ = 1L)
summary(modLp)

beta<-fixef(modLp)                   
sdbeta<-sqrt(diag(vcov(modLp)))                         
var.lap<-getME(modLp, "theta")^2       
round(sqrt(var.lap),2)

#######################################   INLA    ##############################

# INLA default prior
formula= y ~ wsm+wsf+wmf+f(fe, model="iid", hyper=list(theta=list(prior="loggamma",param=c(1,0.0005))))+ f(m, model="iid", hyper=list(theta=list(prior="loggamma",param=c(1,0.0005))))
#INLA GAMMA(0.001,0.001)
formula= y ~ wsm+wsf+wmf+f(fe, model="iid", hyper=list(theta=list(prior="loggamma",param=c(0.001,0.001))))+ f(m, model="iid", hyper=list(theta=list(prior="loggamma",param=c(0.001,0.001))))
#INLA GAMMA(0.5,0.0164)
formula= y ~ wsm+wsf+wmf+f(fe, model="iid", hyper=list(theta=list(prior="loggamma",param=c(0.5,0.0164))))+ f(m, model="iid", hyper=list(theta=list(prior="loggamma",param=c(0.5,0.0164))))
#INLA GAMMA(0.5,0.003737)
formula= y ~ wsm+wsf+wmf+f(fe, model="iid", hyper=list(theta=list(prior="loggamma",param=c(0.5,0.003737))))+ f(m, model="iid", hyper=list(theta=list(prior="loggamma",param=c(0.5,0.003737))))

result = inla(formula,family="binomial", data = data)

estbeta<-result$summary.fixed[,1]                    
sdebeta<-result$summary.fixed[,2]                    

prec.marg<-(result$marginals.hyperpar$'Precision for fe') 
marg.variance <-inla.tmarginal(function(x) 1/x, prec.marg)  
ifelse(sum(marg.variance[,2])==0,var<-0,var<-inla.emarginal(function(x) x, marg.variance))
ifelse(sum(marg.variance[,2])==0,mm<-0,mm<-inla.emarginal(function(x) x^2, marg.variance))
sd_var<-sqrt(mm-(var)^2)  

prec.margB<-(result$marginals.hyperpar$'Precision for m') 
marg.varianceB <-inla.tmarginal(function(x) 1/x, prec.margB)  
ifelse(sum(marg.varianceB[,2])==0,varB<-0,varB<-inla.emarginal(function(x) x, marg.varianceB))
ifelse(sum(marg.varianceB[,2])==0,mmB<-0,mmB<-inla.emarginal(function(x) x^2, marg.varianceB))
sd_varB<-sqrt(mmB-(varB)^2)  

round(sqrt(var),2);round(sqrt(varB),2)

