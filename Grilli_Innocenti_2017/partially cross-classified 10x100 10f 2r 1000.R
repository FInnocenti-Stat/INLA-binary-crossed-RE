##################################################################################################
#######################  Scenarios Partial Cross-Classification  ###################################
#######################      10 Feeders and 2 Receivers        ###################################    
##################################################################################################
rm(list=ls(all=T))
library(INLA)
library(lme4)
set.seed(1234)

sim<-500     

# set fixed effects 
beta<-c(0.100000,0.100000,0.400000,0.400000,0.400000) 
betasize<-length(beta) 
totalsize<-rowcount<-1
fixname<-c("x0","x1","x2","x3","x4")          
dataset<-vector("list",3+length(fixname))
names(dataset)<-c("l2idA","l2idB","y",fixname)

# INLA  Ga(0.001,0.001)
finaloutput<-matrix(0,totalsize,2*betasize)    
sdebeta<-matrix(0,sim,betasize)                
estbeta<-matrix(0,sim,(betasize))              
betaaverage<-matrix("list",betasize)           
names(betaaverage)<-c("b0","b1","b2","b3","b4") 
sdeaverage<-betaaverage                        
sigma<-sigmaU<-var<-prec.marg<-marg.variance<-sd_var<-sd<-mm<-NULL
sigmaB<-sigmaUB<-varB<-prec.margB<-marg.varianceB<-sd_varB<-sdB<-mmB<-NULL

# INLA  Ga(0.5,0.003737)
finaloutput2<-matrix(0,totalsize,2*betasize)    
sdebeta2<-matrix(0,sim,betasize)                
estbeta2<-matrix(0,sim,(betasize))              
betaaverage2<-matrix("list",betasize)           
names(betaaverage2)<-c("b0","b1","b2","b3","b4") 
sdeaverage2<-betaaverage2                        
sigma2<-sigmaU2<-var2<-prec.marg2<-marg.variance2<-sd_var2<-sd2<-mm2<-NULL
sigmaB2<-sigmaUB2<-varB2<-prec.margB2<-marg.varianceB2<-sd_varB2<-sdB2<-mmB2<-conto<-conto2<-NULL

# MLLA
finaloutput.lap<-matrix(0,totalsize,2*betasize)    
sdebeta.lap<-matrix(0,sim,betasize)               
estbeta.lap<-matrix(0,sim,(betasize))             
var.lap<-matrix(0,nrow=sim, ncol=2)
betaaverage.lap<-matrix("list",betasize)           
names(betaaverage.lap)<-c("b0.l","b1.l","b2.l","b3.l","b4.l") 
sdeaverage.lap<-betaaverage.lap                       
var_lapA<-var_lapB<-sigma.lap<-sigmaB.lap<-NULL


freq.receivers<-matrix(0,sim,20)
freq.feeders<-matrix(0,sim,10)

############ start simulations:
# Generate a 3-level hierarchical structure as follows
# 100 observations (i.e. level 1 units) within each level 2 unit 
# 10  level 2 units (Feeders) within each level 3 unit
# 10  level 3 units (Receivers)
# The partially cross-classified structure is composed of
n1=100 # number of observations per feeder
nF=10  # total number of Feeders
nR=10  # total number of Receivers

l1id<-rep(c(1:n1))              # id of level 1 units
l2idF<-rep(c(1:nF),each=n1)     # id of Feeders
l2idR<-rep(c(1:nR),each=nF*n1)  # id of Receivers
cbind(l1id,l2idF,l2idR)         # structure of the cross-classification matrix

length=n1*nF*nR                 # total number of observations in the dataset                  
x<-matrix(1,length,betasize)               

#generate covariates 
x[,2]<-rnorm(length,0,1)
x[,3]<-rbinom(length,1,0.5)
predF<-rbinom(nF*nR,1,0.5)
x[,4]<-rep(predF,length.out=length,each=n1)   
predR<-rbinom(nR,1,0.5)
ifelse(sum(predR==0),predR<-c(1,0,1,0,0,0,0,1,1,1),predR<-predR)
ifelse(sum(predR==1),predR<-c(1,0,1,0,0,0,0,1,1,1),predR<-predR)
x[,5]<-rep(predR,length.out=length,each=nF)   
x0<-x[,1]                                 
x1<-x[,2]
x2<-x[,3]
x3<-x[,4]
x4<-x[,5]

##             Five specifications for the variance components:
sigmau<-matrix(c(0.10),1,1)      # st. dev. random effect Factor A  (0.1, 0.1, 0.5, 0.5, 1.0)
sigmaub<-matrix(c(0.10),1,1)     # st. dev. random effect Factor B  (0.1, 0.5, 0.1, 0.5, 1.0)


for(iter in 1:sim){    
  
 cat(" dataset number", iter ,"\n")
  
  idfact2=sample(c(1:100), 10, replace=FALSE)          # draw 10 Feeders
  
X<-rbind(x[((idfact2[1]-1)*100+1):(idfact2[1]*100),],  # draw the observations corresponding to the 10 sampled Feeders
   x[((idfact2[2]-1)*100+1):(idfact2[2]*100),],
   x[((idfact2[3]-1)*100+1):(idfact2[3]*100),],
   x[((idfact2[4]-1)*100+1):(idfact2[4]*100),],
   x[((idfact2[5]-1)*100+1):(idfact2[5]*100),],
   x[((idfact2[6]-1)*100+1):(idfact2[6]*100),],
   x[((idfact2[7]-1)*100+1):(idfact2[7]*100),],
   x[((idfact2[8]-1)*100+1):(idfact2[8]*100),],
   x[((idfact2[9]-1)*100+1):(idfact2[9]*100),],
   x[((idfact2[10]-1)*100+1):(idfact2[10]*100),]) 
 
   X0<-X[,1]                                 
   X1<-X[,2]
   X2<-X[,3]
   X3<-X[,4]
   X4<-X[,5]
 
  feeders<-rep(idfact2,each=n1)                          # assign each observation to the corresponding feeder
  freq.feeders[iter,]<-idfact2                           # store the sampled feeders for each replicate
  receivers<-feeders
  
  idreceiver1=sample(c(1:10), 10, replace=FALSE)         # draw the first receiver 
  idreceiver2=sample(c(1:10), 10, replace=FALSE)         # draw the second receiver
  
  for(i in 1:10){    # to avoid that, given a feeder, the two selected receivers coincide 
    ifelse(idreceiver1[i]==idreceiver2[i],idreceiver2[i]<-sample(idreceiver2[-i], 1, replace=FALSE),idreceiver1[i]<-idreceiver1[i])
  }
  
  # to ensure that all the 10 receivers are drawn at least one time (no empty column in the cross-classification matrix)
  receivers[1:n1]<-c(rep(idreceiver1[1],each=(n1/2)),rep(idreceiver2[1],each=(n1/2)))
  receivers[(n1+1):(n1*2)]<-c(rep(idreceiver1[2],each=(n1/2)),rep(idreceiver2[2],each=(n1/2)))
  receivers[((n1*2)+1):(n1*3)]<-c(rep(idreceiver1[3],each=(n1/2)),rep(idreceiver2[3],each=(n1/2)))
  receivers[((n1*3)+1):(n1*4)]<-c(rep(idreceiver1[4],each=(n1/2)),rep(idreceiver2[4],each=(n1/2)))
  receivers[((n1*4)+1):(n1*5)]<-c(rep(idreceiver1[5],each=(n1/2)),rep(idreceiver2[5],each=(n1/2)))
  receivers[((n1*5)+1):(n1*6)]<-c(rep(idreceiver1[6],each=(n1/2)),rep(idreceiver2[6],each=(n1/2)))
  receivers[((n1*6)+1):(n1*7)]<-c(rep(idreceiver1[7],each=(n1/2)),rep(idreceiver2[7],each=(n1/2)))
  receivers[((n1*7)+1):(n1*8)]<-c(rep(idreceiver1[8],each=(n1/2)),rep(idreceiver2[8],each=(n1/2)))
  receivers[((n1*8)+1):(n1*9)]<-c(rep(idreceiver1[9],each=(n1/2)),rep(idreceiver2[9],each=(n1/2)))
  receivers[((n1*9)+1):(n1*10)]<-c(rep(idreceiver1[10],each=(n1/2)),rep(idreceiver2[10],each=(n1/2)))
  
  # store the drawn receivers for each replicate
  freq.receivers[iter,]<-recei<-c(idreceiver1[1],idreceiver2[1],idreceiver1[2],idreceiver2[2],
                                  idreceiver1[3],idreceiver2[3],idreceiver1[4],idreceiver2[4],
                                  idreceiver1[5],idreceiver2[5],idreceiver1[6],idreceiver2[6],
                                  idreceiver1[7],idreceiver2[7],idreceiver1[8],idreceiver2[8],
                                  idreceiver1[9],idreceiver2[9],idreceiver1[10],idreceiver2[10])

  
  betaaverage[1:betasize]<-rep(0,betasize)
  sdeaverage<-sdeaverage.lap<-betaaverage.lap<-betaaverage2<-betaaverage                    
  
  
  u<-rnorm(100,0,sigmau)    #   Feeders random effects 
  v<-rnorm(10,0,sigmaub)    #   Receivers random effects 
  
  U<-V<-NULL
  
  for(s in 1:10){
    U[s]<-u[idfact2[s]]   # assign each feeder to the corresponding random effect
  }
  
  for(s in 1:20){         # assign each receiver to the corresponding random effect
    for(t in 1:20){
      ifelse(recei[s]==recei[t],V[s]<-v[recei[s]],V[s]<-v[recei[s]]) 
    }
    
  }
  
  
  rand<-rep(U,each=n1)     #   Feeders random effects
  part<-rep(V,each=(n1/2)) #   Receivers random effects
  randpart<-rand+part

  fixpart<-X%*%beta          # X refers to the partially cross-classified dataset derived from the 3-level hierarchical structure
  
  binomprob<-exp(fixpart+randpart)/(1+exp(fixpart+randpart)) 
  y<-rbinom(length(fixpart),1,binomprob)             
  
  feeders<-as.factor(feeders)                      
  receivers<-as.factor(receivers) 
  dataset<-as.data.frame(cbind(X0,X1,X2,X3,X4,feeders,receivers,y))       
  
  # INLA Ga(0.001,0.001)
  cat(" INLA Ga(0.001,0.001) dataset number", iter ,"\n")
  formula= y ~ X1+X2+X3+X4+f(feeders, model="iid", hyper=list(theta=list(prior="loggamma",param=c(0.001,0.001))))+ f(receivers, model="iid", hyper=list(theta=list(prior="loggamma",param=c(0.001,0.001))))
  result = inla(formula,family="binomial", data = dataset)
  estbeta[iter,]<-result$summary.fixed[,1]                    
  sdebeta[iter,]<-result$summary.fixed[,2]                    
  
  prec.marg<-(result$marginals.hyperpar$'Precision for feeders') 
  marg.variance<-inla.tmarginal(function(x) 1/x, prec.marg)  
  ifelse(sum(marg.variance[,2])==0,var[iter]<-0,var[iter]<-inla.emarginal(function(x) x, marg.variance))
  ifelse(sum(marg.variance[,2])==0,mm[iter]<-0,mm[iter]<-inla.emarginal(function(x) x^2, marg.variance))
  sd_var[iter]<-sqrt(mm[iter]-(var[iter])^2)                  
  
  prec.margB<-(result$marginals.hyperpar$'Precision for receivers') 
  marg.varianceB<-inla.tmarginal(function(x) 1/x, prec.margB)  
  ifelse(sum(marg.varianceB[,2])==0,varB[iter]<-0,varB[iter]<-inla.emarginal(function(x) x, marg.varianceB))
  ifelse(sum(marg.varianceB[,2])==0,mmB[iter]<-0,mmB[iter]<-inla.emarginal(function(x) x^2, marg.varianceB))
  sd_varB[iter]<-sqrt(mmB[iter]-(varB[iter])^2)  
  
  # INLA Ga(0.5,0.003737)
  cat(" INLA Ga(0.5,0.003737) dataset number", iter ,"\n")
  formula2= y ~ X1+X2+X3+X4+f(feeders, model="iid", hyper=list(theta=list(prior="loggamma",param=c(0.5,0.003737))))+ f(receivers, model="iid", hyper=list(theta=list(prior="loggamma",param=c(0.5,0.003737))))
  result2 = inla(formula2,family="binomial", data = dataset)
  estbeta2[iter,]<-result2$summary.fixed[,1]                    
  sdebeta2[iter,]<-result2$summary.fixed[,2]                    
  
  prec.marg2<-(result2$marginals.hyperpar$'Precision for feeders') 
  marg.variance2<-inla.tmarginal(function(x) 1/x, prec.marg2)  
  ifelse(sum(marg.variance2[,2])==0,var2[iter]<-0,var2[iter]<-inla.emarginal(function(x) x, marg.variance2))
  ifelse(sum(marg.variance2[,2])==0,mm2[iter]<-0,mm2[iter]<-inla.emarginal(function(x) x^2, marg.variance2))
  sd_var2[iter]<-sqrt(mm2[iter]-(var2[iter])^2)                  
  
  prec.margB2<-(result2$marginals.hyperpar$'Precision for receivers') 
  marg.varianceB2<-inla.tmarginal(function(x) 1/x, prec.margB2)  
  ifelse(sum(marg.varianceB2[,2])==0,varB2[iter]<-0,varB2[iter]<-inla.emarginal(function(x) x, marg.varianceB2))
  ifelse(sum(marg.varianceB2[,2])==0,mmB2[iter]<-0,mmB2[iter]<-inla.emarginal(function(x) x^2, marg.varianceB2))
  sd_varB2[iter]<-sqrt(mmB2[iter]-(varB2[iter])^2)  
  
  # MLLA
  cat(" MLLA dataset number", iter ,"\n")
  laplace<-glmer(y~ 1+ X1+X2+X3+X4+(1|feeders) + (1|receivers)  , dataset, family=binomial(link=logit),nAGQ = 1L)
  estbeta.lap[iter,]<-fixef(laplace)                   
  sdebeta.lap[iter,]<-sqrt(diag(vcov(laplace)))                         
  var.lap[iter,]<-getME(laplace, "theta")^2                       
  var_lapA[iter]<-var.lap[iter,1]
  var_lapB[iter]<-var.lap[iter,2]
  
  
  rm(dataset)                                                 
  
}

# INLA  Ga(0.001,0.001)
sigma<-mean(var)                                    
sd<-mean(sd_var)                                    
print(sigma)
sigmaB<-mean(varB)                                    
sdB<-mean(sd_varB)                                    
print(sigmaB)

for(i in 1:(betasize)){                                       
  betaaverage[i,]<-mean(estbeta[,i])
  sdeaverage[i]<-mean(sdebeta[,i])
  finaloutput[rowcount,(2*i-1):(2*i)]<-c(as.numeric(betaaverage[i]),as.numeric(sdeaverage[i]))
}  

# INLA  Ga(0.5,0.003737)
sigma2<-mean(var2)                                    
sd2<-mean(sd_var2)                                    
print(sigma2)
sigmaB2<-mean(varB2)                                    
sdB2<-mean(sd_varB2)                                    
print(sigmaB2)

for(i in 1:(betasize)){                                       
  betaaverage2[i,]<-mean(estbeta2[,i])
  sdeaverage2[i]<-mean(sdebeta2[,i])
  finaloutput2[rowcount,(2*i-1):(2*i)]<-c(as.numeric(betaaverage2[i]),as.numeric(sdeaverage2[i]))
}   

# MLLA
sigma.lap<-mean(var_lapA)
print(sigma.lap)
sigmaB.lap<-mean(var_lapB)
print(sigmaB.lap)

for(i in 1:(betasize)){   
  betaaverage.lap[i]<-mean(estbeta.lap[,i])
  sdeaverage.lap[i]<-mean(sdebeta.lap[,i])
  finaloutput.lap[rowcount,(2*i-1):(2*i)]<-c(as.numeric(betaaverage.lap[i]),as.numeric(sdeaverage.lap[i]))
}   


finaloutput<-as.data.frame(round(finaloutput,3))            # Fixed effects Estimates Inla Ga(0.001,0.001)
finaloutput2<-as.data.frame(round(finaloutput2,3))          # Fixed effects Estimates Inla Ga(0.5,0.003737)
finaloutput.lap<-as.data.frame(round(finaloutput.lap,3))    # Fixed effects Estimates MLLA


## relative bias random effects ###
#(net of aberrant estimates: zero estimates for MLLA, estimates>=2 for INLA)

# population parameters
sigmaua2=0.01   #(0.01,0.01,0.25,0.25,1.00)
sigmaub2=0.01   #(0.01,0.25,0.01,0.25,1.00)

rel.bias.reA.Inla<-(mean(subset(var,var<2))-sigmaua2)/sigmaua2     #Inla Ga(0.001,0.001)
rel.bias.reB.Inla<-(mean(subset(varB,varB<2))-sigmaub2)/sigmaub2


(length(subset(var,var>=2))/500)*100           #Percentage of aberrant estimates INLA Ga(0.001,0.001)
(length(subset(varB,varB>=2))/500)*100

rel.bias.reA.Inla2<-(mean(subset(var2,var2<2))-sigmaua2)/sigmaua2   #Inla Ga(0.5,0.003737)
rel.bias.reB.Inla2<-(mean(subset(varB2,varB2<2))-sigmaub2)/sigmaub2


(length(subset(var2,var2>=2))/500)*100         #Percentage of aberrant estimates INLA Ga(0.5,0.003737)
(length(subset(varB2,varB2>=2))/500)*100


rel.bias.reA.Laplace<-(mean(subset(var_lapA,var_lapA!=0))-sigmaua2)/sigmaua2  # MLLA
rel.bias.reB.Laplace<-(mean(subset(var_lapB,var_lapB!=0))-sigmaub2)/sigmaub2


(length(subset(var_lapA,var_lapA==0))/500)*100   #Percentage of Zero estimates for MLLA
(length(subset(var_lapB,var_lapB==0))/500)*100


cbind(rel.bias.reA.Inla2, rel.bias.reA.Inla ,rel.bias.reA.Laplace,
      rel.bias.reB.Inla2, rel.bias.reB.Inla, rel.bias.reB.Laplace )

### relative bias fixed effects ###

# INLA Ga(0.001,0.001)
fe.bugs<-matrix(cbind(finaloutput[1,1],finaloutput[1,3],finaloutput[1,5],finaloutput[1,7], finaloutput[1,9]),1,5)
# INLA Ga(0.5,0.003737)
fe.fong<-matrix(cbind(finaloutput2[1,1],finaloutput2[1,3],finaloutput2[1,5],finaloutput2[1,7], finaloutput2[1,9]),1,5)
# MLLA
fe.laplace<-matrix(cbind(finaloutput.lap[1,1],finaloutput.lap[1,3],finaloutput.lap[1,5],finaloutput.lap[1,7], finaloutput.lap[1,9]),1,5)

rel.bias.fe.bugs<-(fe.bugs-beta)/beta        # INLA Ga(0.001,0.001)
rel.bias.fe.fong<-(fe.fong-beta)/beta        # INLA Ga(0.5,0.003737)
rel.bias.fe.Laplace<-(fe.laplace-beta)/beta  # MLLA

cbind(rel.bias.fe.fong,rel.bias.fe.bugs,rel.bias.fe.Laplace)

### relative bias SE fixed effects ###

# INLA Ga(0.001,0.001)
se.bugs<-matrix(cbind(finaloutput[1,2],finaloutput[1,4],finaloutput[1,6],finaloutput[1,8], finaloutput[1,10]),1,5)
# INLA Ga(0.5,0.003737)
se.fong<-matrix(cbind(finaloutput2[1,2],finaloutput2[1,4],finaloutput2[1,6],finaloutput2[1,8], finaloutput2[1,10]),1,5)
# MLLA
se.laplace<-matrix(cbind(finaloutput.lap[1,2],finaloutput.lap[1,4],finaloutput.lap[1,6],finaloutput.lap[1,8], finaloutput.lap[1,10]),1,5)

sd.betaB<-matrix(cbind(sd(estbeta[,1]),sd(estbeta[,2]),sd(estbeta[,3]),sd(estbeta[,4]), sd(estbeta[,5])),1,5)
sd.betaF<-matrix(cbind(sd(estbeta2[,1]),sd(estbeta2[,2]),sd(estbeta2[,3]),sd(estbeta2[,4]), sd(estbeta2[,5])),1,5)
sd.betaL<-matrix(cbind(sd(estbeta.lap[,1]),sd(estbeta.lap[,2]),sd(estbeta.lap[,3]),sd(estbeta.lap[,4]), sd(estbeta.lap[,5])),1,5)

rel.bias.se.bugs<-(se.bugs-sd.betaB)/sd.betaB        # INLA Ga(0.001,0.001)
rel.bias.se.fong<-(se.fong-sd.betaF)/sd.betaF        # INLA Ga(0.5,0.003737)
rel.bias.se.Laplace<-(se.laplace-sd.betaL)/sd.betaL  # MLLA

cbind(rel.bias.se.fong,rel.bias.se.bugs,rel.bias.se.Laplace)

################################################################################################
######################## Check empty columns: ##############################################
################################################################################################

freq.receivers
freq.feeders

empty.col<-col.vuote<-matrix(0,sim,10)

for(i in 1:500){  #check the empty columns, zero means "empty coloumn"
  
  for(j in 1:20){
    for(k in 1:20){
    
    ifelse(freq.receivers[i,j]==freq.receivers[i,k],empty.col[i,freq.receivers[i,j]]<-freq.receivers[i,j] , empty.col[i,freq.receivers[i,j]]<-freq.receivers[i,j])
  }}
}


freq.emptycol<-matrix(0,sim,1)

for(i in 1:500){  # number of empty columns for the 500 replicates
  for(j in 1:10){
    ifelse(empty.col[i,j]==0, col.vuote[i,j]<-1, col.vuote[i,j]<-0 )
    }
  freq.emptycol[i,]<-sum(col.vuote[i,])
  }
freq.emptycol 
(sum(freq.emptycol)/5000)*100

hist(freq.emptycol,ylim=c(0,300),xlim=c(0,4),main="Empty columns",xlab="empty columns")  
