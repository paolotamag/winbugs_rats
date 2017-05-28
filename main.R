#Paolo Tamagnini
#1536242
#paolotamag@gmail.com

#########################################FUNCTIONS#############################################
histMY = function(chain,val,stringa1,stringa2){
  hist(chain,freq = F, breaks=75, col="orange",xlab = stringa1, main= stringa2)
  abline(v =val,col='red', lw=4) }

calc_moda = function(vattore){
  z <- density(vattore)
  return(z$x[z$y==max(z$y)])}

calc_MoreMode = function(vetofvett){
  foriter=dim(vetofvett)[2]
  modeVect = rep(0, foriter)
  for (i in 1:foriter) {
    modeVect[i]=calc_moda(vetofvett[,i]) }
  return(modeVect)}

threeHistMY = function(chain1,chain2,chainAll,val1,val2,val3,booL,mode,stringa1,stringa2){
  p1 <- hist(chain1,freq = F, breaks=75)
  p2 <- hist(chain2,freq = F, breaks=75)
  p3 <- hist(chainAll,freq = F, breaks=75)
  plot(p3,col=rgb(1,0,0,1/3), main=stringa1,xlab=stringa2, freq = F)  
  plot(p2,col=rgb(1,0.5,0,1/3),add=T,freq = F)  
  plot(p1,col=rgb(0,0,1,1/3),add=T,freq = F)  
  abline(v =val3,col='red',lw=4)
  abline(v =val2,col='orange',lw=4)
  abline(v =val1,col='blue',lw=4)
  if (booL) {
    abline(v =mode,col='green',lw=4) }}

calc_A = function(uai,bi,xB) {
  return(mean(uai)-bi*xB) }

calc_B = function(xi,uaib,xbu,nu) {
  Sxy = (sum(xi*uaib)-xbu*sum(uaib))
  Sxx = (sum(xi^2)-nu*xbu^2)
  return(Sxy/Sxx) }

calc_R2all_lin = function(uailol,a,b,xlol,xBarlol){
  R2tutti<-rep(0, 30)
  R2tuttiNames<-rep(0, 30)
  for (i in 1:30) {
    SSreg = sum((uailol[i,]-(a[i]+b[i]*(xlol-xBarlol)))^2)
    SStot = sum((uailol[i,]-mean(uailol[i,]))^2)
    R2tuttiNames[i] = paste('rat', i, sep = "")
    R2tutti[i] = 1 - SSreg/SStot }
  names(R2tutti) = R2tuttiNames
  return(R2tutti)}

calc_R2all_biv = function(uailol,a,b,xlol,xBarlol){
  R2tutti<-rep(0, 30)
  R2tuttiNames<-rep(0, 30)
  for (i in 1:30) {
    SSreg = sum((uailol[i,]-(a[i]+b[i]*xlol))^2)
    SStot = sum((uailol[i,]-mean(uailol[i,]))^2)
    R2tuttiNames[i] = paste('rat', i, sep = "")
    R2tutti[i] = 1 - SSreg/SStot }
  names(R2tutti) = R2tuttiNames
  return(R2tutti)}

calc_R2all_exp = function(uailol,a,b,g,xlol,xBarlol){
  R2tutti<-rep(0, 30)
  R2tuttiNames<-rep(0, 30)
  for (i in 1:30) {
    SSreg = sum((uailol[i,]-(a[i]-b[i]*g[i]^(xlol-xBarlol)))^2)
    SStot = sum((uailol[i,]-mean(uailol[i,]))^2)
    R2tuttiNames[i] = paste('rat', i, sep = "")
    R2tutti[i] = 1 - SSreg/SStot }
  names(R2tutti) = R2tuttiNames
  return(R2tutti)}

gamma_k <- function(k,vet,Icap) {
  t = length(vet)
  somma = 0
  for (ki in 1:(t-k)) {
    somma = somma + (vet[ki] - Icap)*(vet[ki+k] - Icap) }
  return(somma/(t-k))}


varCap <- function(vet,Icap) {
  somma = 0
  tVar = length(vet)
  print('Summing up all gamma_k, from k = 1 to:')
  print(tVar)
  print('Starting now..')
  print(paste('iter. ', 1, sep = ""))
  for (klol in 1:(tVar-1)) {
    if (klol%%500==0) {
      print(paste('iter. ', klol, sep = "")) }
    somma = somma + gamma_k(klol,vet,Icap) }
  num = gamma_k(0,vet,Icap) + 2*somma
  #print(num)
  return(num/tVar) }


###############################################################################################

nb=50000
nc= 2 #att se impalla poi
ni = 100000
nsim=ni-nb
library(R2jags)
set.seed(123)
setwd(getwd())

############################################################################################################################################
#################################################################### LINEAR ################################################################
############################################################################################################################################

d <- read.jagsdata("rats-data.R")
initsFile <- read.jagsdata("rats-init.R")
inits=function(){
  initsFile }
attributes(initsFile)$names
params = c("tau.c","alpha","tau.beta","beta","tau.alpha","beta.c","alpha.c")
ratJags=jags(data=d,inits=inits,parameters.to.save=params,model.file="rats_paolo.txt", 
             n.chains=nc, n.thin = 1, n.iter=ni, n.burnin =nb)

ratJags$model
ratJags$n.iter
alphaVect = ratJags$BUGSoutput$mean$alpha
betaVect = ratJags$BUGSoutput$mean$beta
alphaC = as.vector(ratJags$BUGSoutput$mean$alpha.c)
betaC = as.vector(ratJags$BUGSoutput$mean$beta.c)
xBar = mean(d$x)

X = t(replicate(30, d$x))
plot(X,d$Y, main='overall lin. reg.', ylab='weigth', xlab='days')
abline(a = alphaC-betaC*xBar, b=betaC, lw=2)


R2all = calc_R2all_lin(d$Y,alphaVect,betaVect,d$x,xBar)


alphaChain1 = ratJags$BUGSoutput$sims.matrix[1:nsim,1:30]
alphaChain2 = ratJags$BUGSoutput$sims.matrix[(nsim+1):(nsim*2),1:30]
alphaChain = ratJags$BUGSoutput$sims.matrix[,1:30]

betaChain1 = ratJags$BUGSoutput$sims.matrix[1:nsim,32:61]
betaChain2 = ratJags$BUGSoutput$sims.matrix[(nsim+1):(nsim*2),32:61]
betaChain = ratJags$BUGSoutput$sims.matrix[,32:61]


#the best fit is achieved for:##################################################
maxR2 = which.max(R2all)
maxR2Value = max(R2all)
i = unname(maxR2)

print('the best fit is achieved for:')
print(paste('rat_',i,sep = ""))


threeHistMY(alphaChain1[,i],alphaChain2[,i],alphaChain[,i],mean(alphaChain1[,i]),mean(alphaChain2[,i]),alphaVect[i],F,0,paste('all chains of alpha_',i,sep=""),paste('alpha_',i,sep=""))


threeHistMY(betaChain1[,i],betaChain2[,i],betaChain[,i],mean(betaChain1[,i]),mean(betaChain2[,i]),betaVect[i],F,0,paste('all chains of beta_',i,sep=""),paste('beta_',i,sep=""))

n = names(maxR2)
plot(d$x,d$Y[i,], main=n, ylab='weigth', xlab='days')
abline(a = alphaVect[i]-betaVect[i]*xBar, b=betaVect[i], col='red', lty='dashed',lw=1)

B=(sum(d$x*d$Y[i,])-xBar*sum(d$Y[i,]))/(sum(d$x^2)-5*xBar^2)
A=mean(d$Y[i,])-B*xBar

abline(a = A, b=B, col='blue', lty='dotted',lw=1)

text(x=15,y=300,label ='red line : MCMC approximation')
text(x=15,y=310,label ='blue line : unbiased estimator')

A
alphaVect[i]-betaVect[i]*xBar

B
betaVect[i]

#the worse fit is achieved for:##################################################
minR2 = which.min(R2all)
minR2Value = min(R2all)
i = unname(minR2)

print('the best fit is achieved for:')
print(paste('rat_',i,sep = ""))

threeHistMY(alphaChain1[,i],alphaChain2[,i],alphaChain[,i],mean(alphaChain1[,i]),mean(alphaChain2[,i]),alphaVect[i],F,0,paste('all chains of alpha_',i,sep=""),paste('alpha_',i,sep=""))

threeHistMY(betaChain1[,i],betaChain2[,i],betaChain[,i],mean(betaChain1[,i]),mean(betaChain2[,i]),betaVect[i],F,0,paste('all chains of beta_',i,sep=""),paste('beta_',i,sep=""))

B = calc_B(d$x,d$Y[i,],xBar,length(d$x))
A = calc_A(d$Y[i,],B,xBar)

nameR = names(minR2)
plot(d$x,d$Y[i,], main=nameR, ylab='weigth', xlab='days')
abline(a = alphaVect[i]-betaVect[i]*xBar, b=betaVect[i], col='red', lty='dashed',lw=1)
abline(a = A, b=B, col='blue', lty='dotted',lw=1)
text(x=15,y=300,label ='red line : MCMC approximation')
text(x=15,y=310,label ='blue line : unbiased estimator')

A
alphaVect[i]-betaVect[i]*xBar

B
betaVect[i]


traceplot(ratJags)

lb = 1
x = alphaChain[,i]
iter = (lb):(length(x)+(lb-1))
runningmeans=cumsum(x)/(1:length(x))
plot(iter,runningmeans,type="l",col='red',main='alpha[3] chain')

x = betaChain[,i]
iter = (lb):(length(x)+(lb-1))
runningmeans=cumsum(x)/(1:length(x))
plot(iter,runningmeans,type="l",col='red',main='beta[3] chain')

tooLong = length(betaChain[,i])
allErrA = rep(0,30)
allErrB = rep(0,30)
lastvalues = 1000
for (i in 1:30) {
  errA = varCap(alphaChain[,i][(tooLong-lastvalues):tooLong],alphaVect[i])
  allErrA[i] = sqrt(abs(errA))/alphaVect[i]
  errB = varCap(betaChain[,i][(tooLong-lastvalues):tooLong],betaVect[i])
  allErrB[i] = sqrt(abs(errB))/betaVect[i] }

mean(allErrA)*100
mean(allErrB)*100

############################################################################################################################################
################################################################## EXPONENTIAL #############################################################
############################################################################################################################################

initsFileExp_i <- read.jagsdata("rats-init_exp_i.R")
initsExp_i=function(){
  initsFileExp_i }
attributes(initsFileExp_i)$names
paramsExp_i = c("tau.c","alpha","tau.beta","beta","tau.alpha","beta.c","gamma","alpha.c")
ratJagsExp_i=jags(data=d,inits=initsExp_i,parameters.to.save=paramsExp_i,model.file="rats_paolo_exp_i.txt",
                  n.chains=nc,n.thin = 1,n.iter=ni,n.burnin =nb)

ratJagsExp_i$model

alphaVectExp_i = ratJagsExp_i$BUGSoutput$mean$alpha
betaVectExp_i = ratJagsExp_i$BUGSoutput$mean$beta
alphaCExp_i = as.vector(ratJagsExp_i$BUGSoutput$mean$alpha.c)
betaCExp_i = as.vector(ratJagsExp_i$BUGSoutput$mean$beta.c)
gammaExp_i = ratJagsExp_i$BUGSoutput$mean$gamma
xBar = mean(d$x)


i = 3

cols = 31
#ratJagsExp_i$BUGSoutput$sims.matrix[1,cols]
alpha_cChain1 = ratJagsExp_i$BUGSoutput$sims.matrix[1:nsim,cols]
alpha_cChain2 = ratJagsExp_i$BUGSoutput$sims.matrix[(nsim+1):(nsim*2),cols]
alpha_cChain = ratJagsExp_i$BUGSoutput$sims.matrix[,cols]
alphaCExp_i_mode = calc_moda(alpha_cChain)

threeHistMY(alpha_cChain1,alpha_cChain2,alpha_cChain,mean(alpha_cChain1),mean(alpha_cChain2),alphaCExp_i,T,alphaCExp_i_mode,'all chains of alpha_c','alpha_c')


cols = 62
#ratJagsExp_i$BUGSoutput$sims.matrix[1,cols]
beta_cChain1 = ratJagsExp_i$BUGSoutput$sims.matrix[1:nsim,cols]
beta_cChain2 = ratJagsExp_i$BUGSoutput$sims.matrix[(nsim+1):(nsim*2),cols]
beta_cChain = ratJagsExp_i$BUGSoutput$sims.matrix[,cols]
betaCExp_i_mode = calc_moda(beta_cChain)

threeHistMY(beta_cChain1,beta_cChain2,beta_cChain,mean(beta_cChain1),mean(beta_cChain2),betaCExp_i,T,betaCExp_i_mode,'all chains of beta_c','beta_c')

cols = 1:30
#ratJagsExp_i$BUGSoutput$sims.matrix[1,cols]
alphaChain1 = ratJagsExp_i$BUGSoutput$sims.matrix[1:nsim,cols]
alphaChain2 = ratJagsExp_i$BUGSoutput$sims.matrix[(nsim+1):(nsim*2),cols]
alphaChain = ratJagsExp_i$BUGSoutput$sims.matrix[,cols]
alphaVectExp_i_mode = calc_MoreMode(alphaChain)

threeHistMY(alphaChain1[,i],alphaChain2[,i],alphaChain[,i],mean(alphaChain1[,i]),mean(alphaChain2[,i]),alphaVectExp_i[i],T,alphaVectExp_i_mode[i],paste('all chains of alpha_',i,sep=""),paste('alpha_',i,sep=""))

text(x=520,y=0.020, labels ="green : mode")
text(x=520,y=0.022, labels ="blue : mean")


cols = 32:61
#ratJagsExp_i$BUGSoutput$sims.matrix[1,cols]
betaChain1 = ratJagsExp_i$BUGSoutput$sims.matrix[1:nsim,cols]
betaChain2 = ratJagsExp_i$BUGSoutput$sims.matrix[(nsim+1):(nsim*2),cols]
betaChain = ratJagsExp_i$BUGSoutput$sims.matrix[,cols]
betaVectExp_i_mode = calc_MoreMode(betaChain)

threeHistMY(betaChain1[,i],betaChain2[,i],betaChain[,i],mean(betaChain1[,i]),mean(betaChain2[,i]),betaVectExp_i[i],T,betaVectExp_i_mode[i],paste('all chains of beta_',i,sep=""),paste('beta_',i,sep=""))

text(x=250,y=0.020, labels ="green : mode")
text(x=250,y=0.022, labels ="blue : mean")

cols = 64:93
#ratJagsExp_i$BUGSoutput$sims.matrix[1,cols]
gammaChain1 = ratJagsExp_i$BUGSoutput$sims.matrix[1:nsim,cols]
gammaChain2 = ratJagsExp_i$BUGSoutput$sims.matrix[(nsim+1):(nsim*2),cols]
gammaChain = ratJagsExp_i$BUGSoutput$sims.matrix[,cols]
gammaExp_i_mode = calc_MoreMode(gammaChain)

threeHistMY(gammaChain1[,i],gammaChain2[,i],gammaChain[,i],mean(gammaChain1[,i]),mean(gammaChain2[,i]),gammaExp_i[i],T,gammaExp_i_mode[i],paste('all chains of gamma_',i,sep=""),paste('gamma_',i,sep=""))

text(x=0.975,y=100, labels ="green : mode")
text(x=0.975,y=110, labels ="blue : mean")



alphaCExp_i
alphaCExp_i_mode
betaCExp_i
betaCExp_i_mode
alphaVectExp_i[i]
alphaVectExp_i_mode[i]
betaVectExp_i[i]
betaVectExp_i_mode[i]
gammaExp_i[i]
gammaExp_i_mode[i]

plot(X,d$Y, main='overall lin. reg.', ylab='weigth', xlab='days')
abline(a = alphaC-betaC*xBar, b=betaC, col='red', lty='dashed',lw=1)
lines(X,alphaCExp_i_mode-betaCExp_i_mode*mean(gammaExp_i_mode)^(X-xBar),col="green",lw=2)
lines(X,alphaCExp_i-betaCExp_i*mean(gammaExp_i)^(X-xBar),col="blue",lty='dotted',lw=2)


i = 3

plot(d$x,d$Y[i,],m=paste("rat_",i,sep=""), ylab='weigth', xlab='days')
abline(a = alphaVect[i]-betaVect[i]*xBar, b=betaVect[i], col='red', lty='dashed',lw=1)
lines(d$x,alphaVectExp_i[i]-betaVectExp_i[i]*gammaExp_i[i]^(d$x-xBar),col="green",lw=2)
lines(d$x,alphaVectExp_i_mode[i]-betaVectExp_i_mode[i]*gammaExp_i_mode[i]^(d$x-xBar),col="blue",lty='dotted',lw=2)

text(x=15,y=300, labels ="green line : MCMC exp model - means - R2 = 0.9867631")
text(x=15,y=290,label ='red line : MCMC linear model - R2 = 0.9632383')
text(x=15,y=310,label ='blue line : MCMC exp model - modes - R2 = 0.9830632')


R2all_Exp_i = calc_R2all_exp(d$Y,alphaVectExp_i,betaVectExp_i,gammaExp_i,d$x,xBar)
R2all_Exp_i_mode = calc_R2all_exp(d$Y,alphaVectExp_i_mode,betaVectExp_i_mode,gammaExp_i_mode,d$x,xBar)

R2all[3]
R2all_Exp_i[3]
R2all_Exp_i_mode[3]

mean(R2all)
mean(R2all_Exp_i)
mean(R2all_Exp_i_mode)

lb = 1
x = alphaChain[,i]
iter = (lb):(length(x)+(lb-1))
runningmeans=cumsum(x)/(1:length(x))
plot(iter,runningmeans,type="l",col='red',main='alpha[3] chain')

x = betaChain[,i]
iter = (lb):(length(x)+(lb-1))
runningmeans=cumsum(x)/(1:length(x))
plot(iter,runningmeans,type="l",col='red',main='beta[3] chain')

x = gammaChain[,i]
iter = (lb):(length(x)+(lb-1))
runningmeans=cumsum(x)/(1:length(x))
plot(iter,runningmeans,type="l",col='red',main='gamma[3] chain')


################################################# R2 and DIC plots #####################################

db <- read.jagsdata("birats-data.R")
initsFileb <- read.jagsdata("birats-inits.R")
initsb=function(){
  initsFileb }
attributes(initsFileb)$names
paramsb = c("tau.c", "tau.beta", "Omega.beta", "beta", "mu.beta")

maxIter = 50000
nR = 100
nPerc = 3
R2allExp_i_M = list()
R2all_Exp_i_mode_M = list()
DICmodelExp_M = list()
DICmodelLin_M = list()
R2allLin_M = list()
R2allLinBiv_M = list()
DICmodelBiv_M = list()

for (j in 1:3) {
  print('j:')
  print(j)
  listPerc = c(10,25,50)
  R2allExp_i_O = rep(0,nR)
  R2all_Exp_i_mode_O = rep(0,nR)
  DICmodelExp = rep(0,nR)
  DICmodelLin = rep(0,nR)
  R2allLin_O = rep(0,nR)
  R2allLinBiv_O = rep(0,nR)
  DICmodelBiv = rep(0,nR)
  perc = listPerc[j]/100
  
  for (i in 1:nR) {
    
    print('i:')
    print(i)
    
    ni=maxIter*i/nR
    
    nb = perc*ni
    
    ratJags=jags(data=d,inits=inits,parameters.to.save=params,model.file="rats_paolo.txt", 
                 n.chains=nc, n.thin = 1, n.iter=ni, n.burnin =nb)
    
    a = ratJags$BUGSoutput$mean$alpha
    b = ratJags$BUGSoutput$mean$beta
    
    R2allLin_O[i]=mean(calc_R2all_lin(d$Y,a,b,d$x,xBar))
    
    ratJagsExp_i=jags(data=d,inits=initsExp_i,parameters.to.save=paramsExp_i,model.file="rats_paolo_exp_i.txt",
                      n.chains=nc,n.thin = 1,n.iter=ni,n.burnin =nb)
    
    a = ratJagsExp_i$BUGSoutput$mean$alpha
    b = ratJagsExp_i$BUGSoutput$mean$beta
    g = ratJagsExp_i$BUGSoutput$mean$gamma
    
    R2allExp_i_O[i]=mean(calc_R2all_exp(d$Y,a,b,g,d$x,xBar))
    
    cols = 1:30
    alphaChain = ratJagsExp_i$BUGSoutput$sims.matrix[,cols]
    cols = 32:61
    betaChain = ratJagsExp_i$BUGSoutput$sims.matrix[,cols]
    cols = 64:93
    gammaChain = ratJagsExp_i$BUGSoutput$sims.matrix[,cols]
    
    am = calc_MoreMode(alphaChain)
    bm = calc_MoreMode(betaChain)
    gm = calc_MoreMode(gammaChain)
    
    R2all_Exp_i_mode_O[i]=mean(calc_R2all_exp(d$Y,am,bm,gm,d$x,xBar))
    
    ratJagsb=jags(data=db,inits=initsb,parameters.to.save=paramsb,model.file="birats4.bug",
                  n.chains=nc, n.thin = 1, n.iter=ni, n.burnin =nb)
    
    a = ratJagsb$BUGSoutput$mean$beta[,1]
    b = ratJagsb$BUGSoutput$mean$beta[,2]
    
    R2allLinBiv_O[i]=mean(calc_R2all_biv(d$Y,a,b,d$x,xBar))
    
    DICmodelBiv[i] = ratJagsb$BUGSoutput$DIC
    DICmodelExp[i] = ratJagsExp_i$BUGSoutput$DIC
    DICmodelLin[i] = ratJags$BUGSoutput$DIC }  
  
  R2allExp_i_M = c(R2allExp_i_M,list(R2allExp_i_O))
  R2allLinBiv_M = c(R2allLinBiv_M,list(R2allLinBiv_O))
  R2allLin_M = c(R2allLin_M,list(R2allLin_O))
  R2all_Exp_i_mode_M = c(R2all_Exp_i_mode_M,list(R2all_Exp_i_mode_O))
  
  DICmodelExp_M = c(DICmodelExp_M,list(DICmodelExp))
  DICmodelLin_M = c(DICmodelLin_M,list(DICmodelLin))
  DICmodelBiv_M = c(DICmodelBiv_M,list(DICmodelBiv)) }


titoloR = c('R2 at burnin = 10 %',"R2 at burnin = 25 %","R2 at burnin = 50 %")
titoloDIC = c('DIC at burnin = 10 %',"DIC at burnin = 25 %","DIC at burnin = 50 %")

klol = 3

RLIN = R2allLin_M[[klol]]

REXP = R2allExp_i_M[[klol]]

REXP_MODE = R2all_Exp_i_mode_M[[klol]]

RLIN_BIV =R2allLinBiv_M[[klol]]

DIC_LIN = DICmodelLin_M[[klol]]

DIC_EXP = DICmodelExp_M[[klol]]

DIC_LIN_BIV = DICmodelBiv_M[[klol]]

xfit = (1:nR)*maxIter/nR

yrange<-c(0.8,1)
xch = 20000
plot(xfit,RLIN,type="l",ylim=yrange,col='black',ylab='R2', xlab = '# of iterations', main=titoloR[klol])
lines(xfit,REXP,type="l",col='blue',lty='dashed')
lines(xfit,REXP_MODE,type="l",col='green',lty='dashed')
lines(xfit,RLIN_BIV,type="l",col='red',lty='dotted',lw=2)

text(x=xch, y=0.9, pos=4, labels=c('black line : avg R2 Linear Model'))
text(x=xch, y=0.88, pos=4, labels=c('blue line : avg R2 Exponential Model means'))
text(x=xch, y=0.87, pos=4, labels=c('green line : avg R2 Exponential Model modes'))
text(x=xch, y=0.89, pos=4, labels=c('red line : avg R2 Bivariate Model'))


yrange<-c(min(DIC_LIN,DIC_EXP,DIC_LIN_BIV),10000)
xch = 35000
plot(xfit,DIC_LIN,type="l",ylim=yrange,col='black',ylab='DIC', xlab = '# of iterations', main=titoloDIC[klol])
lines(xfit,DIC_EXP,type="l",col='blue',lty='dashed')
lines(xfit,DIC_LIN_BIV,type="l",col='red',lty='dotted',lw=2)

text(x=xch, y=6000, pos=4, labels=c('black line : DIC Linear Model'))
text(x=xch, y=6500, pos=4, labels=c('blue line : DIC Exponential Model'))
text(x=xch, y=7000, pos=4, labels=c('red line : DIC Bivariate Model'))
