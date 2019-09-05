#####SMC####

D = as.matrix((read.table("DistMatrix.txt")))[18:20,18:20]
N = as.matrix.data.frame(read.csv("UKpop.csv"))[,c(1,19:21)]
B = as.matrix.data.frame(read.csv("UKbirths.csv"))[,c(1,19:21)]
y = data.frame(read.csv("UKmeasles.csv"))[,c(1,19:21)]
year = y[,1]

school = 4
M = 130 #times (biweekly reporting so it is 130 biweeks)
J = 500 #particles
X = array(0,dim = c((M+1),J,3*5))
X.true = array(0,dim = c((M+1),3*5))
proposed_y = array(0,dim = c(M,3))

# Fixed parameters
sigma2 = 0.08
c = 0.4
year_0 = 1944
entryage = school
school_day = 251
tol = 10^(-6)

T = ((entryage+1) * 52 + 2)/2 # start_year:1949 (there are 262 weeks from 1944 to 1949)

# holiday effect
p = 0.739
Beta.s = function(R0,v_ir,a) (1+2*(1-p)*a)*R0*v_ir
h = function(a) (1-2*p*a)/(1+2*(1-p)*a) #holiday coefficient
h_eff = function(t){
  d = (year[t]-year_0) %% 1 * 365.25
  if (d<7){1}
  else if (d>99 & d<116){1}
  else if (d>198 & d<253){1}
  else if (d>299 & d<309){1}
  else if (d>355){1}
  else {0}
}

# birth effect
r = function(k,t) { 
  d = function(t) (year[t]-year_0) %%  1*365.25
  if(d(t) < school_day && school_day < d(t+1)){
    c * as.numeric(B[round(year[t]-entryage-year_0+0.51),(k+1)])
  } else {(1-c) * as.numeric(B[round(year[t]-year_0-entryage+0.51),(k+1)]/26)}
}

# measurement model 
dmeas <- function(y,rho,phi, H) {
  if(y>=0.5){
    pnorm(y+0.5,mean = rho*H,sd = sqrt(rho*(1-rho)*H + (phi*rho*H)^2+1)) - pnorm(y-0.5,mean = rho*H,sd = sqrt(rho*(1-rho)*H + (phi*rho*H)^2+1))
  }else {pnorm(0.5,mean = rho*H,sd = sqrt(rho*(1-rho)*H + (phi*rho*H)^2+1))}
}

rmeas <- function(rho,phi, H) {r=mean(round(rnorm(100,mean = rho*H,sd = sqrt(rho*(1-rho)*H + (phi*rho*H)^2+1))))
if(r<0){0
}else{r}}


#### transmission model ####
seir <- function(t,R0,a,alpha,mu,v_ei,v_ir,sigma2,rho,phi,G,X){
  for (k in 1:3) {
    vkl = numeric(3)
    ss = numeric(3)
    population_k = as.numeric(N[round(year[t]-year_0+0.51),(k+1)])
    # N[round(year[t]-year_0+0.5),]
    # time  LIVERPOOL     LONDON MANCHESTER 
    # 1949     800921    3375914     699712
    for(l in 1:3){
      population_l = as.numeric(N[round(year[t]-year_0+0.51),(l+1)])
      if(l==k){
        vkl[l] = 0
        ss[l] = 0
      }else{
        vkl[l] = round(G * (mean(D)/((mean(N[round(year[t]-year_0+0.51),]))^2)) *  
                         (population_k * population_l/D[k,l]))
        ss[l] = vkl[l]/population_k*(
          (X[m,j,3+(l-1)*5]/population_l)^alpha - (X[m,j,3+(k-1)*5]/population_k)^alpha)
        
      }}
    # tot_travel = sum(vkl)
    # travel_prob = 1- exp(-tot_travel/as.numeric(N[round(year[t]-year_0+0.5),(k+1)]))
    # 
    # for (l in 1:3) {
    #   if(l==k){
    #     ss[l] = 0
    #   }else{
    #     ss[l] = travel_prob * vkl[l]/tot_travel *
    #       ((X[m,j,3+(l-1)*5]/as.numeric(N[round(year[t]-year_0+0.5),(l+1)]))^alpha - 
    #          (X[m,j,3+(k-1)*5]/as.numeric(N[round(year[t]-year_0+0.5),(k+1)]))^alpha)
    #   }
    # }
    
    S = X[m,j,1+(k-1)*5]
    E = X[m,j,2+(k-1)*5]
    I = X[m,j,3+(k-1)*5]
    N_ir = X[m,j,4+(k-1)*5]
    P = X[m,j,5+(k-1)*5]
    
    v_se = ifelse(h_eff(t)==1,h(a),1) * Beta.s(R0,v_ir,a)  * ((I /P)^alpha + sum(ss))
    
    gamma_SE = rgamma(1,shape = 1/sigma2,scale = sigma2)
    gamma_SD = rgamma(1,shape = 1/sigma2,scale = sigma2)
    gamma_EI = rgamma(1,shape = 1/sigma2,scale = sigma2)
    gamma_ED = rgamma(1,shape = 1/sigma2,scale = sigma2)
    gamma_IR = rgamma(1,shape = 1/sigma2,scale = sigma2)
    gamma_ID = rgamma(1,shape = 1/sigma2,scale = sigma2)
    
    S_out = v_se * gamma_SE + mu * gamma_SD
    prob_se = v_se * gamma_SE/  Sout * (1-exp(-S_out))
    prob_sd = mu * gamma_SD/ S_out * (1-exp(-S_out))
    dN_se = rbinom(1, S, prob_se)
    dN.sd= ifelse(S-dN_se==0,0,rbinom(1, S-dN_se, prob_sd/(exp(-S_out)+prob_sd)))
    E_out = v_ei * gamma_EI + mu * gamma_ED
    p_ei = v_ei * gamma_EI/  E_out * (1-exp(-E_out))
    prob_ed = mu * gamma_ED/ E_out * (1-exp(-E_out))
    dN_ei = rbinom(1, E, p_ei)
    dN.ed = ifelse(E-dN_ei==0,0,rbinom(1, E-dN_ei, prob_ed/(exp(-E_out)+prob_ed))) 
    I_out = v_ir * gamma_IR + mu * gamma_ID
    prob_ir = v_ir * gamma_IR/  I_out * (1-exp(-I_out))
    prob_id = mu * gamma_ID/ I_out * (1-exp(-I_out))
    dN_ir = rbinom(1, I, prob_ir)
    dN.id = ifelse(I-dN_ir==0,0,rbinom(1, I-dN_ir, prob_id/(exp(-I_out)+prob_id))) 
    
    S_new = S + round(r(k,t)) - dN_se - dN.sd
    E_new = E + dN_se - dN_ei - dN.ed
    I_new = I + dN_ei - dN_ir - dN.id
    N_ir_new = dN_ir
    P_new  = N[round(year[T]-year_0+0.5),(k+1)]
    
    X[m+1,j,1+(k-1)*5] = S_new
    X[m+1,j,2+(k-1)*5] = E_new
    X[m+1,j,3+(k-1)*5] = I_new
    X[m+1,j,4+(k-1)*5] = N_ir_new
    X[m+1,j,5+(k-1)*5] = P_new
  }
  X[m+1,j,]
}

# define parameters and initial parameters
R0 = 25
a = 0.163
G = 500   
alpha = 0.97 
v_ei = 1
v_ir = 1     
mu = 0.0166/52
rho = 0.5  
phi = 0.25 
s = 0.05
e = 0.002
i = 0.003

####proposed X and y####
j = 1
X.true = array(0,dim = c((M+1),1,15))
for (k in 1:3) {
  X.true[1,1,1+(k-1)*5] = round(s*N[round(year[T]-year_0+0.5),(k+1)]) #S
  X.true[1,1,2+(k-1)*5] = round(e*N[round(year[T]-year_0+0.5),(k+1)]) #E
  X.true[1,1,3+(k-1)*5] = round(i*N[round(year[T]-year_0+0.5),(k+1)]) #I
  X.true[1,1,4+(k-1)*5] = 0   #delta N_ir(Biweekly)
  X.true[1,1,5+(k-1)*5] = N[round(year[T]-year_0+0.5),(k+1)] #P
}

for (m in 1:M) {
  X.true[m+1,1,] = seir(T + m -1,R0,a,alpha,mu,v_ei,v_ir,sigma2,rho,phi,G,X.true)
  proposed_y[m,] = sapply(1:3, function(k) round(rmeas(rho,phi,X.true[m+1,1,4+(k-1)*5])))
}


set.seed(2018)
w = array(0,dim = c(M,J))
w.o = array(0,dim = c(M,J))
tolerance = 10^(-30)

####initial X####
for (k in 1:3) {
  X[1,,1+(k-1)*5] = rep(round(s*N[round(year[T]-year_0+0.5),(k+1)]),J) #S
  X[1,,2+(k-1)*5] = rep(round(e*N[round(year[T]-year_0+0.5),(k+1)]),J) #E
  X[1,,3+(k-1)*5] = rep(round(i*N[round(year[T]-year_0+0.5),(k+1)]),J) #I
  X[1,,4+(k-1)*5] = rep(0,J)   #delta N_ir(Biweekly)
  X[1,,5+(k-1)*5] = rep(N[round(year[T]-year_0+0.5),(k+1)],J) #P 
}
# N[6,] population at 1949

#### iteration ####
for(m in 1:M){
  t = T+m-1 #biweekly reported round(T+m*0.5-0.5)
  for (j in 1:J) {
    #Proposed X using proposed parameters
    X[m+1,j,] = seir(t,R0,a,alpha,mu,v_ei,v_ir,sigma2,rho,phi,G,X)
    w.o[m,j] = prod(sapply(1:3,function(k) dmeas(proposed_y[m,k],rho,phi, X[m+1,j,4+(k-1)*5])))}
  #calculate true weight
  for (j in 1:J) {
    w[m,j] = ifelse(w.o[m,j]==0,0,w.o[m,j]/sum(w.o[m,]))
    #Filtering failure
    if(w[m,j] > tolerance){
      w[m,j] = w[m,j]
    }else{
      w[m,j] = tolerance
    }}
  #Resample
  number = as.numeric(rmultinom(1,J,w[m,]))
  number[0] = 0
  for(ii in 1:J){
    if((number[ii]==0)==FALSE){
      for(iii in 1:number[ii]){
        n = sum(number[0:(ii-1)])
        for (k in 1:3) {
          X[m+1,iii+n,1+(k-1)*5]=X[m+1,ii,1+(k-1)*5]
          X[m+1,iii+n,2+(k-1)*5]=X[m+1,ii,2+(k-1)*5]
          X[m+1,iii+n,3+(k-1)*5]=X[m+1,ii,3+(k-1)*5]
          X[m+1,iii+n,4+(k-1)*5]=X[m+1,ii,4+(k-1)*5]
          X[m+1,iii+n,5+(k-1)*5]=X[m+1,ii,5+(k-1)*5]
        }
      }}}
}

log_likelihood = sum(log(rowMeans(w.o)))

####London####
par(mfrow=c(2,2))
plot(rowMeans(X[,,6]),xlim=c(0,140),ylim=c(0,200000),type = 'l',main = 'London',xlab = 'Times(Biweek)',ylab='Susceptible Population')
lines(sapply(1:M, function(m) median(X[m,,6])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,,6],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,,6],0.9)),col=3,lty=2)
lines(X.true[,,6],type='l',col=4,lty=1)
legend("topright",inset=0.02,bty="n",legend=c("Mean","Median",'10%,90% quantile',"True S"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.5)

plot(rowMeans(X[,,7]),xlim=c(0,140),ylim=c(0,40000),type = 'l',main = 'London',xlab = 'Times(Biweek)',ylab='Exposed Population',cex=0.1)
lines(sapply(1:M, function(m) median(X[m,,7])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,,7],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,,7],0.9)),col=3,lty=2)
lines(X.true[,,7],type='l',col=4,lty=1)
legend("topleft", bty="n",legend=c("Mean","Median",'10%,90% quantile',"True E"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.75)

plot(rowMeans(X[,,8]),xlim=c(0,140),ylim=c(0,35000),type = 'l',main = 'London',xlab = 'Times(Biweek)',ylab='Infectious Population')
lines(sapply(1:M, function(m) median(X[m,,8])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,,8],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,,8],0.9)),col=3,lty=2)
lines(X.true[,,8],type='l',col=4,lty=1)
legend("topleft", bty="n",legend=c("","Median",'10%,90% quantile',"True I"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)

plot(rowMeans(X[,,9]),xlim=c(0,140),ylim = c(0,25000),type = 'l',main = 'London',xlab = 'Times(Biweek)',ylab='BiWeekly Recovered')
lines(sapply(1:M, function(m) median(X[m,,9])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,,9],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,,9],0.9)),col=3,lty=2)
lines(X.true[,,9],type='l',col=4,lty=1)
legend("topleft", bty="n",legend=c("Mean","Median",'10%,90% quantile',"True N_ir"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)

par(mfrow=c(2,2))
plot(sapply(1:M,function(m) var(X[m,,6])),type='l',main='London',xlab = 'Times(Biweek)',ylab='Variance of Susceptible Population')
plot(sapply(1:M,function(m) var(X[m,,7])),type='l',main='London',xlab = 'Times(Biweek)',ylab='Variance of Exposed Population')
plot(sapply(1:M,function(m) var(X[m,,8])),type='l',main='London',xlab = 'Times(Biweek)',ylab='Variance of Infectious Population')
plot(sapply(1:M,function(m) var(X[m,,9])),type='l',main='London',xlab = 'Times(Biweek)',ylab='Variance of Biweekly Recovered')

####Liverpool####
par(mfrow=c(2,2))
plot(rowMeans(X[,,1]),xlim=c(0,140),type = 'l',main = 'Liverpool',xlab = 'Times(Biweek)',ylab='Susceptible Population')
lines(sapply(1:M, function(m) median(X[m,,1])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,,1],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,,1],0.9)),col=3,lty=2)
lines(X.true[,,1],type='l',col=4,lty=1)
legend("topleft", bty="n",legend=c("Mean","Median",'10%,90% quantile',"True S"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)

plot(rowMeans(X[,,2]),xlim=c(0,140),ylim=c(0,12000),type = 'l',main = 'Liverpool',xlab = 'Times(Biweek)',ylab='Exposed Population')
lines(sapply(1:M, function(m) median(X[m,,2])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,,2],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,,2],0.9)),col=3,lty=2)
lines(X.true[,,2],type='l',col=4,lty=1)
legend("topleft", bty="n",legend=c("Mean","Median",'10%,90% quantile',"True E"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)

plot(rowMeans(X[,,3]),xlim=c(0,140),ylim=c(0,12000),type = 'l',main = 'Liverpool',xlab = 'Times(Biweek)',ylab='Infectious Population')
lines(sapply(1:M, function(m) median(X[m,,3])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,,3],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,,3],0.9)),col=3,lty=2)
lines(X.true[,,3],type='l',col=4,lty=1)
legend("topleft", bty="n",legend=c("Mean","Median",'10%,90% quantile',"True I"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)

plot(rowMeans(X[,,4]),xlim=c(0,140),ylim=c(0,8000),type = 'l',main = 'Liverpool',xlab = 'Times(Biweek)',ylab='BiWeekly Recovered')
lines(sapply(1:M, function(m) median(X[m,,4])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,,4],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,,4],0.9)),col=3,lty=2)
lines(X.true[,,4],type='l',col=4,lty=1)
legend("topleft", bty="n",legend=c("Mean","Median",'10%,90% quantile',"True N_ir"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)

par(mfrow=c(2,2))
plot(sapply(1:M,function(m) var(X[m,,1])),type='l',main='Liverpool',xlab = 'Times(Biweek)',ylab='Variance of Susceptible Population')
plot(sapply(1:M,function(m) var(X[m,,2])),type='l',main='Liverpool',xlab = 'Times(Biweek)',ylab='Variance of Exposed Population')
plot(sapply(1:M,function(m) var(X[m,,3])),type='l',main='Liverpool',xlab = 'Times(Biweek)',ylab='Variance of Infectious Population')
plot(sapply(1:M,function(m) var(X[m,,4])),type='l',main='Liverpool',xlab = 'Times(Biweek)',ylab='Variance of Biweekly Recovered')


####Manchester####
par(mfrow=c(2,2))
plot(rowMeans(X[,,11]),xlim=c(0,140),type = 'l',main = 'Manchester',xlab = 'Times(Biweek)',ylab='Susceptible Population')
lines(sapply(1:M, function(m) median(X[m,,11])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,,11],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,,11],0.9)),col=3,lty=2)
lines(X.true[,,11],type='l',col=4,lty=1)
legend("topleft", bty="n",legend=c("Mean","Median",'10%,90% quantile',"True S"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)

plot(rowMeans(X[,,12]),xlim=c(0,140),ylim=c(0,8000),type = 'l',main = 'Manchester',xlab = 'Times(Biweek)',ylab='Exposed Population')
lines(sapply(1:M, function(m) median(X[m,,12])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,,12],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,,12],0.9)),col=3,lty=2)
lines(X.true[,,12],type='l',col=4,lty=1)
legend("topleft", bty="n",legend=c("Mean","Median",'10%,90% quantile',"True E"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)

plot(rowMeans(X[,,13]),xlim=c(0,140),ylim=c(0,7000),type = 'l',main = 'Manchester',xlab = 'Times(Biweek)',ylab='Infectious Population')
lines(sapply(1:M, function(m) median(X[m,,13])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,,13],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,,13],0.9)),col=3,lty=2)
lines(X.true[,,13],type='l',col=4,lty=1)
legend("topleft", bty="n",legend=c("Mean","Median",'10%,90% quantile',"True I"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)

plot(rowMeans(X[,,14]),xlim=c(0,140),ylim=c(0,5000),type = 'l',main = 'Manchester',xlab = 'Times(Biweek)',ylab='BiWeekly Recovered')
lines(sapply(1:M, function(m) median(X[m,,14])),col=2,lty=1)
lines(sapply(1:M, function(m) quantile(X[m,,14],0.1)),col=3,lty=2)
lines(sapply(1:M, function(m) quantile(X[m,,14],0.9)),col=3,lty=2)
lines(X.true[,,14],type='l',col=4,lty=1)
legend("topleft", bty="n",legend=c("Mean","Median",'10%,90% quantile',"True N_ir"),col=c(1,2,3,4), lty=c(1,1,2,1), cex=0.7)

par(mfrow=c(2,2))
plot(sapply(1:M,function(m) var(X[m,,11])),type='l',main='Manchester',xlab = 'Times(Biweek)',ylab='Variance of Susceptible Population')
plot(sapply(1:M,function(m) var(X[m,,12])),type='l',main='Manchester',xlab = 'Times(Biweek)',ylab='Variance of Exposed Population')
plot(sapply(1:M,function(m) var(X[m,,13])),type='l',main='Manchester',xlab = 'Times(Biweek)',ylab='Variance of Infectious Population')
plot(sapply(1:M,function(m) var(X[m,,14])),type='l',main='Manchester',xlab = 'Times(Biweek)',ylab='Variance of Biweekly Recovered')

par(mfrow=c(1,1))
plot(proposed_y[,1],xlim=c(0,140),ylim=c(0,15000),type = 'l',main = 'Proposed Data',xlab = 'Times(Biweek)',ylab='Biweekly reported cases',col=1)
lines(proposed_y[,2],col=2)
lines(proposed_y[,3],col=3)
legend("topright", legend=c("Liverpool",'London','Manchester'),col=c(1,2,3), lty=c(1,1,1), cex=0.7)


