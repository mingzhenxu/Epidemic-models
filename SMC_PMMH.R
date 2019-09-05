#####SMC PMMH####

D= as.matrix((read.table("DistMatrix.txt")))[18:20,18:20]
N = as.matrix.data.frame(read.csv("UKpop.csv"))[,c(1,19:21)]
B = as.matrix.data.frame(read.csv("UKbirths.csv"))[,c(1,19:21)]
y = data.frame(read.csv("UKmeasles.csv"))[,c(1,19:21)]
year = y[,1]

school = 4
M = 10 #times (biweekly reporting so it is 130 biweeks)
J = 50 #particles
I = 50 #number of iterations
X = array(0,dim = c((I+1),(M+1),J,3*5))
#X.true = array(0,dim = c((M+1),3*5))
#proposed_y = array(0,dim = c(M,3))

# Fixed parameters
sigma2 = 0.08
c = 0.4
year_0 = 1944
entry_age = school
school.start.day = 251
tol = 10^(-6)

T = ((entry_age+1) * 52 + 2)/2 # start_year:1949 (there are 262 weeks from 1944 to 1949)

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
  if(d(t) < school.start.day && school.start.day < d(t+1)){
    c * as.numeric(B[round(year[t]-entry_age-year_0+0.51),(k+1)])
  } else {(1-c) * as.numeric(B[round(year[t]-year_0-entry_age+0.51),(k+1)]/26)}
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
seir <- function(i,t,R0,a,alpha,mu,v_ei,v_ir,sigma2,rho,phi,G,X){
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
          (X[i,m,j,3+(l-1)*5]/population_l)^alpha - (X[i,m,j,3+(k-1)*5]/population_k)^alpha)
        
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
    
    S = X[i,m,j,1+(k-1)*5]
    E = X[i,m,j,2+(k-1)*5]
    I = X[i,m,j,3+(k-1)*5]
    N_ir = X[i,m,j,4+(k-1)*5]
    P = X[i,m,j,5+(k-1)*5]
    
    v_se = ifelse(h_eff(t)==1,h(a),1) * Beta.s(R0,v_ir,a)  * ((I /P)^alpha + sum(ss))
    
    gamma_SE = rgamma(1,shape = 1/sigma2,scale = sigma2)
    gamma_SD = rgamma(1,shape = 1/sigma2,scale = sigma2)
    gamma_EI = rgamma(1,shape = 1/sigma2,scale = sigma2)
    gamma_ED = rgamma(1,shape = 1/sigma2,scale = sigma2)
    gamma_IR = rgamma(1,shape = 1/sigma2,scale = sigma2)
    gamma_ID = rgamma(1,shape = 1/sigma2,scale = sigma2)
    
    S_out = v_se * gamma_SE + mu * gamma_SD
    p_se = v_se * gamma_SE/  S_out * (1-exp(-S_out))
    prob_sd = mu * gamma_SD/ S_out * (1-exp(-S_out))
    dN_se = rbinom(1, S, p_se)
    dN.sd= ifelse(S-dN_se==0,0,rbinom(1, S-dN_se, prob_sd/(exp(-S_out)+prob_sd))) 
    
    E_out = v_ei * gamma_EI + mu * gamma_ED
    p_ei = v_ei * gamma_EI/  E_out * (1-exp(-E_out))
    prob_ed = mu * gamma_ED/ E_out * (1-exp(-E_out))
    dN_ei = rbinom(1, E, p_ei)
    dN.ed = ifelse(E-dN_ei==0,0,rbinom(1, E-dN_ei, prob_ed/(exp(-E_out)+prob_ed))) 
    
    I_out = v_ir * gamma_IR + mu * gamma_ID
    p_ir = v_ir * gamma_IR/  I_out * (1-exp(-I_out))
    prob_id = mu * gamma_ID/ I_out * (1-exp(-I_out))
    dN_ir = rbinom(1, I, p_ir)
    dN.id = ifelse(I-dN_ir==0,0,rbinom(1, I-dN_ir, prob_id/(exp(-I_out)+prob_id))) # because prob_id/(1-p_ir) = prob_id/(exp(-I_out)+prob_id)
    
    S_new = S + round(r(k,t)) - dN_se - dN.sd
    E_new = E + dN_se - dN_ei - dN.ed
    I_new = I + dN_ei - dN_ir - dN.id
    N_ir_new = dN_ir
    prob_new  = N[round(year[T]-year_0+0.5),(k+1)]
    
    X[i,m+1,j,1+(k-1)*5] = S_new
    X[i,m+1,j,2+(k-1)*5] = E_new
    X[i,m+1,j,3+(k-1)*5] = I_new
    X[i,m+1,j,4+(k-1)*5] = N_ir_new
    X[i,m+1,j,5+(k-1)*5] = prob_new
  }
  X[i,m+1,j,]
}

set.seed(2019)

#### create theta matrix ####
theta.1 = array(0,dim = c((I+1),(M+1),J))
theta.2 = array(0,dim = c((I+1),(M+1),J))
theta.3 = array(0,dim = c((I+1),(M+1),J))

#### Parametrize ####
theta.1[1,1,] = rnorm(J,-3,1)
theta.2[1,1,] = sapply(1:J,function(j) min(rnorm(1,-8,1),-theta.1[1,1,j]))
theta.3[1,1,] = sapply(1:J,function(j) min(rnorm(1,-8,1),
                                            (log(1-exp(theta.1[1,1,j]+theta.2[1,1,j])) - 
                                            log(exp(theta.1[1,1,j])+exp(theta.2[1,1,j])+2*exp(theta.1[1,1,j]+theta.2[1,1,j])) )))

logit = function(x) 1/(1+exp(-x))

# Initializations of parameters
R0 = 25
a = 0.163
G = 500   
alpha = 0.97 
v_ei = 1
v_ir = 1     
mu = 0.0166/52
rho = 0.5  
phi = 0.25 

w = array(0,dim = c(I,M,J))
w.o = array(0,dim = c(I,M,J))
tolerance = 10^(-30)
# s = 0.05
# e = 0.002
# i = 0.003

#### PMMH ####
####initialize####
for (iter in 1) {
  s = sapply(1:J,function(j) logit(theta.1[iter,1,j]))
  e = sapply(1:J,function(j) logit(theta.2[iter,1,j]))
  i = sapply(1:J,function(j) logit(theta.3[iter,1,j]))
  log_likelihood = rep(0,I+1)
  
  for(m in 1:M){ 
    #Propsed particles
    for (j in 1:J) {
      t = T + m -1
      if(m==1){
        # initial X
        for (k in 1:3) {
          X[iter,1,j,1+(k-1)*5] = round(s[j]*N[round(year[T]-year_0+0.5),(k+1)]) #S
          X[iter,1,j,2+(k-1)*5] = round(e[j]*N[round(year[T]-year_0+0.5),(k+1)]) #E
          X[iter,1,j,3+(k-1)*5] = round(i[j]*N[round(year[T]-year_0+0.5),(k+1)]) #I
          X[iter,1,j,4+(k-1)*5] = 0   #delta N_ir(Biweekly)
          X[iter,1,j,5+(k-1)*5] = N[round(year[T]-year_0+0.5),(k+1)] #P 
        }
      }
      #density for parameters is Norm with sd 
      theta.1[iter,m+1,] = sapply(1:J,function(j) theta.1[iter,m+1,j] = theta.1[iter,m,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) 
                                                                                                          * ifelse(iter > 50,0.0008,0.003))))
      theta.2[iter,m+1,] = sapply(1:J,function(j) theta.2[iter,m+1,j] = min(theta.2[iter,m,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) 
                                                                                                              * ifelse(iter>50,0.0005,0.002))),-theta.1[iter,m+1,j]))
      theta.3[iter,m+1,] = sapply(1:J,function(j) theta.3[iter,m+1,j] = min(theta.3[iter,m,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) 
                                                                                                              * ifelse(iter>50,0.0005,0.002))),
                                                                        (log(1-exp(theta.1[iter,m+1,j]+theta.2[iter,m+1,j])) - 
                                                                           log(exp(theta.1[iter,m+1,j])+exp(theta.2[iter,m+1,j]) 
                                                                               + 2*exp(theta.1[iter,m+1,j]+theta.2[iter,m+1,j])))))
      
      #Proposed X using proposed parameters
      X[iter,m+1,j,] = seir(iter,t,R0,a,alpha,mu,v_ei,v_ir,sigma2,rho,phi,G,X)
      w.o[iter,m,j] = prod(sapply(1:3,function(k) dmeas(proposed_y[m,k],rho,phi, X[iter,m+1,j,4+(k-1)*5])))
      }
    #calculate true weight
    for (j in 1:J) {
      w[iter,m,j] = ifelse(w.o[iter,m,j]==0,0,w.o[iter,m,j]/sum(w.o[iter,m,]))
      #Filtering failure: at some time point, the conditional likelihood of 
      if(w[iter,m,j] > tolerance){
        w[iter,m,j] = w[iter,m,j]
      }else{
        w[iter,m,j] = tolerance
      }}
    #Resample
    number = as.numeric(rmultinom(1,J,w[iter,m,]))
    number[0] = 0
    for(ii in 1:J){
      if((number[ii]==0)==FALSE){
        for(iii in 1:number[ii]){
          n = sum(number[0:(ii-1)])
          for (k in 1:3) {
            X[iter,m+1,iii+n,1+(k-1)*5]=X[iter,m+1,ii,1+(k-1)*5]
            X[iter,m+1,iii+n,2+(k-1)*5]=X[iter,m+1,ii,2+(k-1)*5]
            X[iter,m+1,iii+n,3+(k-1)*5]=X[iter,m+1,ii,3+(k-1)*5]
            X[iter,m+1,iii+n,4+(k-1)*5]=X[iter,m+1,ii,4+(k-1)*5]
            X[iter,m+1,iii+n,5+(k-1)*5]=X[iter,m+1,ii,5+(k-1)*5]
          }
          theta.1[iter,m+1,iii+n]=theta.1[iter,m+1,ii]
          theta.2[iter,m+1,iii+n]=theta.2[iter,m+1,ii]
          theta.3[iter,m+1,iii+n]=theta.3[iter,m+1,ii]
        }}}
  }
  log_likelihood[iter] = sum(log(rowMeans(w.o[iter,,])))
  theta.1[iter+1,1,] = sapply(1:J,function(j) theta.1[iter+1,1,j] = theta.1[iter,M+1,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) * ifelse(iter > 50,0.001,0.005))))
  theta.2[iter+1,1,] = sapply(1:J,function(j) theta.2[iter+1,1,j] = min(theta.2[iter,M+1,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) * ifelse(iter > 50,0.0001,0.001))),-theta.1[iter,M+1,j]))
  theta.3[iter+1,1,] = sapply(1:J,function(j) theta.3[iter+1,1,j] = min(theta.3[iter,M+1,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) * ifelse(iter > 50,0.0001,0.001))),
                                                                    (log(1-exp(theta.1[iter,M+1,j]+theta.2[iter,M+1,j])) - 
                                                                       log(exp(theta.1[iter,M+1,j])+exp(theta.2[iter,M+1,j]) + 2*exp(theta.1[iter,M+1,j]+theta.2[iter,M+1,j])))))
  proposed_10= matrix(0,I+1,J)
  proposed_11 = matrix(0,I+1,J)
  proposed_12 = matrix(0,I+1,J)
  
  proposed_10[iter,] = sapply(1:J,function(j) proposed_10[iter,j] = theta.1[iter,M+1,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) * ifelse(iter > 50,0.001,0.005))))
  proposed_11[iter,] = sapply(1:J,function(j) proposed_11[iter,j] = min(theta.2[iter,M+1,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) * ifelse(iter > 50,0.0001,0.001))),-theta.1[iter,M+1,j]))
  proposed_12[iter,] = sapply(1:J,function(j) proposed_12[iter,j] = min(theta.3[iter,M+1,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) * ifelse(iter > 50,0.0001,0.001))),
                                                                  (log(1-exp(theta.1[iter,M+1,j]+theta.2[iter,M+1,j])) 
                                                                   -  log(exp(theta.1[iter,M+1,j])+exp(theta.2[iter,M+1,j]) + 2*exp(theta.1[iter,M+1,j]+theta.2[iter,M+1,j])))))
  }
                                                                               
for (iter in 2:I) {
  s = sapply(1:J,function(j) logit(theta.1[iter,1,j]))
  e = sapply(1:J,function(j) logit(theta.2[iter,1,j]))
  i = sapply(1:J,function(j) logit(theta.3[iter,1,j]))
  # iterations
  for(m in 1:M){ 
    #Propsed particles
    for (j in 1:J) {
      t = T + m -1
      
      if(m==1){
        # initial X
        for (k in 1:3) {
          X[iter,1,j,1+(k-1)*5] = round(s[j]*N[round(year[T]-year_0+0.5),(k+1)]) #S
          X[iter,1,j,2+(k-1)*5] = round(e[j]*N[round(year[T]-year_0+0.5),(k+1)]) #E
          X[iter,1,j,3+(k-1)*5] = round(i[j]*N[round(year[T]-year_0+0.5),(k+1)]) #I
          X[iter,1,j,4+(k-1)*5] = 0   #delta N_ir(Biweekly)
          X[iter,1,j,5+(k-1)*5] = N[round(year[T]-year_0+0.5),(k+1)] #P 
        }
      }
      
      #density for parameters is Norm with sd based on the geometric cooling scheme
      theta.1[iter,m+1,] = sapply(1:J,function(j) theta.1[iter,m+1,j] = theta.1[iter,m,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) 
                                                                                                             * ifelse(iter > 50,0.0008,0.003))))
      theta.2[iter,m+1,] = sapply(1:J,function(j) theta.2[iter,m+1,j] = min(theta.2[iter,m,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) 
                                                                                                                 * ifelse(iter>50,0.0005,0.002))),-theta.1[iter,m+1,j]))
      theta.3[iter,m+1,] = sapply(1:J,function(j) theta.3[iter,m+1,j] = min(theta.3[iter,m,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) 
                                                                                                                 * ifelse(iter>50,0.0005,0.002))),
                                                                          (log(1-exp(theta.1[iter,m+1,j]+theta.2[iter,m+1,j])) - 
                                                                             log(exp(theta.1[iter,m+1,j])+exp(theta.2[iter,m+1,j]) 
                                                                                 + 2*exp(theta.1[iter,m+1,j]+theta.2[iter,m+1,j])))))
      
      #Proposed X using proposed parameters
      X[iter,m+1,j,] = seir(iter,t,R0,a,alpha,mu,v_ei,v_ir,sigma2,rho,phi,G,X)
      w.o[iter,m,j] = prod(sapply(1:3,function(k) dmeas(proposed_y[m,k],rho,phi, X[iter,m+1,j,4+(k-1)*5])))
      }
    #calculate true weight
    for (j in 1:J) {
      w[iter,m,j] = ifelse(w.o[iter,m,j]==0,0,w.o[iter,m,j]/sum(w.o[iter,m,]))
      #Filtering failure: at some time point, the conditional likelihood of 
      if(w[iter,m,j] > tolerance){
        w[iter,m,j] = w[iter,m,j]
      }else{
        w[iter,m,j] = tolerance
      }}
    #Resample
    number = as.numeric(rmultinom(1,J,w[iter,m,]))
    number[0] = 0
    for(ii in 1:J){
      if((number[ii]==0)==FALSE){
        for(iii in 1:number[ii]){
          n = sum(number[0:(ii-1)])
          for (k in 1:3) {
            X[iter,m+1,iii+n,1+(k-1)*5]=X[iter,m+1,ii,1+(k-1)*5]
            X[iter,m+1,iii+n,2+(k-1)*5]=X[iter,m+1,ii,2+(k-1)*5]
            X[iter,m+1,iii+n,3+(k-1)*5]=X[iter,m+1,ii,3+(k-1)*5]
            X[iter,m+1,iii+n,4+(k-1)*5]=X[iter,m+1,ii,4+(k-1)*5]
            X[iter,m+1,iii+n,5+(k-1)*5]=X[iter,m+1,ii,5+(k-1)*5]
          }
          theta.1[iter,m+1,iii+n]=theta.1[iter,m+1,ii]
          theta.2[iter,m+1,iii+n]=theta.2[iter,m+1,ii]
          theta.3[iter,m+1,iii+n]=theta.3[iter,m+1,ii]
        }}}
  }
  log_likelihood[iter] = sum(log(rowMeans(w.o[iter,,])))
  
  proposed_10[iter,] = sapply(1:J,function(j) proposed_10[iter,j] = theta.1[iter,M+1,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) * ifelse(iter > 50,0.001,0.005))))
  proposed_11[iter,] = sapply(1:J,function(j) proposed_11[iter,j] = min(theta.2[iter,M+1,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) * ifelse(iter > 50,0.0001,0.001))),-theta.1[iter,M+1,j]))
  proposed_12[iter,] = sapply(1:J,function(j) proposed_12[iter,j] = min(theta.3[iter,M+1,j] + rnorm(1,0,sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) * ifelse(iter > 50,0.0001,0.001))),
                                                                      (log(1-exp(theta.1[iter,M+1,j]+theta.2[iter,M+1,j])) 
                                                                       -  log(exp(theta.1[iter,M+1,j])+exp(theta.2[iter,M+1,j]) + 2*exp(theta.1[iter,M+1,j]+theta.2[iter,M+1,j])))))
  likelihoodrate <- (log_likelihood[iter]/log_likelihood[iter-1])
  acceptrate_10 <- likelihoodrate*(dnorm(proposed_10[iter,],-3,1)/dnorm(proposed_10[iter-1,],-3,1))*(dnorm(proposed_10[iter-1,],proposed_10[iter,],sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) * ifelse(iter > 50,0.001,0.005))) /dnorm(proposed_10[iter,],proposed_10[iter-1,],sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) 
                                                                     * ifelse(iter > 50,0.001,0.005))))
  acceptrate_11 <- likelihoodrate*(dnorm(proposed_11[iter,],-8,1)/dnorm(proposed_11[iter-1,],-8,1))*(dnorm(proposed_11[iter-1,],proposed_11[iter,],sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) 
                                                      * ifelse(iter > 50,0.001,0.005))) 
    /dnorm(proposed_11[iter,],proposed_11[iter-1,],sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) 
                                                       * ifelse(iter > 50,0.001,0.005))))
  acceptrate_12 <- likelihoodrate*(dnorm(proposed_12[iter,],-8,1)/dnorm(proposed_12[iter-1,],-8,1)) *(dnorm(proposed_12[iter-1,],proposed_12[iter,],sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) 
                                                      * ifelse(iter > 50,0.001,0.005))) 
    /dnorm(proposed_12[iter,],proposed_12[iter-1,],sd = sqrt(0.92^(((m-1)+(iter-1)*M)/(50*M)) 
                                                       * ifelse(iter > 50,0.001,0.005))))
  for (j in 1:J){
   if (runif(1)<acceptrate_10[j]){  
     theta.1[iter+1,1,j] = proposed_10[iter,j]
   }else{
     theta.1[iter+1,1,j] = theta.1[iter,1,j]
    }
   
  if (runif(1)<acceptrate_11[j]){
    theta.2[iter+1,1,j] = proposed_11[iter,j]
  }else{
    theta.2[iter+1,1,j] = theta.2[iter,1,j]
  }
  if (runif(1)<acceptrate_12[j]){
    theta.3[iter+1,1,j] = proposed_12[iter,j]
  }else{
    theta.3[iter+1,1,j] = theta.3[iter,1,j]
  }
  }
}
  
  
  






















####plot####

par(mfrow=c(1,2))
plot(rowMeans(logit(theta.1[,1,])),ylim=c(min(logit(theta.1)),0.1),type='l',col=1,main = 'Initial proportion of Susceptible',xlab = 'Iteraion times',ylab = 'I(S)')
lines(sapply(1:I,function(iter) median(logit(theta.1[iter,1,]),0.1)),col=2,lty=1)
lines(sapply(1:I,function(iter) quantile(logit(theta.1[iter,1,]),0.1)),col=3,lty=2)
lines(sapply(1:I,function(iter) quantile(logit(theta.1[iter,1,]),0.9)),col=3,lty=2)
lines(sapply(1:I,function(iter) 0.04),col=4,lty=2)
legend("topright", bty="n",legend=c("Est. I(S):mean","Est. median",'Est. 10%,90% quantile',"True Value"),col=c(1,2,3,4), lty=c(1,1,2,2), cex=0.7)

plot(rowMeans(logit(theta.2[,1,])),ylim=c(min(logit(theta.2)),0.001),type='l',col=1,main = 'Initial proportion of Exposed',xlab = 'Iteraion times',ylab = 'I(E)')
lines(sapply(1:I,function(iter) median(logit(theta.2[iter,1,]),0.1)),col=2,lty=1)
lines(sapply(1:I,function(iter) quantile(logit(theta.2[iter,1,]),0.1)),col=3,lty=2)
lines(sapply(1:I,function(iter) quantile(logit(theta.2[iter,1,]),0.9)),col=3,lty=2)
lines(sapply(1:I,function(iter) 0.00027),col=4,lty=2)
legend("topright", bty ="n",legend=c("Est. I(S):mean","Est. median",'Est. 10%,90% quantile',"True Value"),col=c(1,2,3,4), lty=c(1,1,2,2), cex=0.7)

plot(rowMeans(logit(theta.3[,1,])),ylim=c(min(logit(theta.3)),0.0015),type='l',col=1,main = 'Initial proportion of Infectious',xlab = 'Iteraion times',ylab = 'I(I)')
lines(sapply(1:I,function(iter) median(logit(theta.3[iter,1,]),0.1)),col=2,lty=1)
lines(sapply(1:I,function(iter) quantile(logit(theta.3[iter,1,]),0.1)),col=3,lty=2)
lines(sapply(1:I,function(iter) quantile(logit(theta.3[iter,1,]),0.9)),col=3,lty=2)
lines(sapply(1:I,function(iter) 0.00032),col=4,lty=2)
legend("topright", bty ="n",legend=c("Est. I(S):mean","Est. median",'Est. 10%,90% quantile',"True Value"),col=c(1,2,3,4), lty=c(1,1,2,2), cex=0.7)

