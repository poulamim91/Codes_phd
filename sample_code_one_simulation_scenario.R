rm(list=ls())
setwd("/nas/longleaf/home/poulamim/gen_recur6/") #### Working Dir ####
require(survival)
runcase<-function(ln){
  setwd("/nas/longleaf/home/poulamim/gen_recur6/") #### Working Dir ####
  #obs_size = c(1000,1000,1000,2000,2000,2000,4000,4000,4000)
  #g_list = c(0.25,0.4,0.65,0.25,0.4,0.65,0.25,0.4,0.65)
  #cenn = c(1.33,1.33,1.9,1.33,1.33,1.9,1.33,1.33,1.9)
  obs_size = c(1000,2000,4000,6000)
  g_list = rep(0.2,4)
  cenn = rep(4.5, 4)
  ob_num=obs_size[ln]
  g_val = g_list[ln]
  c11 = cenn[ln]
  cp = rep(60,4)
  log_file <- file(paste0("gen_",ob_num,"n",cp[ln],"recur6.log"), open ="w") # Change log file location #
  sink(log_file, append=TRUE)
  sink(log_file, append=TRUE, type="message")
  data_gen_new_fixed <- function(n=1000,r=20,beta_1=0.5,g1=0.05,g2=0.05,cens)
  {
    z <- runif(n,0,1)
    R <- matrix(NA,nrow=n,ncol=r)
    v = rgamma(n,g1,g2)
    amat<- matrix(runif(n*r,0,1),ncol=r,nrow=n)
    R[,1]=-log(amat[,1])/v/exp(beta_1*z)
    for(i in 2:20)
      R[,i] <- R[,i-1] -log(amat[,i])/v/exp(beta_1*z)
    flag_2=0
    C <- runif(n,0,cens)
    del_ <- mean(R[,1]<= C)
    t1<- R[which(R[,1]<=C),]
    C1 <- C[which(R[,1]<=C)]
    cnt=0
    for(k in 1:nrow(t1))
      cnt = cnt+sum(t1[k,]<=C1[k])
    mean_rec_ <- cnt/nrow(t1)
    return(list(c(mean_rec_,del_,mean(z),cens),R,C,z))
 }
   create_r_u <- function(R_obsd,Cens)
   {
     temp <- rep(0,1)
    for(i in 1:length(Cens))
   {
       t = which(R_obsd[i,]<Cens[i])
       if(length(t)>0) temp=c(temp,R_obsd[i,t])
     }
     temp=temp[-1]
     return(sort(unique(temp)))
   }
  create_y <- function(R_u,Cens)
   {
     temp <- matrix(NA,nrow=length(R_u),ncol=length(Cens))
     for( i in 1:length(Cens))
       temp[,i] = as.numeric(R_u<=Cens[i])
     return(temp)
   }
   create_alpha <- function(i,zi,Y_mat,Delta_)
   {
     temp=sum((1-Delta_)*Y_mat[i,])
     out=ifelse(temp==0,0,sum((1-Delta_)*zi*Y_mat[i,])/temp)
     out
   }
   create_q <- function(i,zi,eta,Y_mat,Delta_)
   {
    temp=sum(Delta_*(1-zi)*Y_mat[i,])
    out=ifelse(temp==0,0,sum(Delta_*(1-zi)*eta*Y_mat[i,])/temp)
    out
   }
 create_wt1 <- function(index,alpha_hat,q_hat,Delta_,zi,eta)
   {
     if(alpha_hat[index]==0){t1=0}else{t1=(1-Delta_)*zi/alpha_hat[index]}
     if(q_hat[index]==0){t2=0}else{t2=Delta_*(1-zi)*eta/q_hat[index]}
     return(Delta_*zi+t1+t2)
  }
  create_wt <- function(R_u,Delta_,Y_mat,zi,eta)
    {
      alpha_hat <- sapply(1:length(R_u),create_alpha,zi,Y_mat,Delta_)
      q_hat <- sapply(1:length(R_u),create_q,zi,eta,Y_mat,Delta_)
      wt<-matrix(unlist(lapply(1:length(R_u),create_wt1,alpha_hat,q_hat,Delta_,zi,eta)),nrow=length(R_u),ncol=length(Delta_),byrow=T)
      return(wt)
    }
  create_u1 <- function(i,R,R_u,Cens,Y_mat,w_mat,z_obsd,beta0,s1,s0inv)
   {
    times <- which(R[i,]<Cens[i])
    if(length(times)==0)
     {output=0}else
     {
       output=0
       for(k in 1:length(times))
       {
         r1=which(R_u==R[i,times[k]])
         output= output+w_mat[r1,i]*(z_obsd[i] - s1[r1]*s0inv[r1])
       }
     }
     return(output)
   }
  create_u <- function(beta0,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv)
  {
  out <- sapply(1:length(Cens),create_u1,R,R_u,Cens,Y_mat,w_mat,z_obsd,beta0,s1,s0inv)
  sum(out)
  }
  create_uprime <- function(beta0,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv,s2)
  {
  out <- sapply(1:length(Cens),create_u1prime,R,R_u,Cens,Y_mat,w_mat,z_obsd,beta0,s1,s0inv,s2)
  sum(out)
  }
 create_u1prime <- function(i,R,R_u,Cens,Y_mat,w_mat,z_obsd,beta0,s1,s0inv,s2)
 {
  times <- which(R[i,]<Cens[i])
  if(length(times)==0)
  {output=0}else
  {
   output=0
   for(k in 1:length(times))
    {
     r1=which(R_u==R[i,times[k]])
     output= output - w_mat[r1,i]*(s2[r1]*s0inv[r1]-(s1[r1]*s0inv[r1])^2)
     }
   }
  return(output)
  }
 create_M <- function(i,betaest,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv)
  {
   N=nrow(R)
   times <- which(R[i,]<Cens[i])
   if(length(times)==0)
   {output=0}else
   {
    output=0
    for(k in 1:length(times))
    {
    r1=which(R_u==R[i,times[k]])
    output= output+z_obsd[i] - s1[r1]*s0inv[r1]
    }
  }
  for(j in 1:nrow(R))
  {
  times <- which(R[j,]<Cens[j])
  if(length(times)>0)
  {
  for(k in 1:length(times))
  {
  r1=which(R_u==R[j,times[k]])
  output= output-(z_obsd[i]-s1[r1]*s0inv[r1])*exp(betaest*z_obsd[i])*Y_mat[r1,i]*w_mat[r1,j]*s0inv[r1]/N
  }
  }
  }
  return(output)
}
create_s1 <- function(i,w_mat,Y_mat,z_obsd,betaest)
{
return(mean(w_mat[i,]*Y_mat[i,]*z_obsd*exp(betaest*z_obsd)))
}
create_s2 <- function(i,w_mat,Y_mat,z_obsd,betaest)
  {
  return(mean(w_mat[i,]*Y_mat[i,]*z_obsd^2*exp(betaest*z_obsd)))
  }
  create_s0inv <- function(i,w_mat,Y_mat,z_obsd,betaest)
  {
  temp=mean(w_mat[i,]*Y_mat[i,]*exp(betaest*z_obsd))
  temp=ifelse(temp==0,0,1/temp)
  return(temp)
  }
create_v1 <- function(j,index,R,Cens,R_u,Y_mat,Ri_bt,s0inv,E_R,E_Y)
  {
  times <- which(R[j,]<Cens[j])
  out=0
  if(length(times)>0)
  {
    for(k in 1:length(times))
    {
      r1=which(R_u==R[j,times[k]])
      out=out+ sum((Ri_bt[r1,index]-Y_mat[r1,index]*E_R[r1]/E_Y[r1])*s0inv[r1])
    }
  }
  return(out)
}
create_v <- function(index,zi,Delta,R,Cens,R_u,Y_mat,Ri_bt,s0inv,E_R,E_Y)
  {
  if(zi[index]==0 | Delta[index]==1)
    {output=0}else{
    p1<-unlist(lapply(1:length(Delta),create_v1,index,R,Cens,R_u,Y_mat,Ri_bt,s0inv,E_R,E_Y))
    output=(mean(p1))^2
    }
  return(output)
  }
create_v2_1 <- function(j,index,R,Cens,R_u,Y_mat,Mhat,E_M,E_Y2)
  {
    times <- which(R[j,]<Cens[j])
    out=0
    if(length(times)>0)
    {
    for(k in 1:length(times))
    {
    r1=which(R_u==R[j,times[k]])
    out=out+ sum(Y_mat[r1,index]*E_M[r1]/E_Y2[r1])
    }
    }
    return(out)
  }
create_v2<-function(index,zi,eta,Delta,R,Cens,R_u,Y_mat,Mhat,E_M,E_Y2)
   {
    if(eta[index]==0)
     {output=0}else{
      p1<- unlist(lapply(1:length(Delta),create_v2_1,index,R,Cens,R_u,Y_mat,Mhat,E_M,E_Y2))
      output<- (Mhat[index]-sum(p1))^2
     }
     output
  }
 create_Ribt1<-function(index,z_obsd,Y_mat,s1,s0inv,betaest)
  {
  return((z_obsd[index]-s1*s0inv)*Y_mat[,index]*exp(betaest*z_obsd[index]))
  }
 create_Ribt<-function(z_obsd,Y_mat,s1,s0inv,betaest)
  {
   out<- unlist(lapply(1:length(z_obsd),create_Ribt1,z_obsd,Y_mat,s1,s0inv,betaest))
   return(matrix(out,nrow=nrow(Y_mat),ncol=ncol(Y_mat),byrow = F))
  }
 create_v2_1_new <- function(j,index,R,Cens,R_u,Y_mat,E_M,E_Y2)
  {
   times <- which(R[j,]<Cens[j])
   out=0
   if(length(times)>0)
   {
   for(k in 1:length(times))
   {
   r1=which(R_u==R[j,times[k]])
   out=out+ sum(Y_mat[r1,index]*E_M[r1]/E_Y2[r1])
   }
   }
   return(out)
 }
 create_v2_new<-function(index,zi,Delta,R,Cens,R_u,Y_mat,Mhat,E_M,E_Y2,t2)
 {
  p1<- unlist(lapply(1:length(Delta),create_v2_1_new,t2[index],R,Cens,R_u,Y_mat,E_M,E_Y2))
  output<- (Mhat[index]-sum(p1))^2
  output
 }
 create_dM1 <- function(i,betaest,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv,Delta){
   N=nrow(R)
   out=rep(0,length(R_u))
   times <- which(R[i,]<Cens[i])
   if(length(times)>0)
   {
    for(k in 1:length(times)){
     r1=which(R_u==R[i,times[k]])
     out[r1]= out[r1]+z_obsd[i] - s1[r1]*s0inv[r1]
     }
   }
    for(j in 1:N)
    {
    times <- which(R[j,]<Cens[j])
    if(length(times)>0)
    {
    for(k in 1:length(times)){
    r1=which(R_u==R[j,times[k]])
    out[r1]= out[r1]-(z_obsd[i]-s1[r1]*s0inv[r1])*exp(betaest*z_obsd[i])*Y_mat[r1,i]*w_mat[r1,j]*s0inv[r1]/N
    }
    }
  }
 return(out)
}
create_dM<-function(betaest,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv,Delta)
 {
  out<- unlist(lapply(1:length(z_obsd),create_dM1,betaest,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv,Delta))
  return(matrix(out,nrow=nrow(Y_mat),ncol=ncol(Y_mat),byrow = F))
}
create_dM_new<-function(betaest,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv,Delta,t1)
{
out<- unlist(lapply(t1,create_dM1,betaest,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv,Delta))
return(matrix(out,nrow=nrow(Y_mat),ncol=length(t1),byrow = F))
}
exp_R <-function(index,alphatilde,zi,Delta,Ri_bt)
{
  mean(zi/alphatilde*(1-Delta)*Ri_bt[index,])
}
exp_Y<-function(index,Delta,Y_mat)
  {
    mean((1-Delta)*Y_mat[index,])
  }
exp_M <-function(index,qtilde,zi,eta,Delta,dM)
  {
  sum(eta/qtilde*dM[index,])/sum(Delta*(1-zi))
  }
exp_M_new <-function(index,qtilde,zi,Delta,dM)
{
    sum(dM[index,]/qtilde)/sum(Delta*(1-zi))
}       
exp_Y2<-function(index,Delta,Y_mat)
{
  sum(Delta*Y_mat[index,])/sum(Delta)
}
beta_estimation <- function(R,R_u,Cens,Y_mat,z_obsd,Delta,wt_mat,initial_beta=0,tol=0.01,max.iter=50)
{
beta0<- initial_beta
flag=0
m=0
while(flag==0 & m<max.iter)
{
  s1 <- sapply(1:length(R_u),create_s1,w_mat,Y_mat,z_obsd,beta0)
  s1 <- ifelse(is.na(s1),0,s1)
  s0inv <- sapply(1:length(R_u),create_s0inv,w_mat,Y_mat,z_obsd,beta0)
  s0inv <- ifelse(is.na(s0inv),0,s0inv)
  s2 <- sapply(1:length(R_u),create_s2,w_mat,Y_mat,z_obsd,beta0)
  s2 <- ifelse(is.na(s2),0,s2)
  u_beta <- create_u(beta0,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv)
  uprime_beta <- create_uprime(beta0,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv,s2)
  beta1 <- beta0 - u_beta/uprime_beta
  if(abs(beta1-beta0)<tol) flag=1
  beta0=beta1
  m=m+1
  }
  beta0=ifelse(m>=max.iter,NA,beta0)
  return(beta0)
}
var_estimation <- function(R,R_u,Cens,Y_mat,z_obsd,Delta,wt_mat,zi,eta,betaest)
{
N=nrow(R)
alphatilde<- mean(zi)
qtilde <- mean(eta)/mean(Delta*(1-zi))
s1 <- sapply(1:length(R_u),create_s1,w_mat,Y_mat,z_obsd,betaest)
s1 <- ifelse(is.na(s1),0,s1)
s0inv <- sapply(1:length(R_u),create_s0inv,w_mat,Y_mat,z_obsd,betaest)
s0inv <- ifelse(is.na(s0inv),0,s0inv)
s2 <- sapply(1:length(R_u),create_s2,w_mat,Y_mat,z_obsd,betaest)
s2 <- ifelse(is.na(s2),0,s2)
Ahat <- -create_uprime(betaest,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv,s2)/nrow(R)
t1<- which(zi==1)
d1M<-create_dM_new(betaest,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv,Delta,t1)
Mhat <- apply(d1M,2,sum)
Qhat <- sum(Mhat^2)/alphatilde/N
Ri_bt<-create_Ribt(z_obsd,Y_mat,s1,s0inv,betaest)
E_R <-unlist(lapply(1:nrow(Ri_bt),exp_R,alphatilde,zi,Delta,Ri_bt))
E_Y <-unlist(lapply(1:nrow(Ri_bt),exp_Y,Delta,Y_mat))
Vhat <- mean(unlist(lapply(1:ncol(Y_mat),create_v,zi,Delta,R,Cens,R_u,Y_mat,Ri_bt,s0inv,E_R,E_Y)))/alphatilde
t2<- which(eta==1)
d2M<-create_dM_new(betaest,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv,Delta,t2)
M2hat <- apply(d2M,2,sum)
E2_M<-unlist(lapply(1:nrow(Ri_bt),exp_M_new,qtilde,zi,Delta,d2M))
E_Y2 <- unlist(lapply(1:nrow(Ri_bt),exp_Y2,Delta,Y_mat))
V2hat <- sum(unlist(lapply(1:length(t2),create_v2_new,zi,Delta,R,Cens,R_u,Y_mat,M2hat,E2_M,E_Y2,t2)))/sum(Delta*(1-zi))/qtilde
output <- Ahat^(-2)*(Qhat+(1-alphatilde)/alphatilde*Vhat+(1-alphatilde)*mean(Delta)*(1-qtilde)/qtilde*V2hat)/N
return(output)
}
set.seed(555)
sim.num=1000
beta_ests <- rep(NA,sim.num)
beta_ests_cox <- rep(NA,sim.num)
var_ests <- rep(NA,sim.num)
for( loop.num in 1:sim.num)
{
  t1 <- proc.time()
  coh.size = ob_num
  dat_list <-data_gen_new_fixed(n = coh.size,beta_1=0.5,g1=g_val,g2=g_val,cens=c11)
  avg_atleast1_recur = dat_list[[1]][2]
  avg_recur = dat_list[[1]][1]
  sigma = dat_list[[1]][3]
  R_obsd = dat_list[[2]]
  Cens = dat_list[[3]]
  z_obsd = dat_list[[4]]
  R_u <-create_r_u(R_obsd,Cens)
  Y_mat <- create_y(R_u,Cens)
  Delta <- as.numeric(R_obsd[,1]<Cens)
  zi <- rbinom(ncol(Y_mat),1,0.075)
  eta<- rep(0,ncol(Y_mat))
  sub<- which(Delta*(1-zi)==1)
  eta[sub]<-rbinom(length(sub),1,0.1)
  w_mat <- create_wt(R_u,Delta,Y_mat,zi,eta)
  temp_out <- try(beta_estimation(R_obsd,R_u,Cens,Y_mat,z_obsd,Delta,w_mat))
  if("try-error" %in% class(temp_out) | is.na(temp_out))
    {
      beta_ests[loop.num]=NA
      var_ests[loop.num]=NA
    }else{
      beta_ests[loop.num]=temp_out
      var_ests[loop.num] <- var_estimation(R_obsd,R_u,Cens,Y_mat,z_obsd,Delta,w_mat,zi,eta,beta_ests[loop.num])
    }
  t2 <- proc.time()
  message("Total time taken for ",loop.num)
  print(t2-t1)
}
obs_points<- which(!is.na(beta_ests))
print(beta_ests)
print(var_ests)
print(mean(beta_ests[obs_points]))
print(var(beta_ests[obs_points]))
print(mean(var_ests[obs_points]))
print(mean(abs(beta_ests[obs_points])<1.96*sqrt(var_ests[obs_points])))
sink()     
}

