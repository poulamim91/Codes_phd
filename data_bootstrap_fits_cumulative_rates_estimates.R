rm(list=ls())

require(survival)
library(MASS)
library(ggplot2)
library(gridExtra)
library(grid)

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
    sum((1-Delta_)*zi*Y_mat[i,])/sum((1-Delta_)*Y_mat[i,])
  }
  
  
  create_wt <- function(R_u,Delta_,Y_mat,zi)
  {
    
    
    alpha_hat <- sapply(1:length(R_u),create_alpha,zi,Y_mat,Delta_)
    
    wt<- matrix(0,nrow=nrow(Y_mat),ncol=ncol(Y_mat))
    for (i in 1:length(R_u))
    {
      for(j in 1:length(Delta_))
      {
        t1<- ifelse(alpha_hat[i]==0,0,(1-Delta_[j])*zi[j]/alpha_hat[i])
        wt[i,j]=Delta_[j]+t1
      }
      
    }
    
    return(wt)
  }

  
  create_s1 <- function(i,w_mat,Y_mat,z_obsd,betaest)
  {
    return(apply(as.vector(w_mat[i,]*Y_mat[i,]*exp(z_obsd%*%betaest))*z_obsd, 2, mean))
    
  }
    
  
  create_s2 <- function(i,w_mat,Y_mat,z_obsd,betaest)
  {
    n = nrow(z_obsd)
    z_mat = matrix(0, ncol = ncol(z_obsd), nrow = ncol(z_obsd))
    for(j in 1:n){
      c = (w_mat[i,j]*Y_mat[i,j]*exp(z_obsd[j,]%*%betaest))
      z_mat = z_mat + as.vector(c)*(z_obsd[j,]%*%t(z_obsd[j,]))
    }
    z_mat = z_mat/n
    return(list(z_mat))
  }
  
  create_s0inv <- function(i,w_mat,Y_mat,z_obsd,betaest)
  {
    temp=mean(w_mat[i,]*Y_mat[i,]*exp(z_obsd%*%betaest))
    temp=ifelse(temp==0,0,1/temp)
    return(temp)
    
  }
  
  create_u1 <- function(i,R,R_u,Cens,Y_mat,w_mat,z_obsd,beta0,s1,s0inv)
  {
    times <- which(R[i,]<Cens[i])
    output = rep(0, ncol(z_obsd))
    if(length(times)>0)
    { 
      for(k in 1:length(times))
      {
        r1=which(R_u==R[i,times[k]])
        output= output+ (z_obsd[i,] - unlist(s1[[r1]])*s0inv[r1])
      }
    }
    return(output) 
    
  }
  
  create_u <- function(beta0,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv)
  {
    out <- lapply(1:length(Cens),create_u1,R,R_u,Cens,Y_mat,w_mat,z_obsd,beta0,s1,s0inv)
    out1 = Reduce('+', out)
   return(out1) 
  }
  
  create_uprime <- function(beta0,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv,s2)
  {
    out <- lapply(1:length(Cens),create_u1prime,R,R_u,Cens,Y_mat,w_mat,z_obsd,beta0,s1,s0inv,s2)
    out1 = Reduce('+', out) 
  }
  
  create_u1prime <- function(i,R,R_u,Cens,Y_mat,w_mat,z_obsd,beta0,s1,s0inv,s2)
  {
    times <- which(R[i,]<Cens[i])
    output = matrix(0, nrow = ncol(z_obsd), ncol = ncol(z_obsd))
    if(length(times)>0)
    {
      for(k in 1:length(times))
      {
        r1=which(R_u==R[i,times[k]])
        output= output - unlist(s2[[r1]])*s0inv[r1] + unlist(s1[[r1]])%*%t(unlist(s1[[r1]]))*s0inv[r1]^2
      }
    }
    return(output) 
    
  }
  
  create_M <- function(i,betaest,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv)
  {
    N=nrow(R)
    times <- which(R[i,]<Cens[i])
    output=rep(0, ncol(z_obsd))
    if(length(times)>0)
    {
      for(k in 1:length(times))
      {
        r1=which(R_u==R[i,times[k]])
        output= output+z_obsd[i,] - unlist(s1[[r1]])*s0inv[r1]
      }
    }
    
    for(j in 1:N)
    {
      times <- which(R[j,]<Cens[j])
      if(length(times)>0)
      {
        for(k in 1:length(times))
        {
          r1=which(R_u==R[j,times[k]])
          output= output-(z_obsd[i,]-unlist(s1[[r1]])*s0inv[r1])*as.vector(exp(z_obsd[i,]%*%betaest)*Y_mat[r1,i]*s0inv[r1]/N)
        }
      }
      
    }
    return(output%*%t(output))
    
  }

  create_Ribt1<-function(index,z_obsd,Y_mat,s1,s0inv,betaest)
  {
    r = matrix(0, ncol = length(s0inv), nrow = ncol(z_obsd))
    for(i in 1:ncol(z_obsd)){
      r[i,] =(z_obsd[index,i]-as.vector(matrix(unlist(s1), ncol = ncol(z_obsd), byrow = T)[,i]*s0inv))*as.vector(Y_mat[,index]*exp(betaest[i]*z_obsd[index,i]))
    }
    return(r)
  }
  
  create_Ribt<-function(z_obsd,Y_mat,s1,s0inv,betaest)
  {
    
    out<- lapply(1:nrow(z_obsd),create_Ribt1,z_obsd,Y_mat,s1,s0inv,betaest)
    return(out)
  }
  
  
  create_v1 <- function(j,index,R,Cens,R_u,Y_mat,Ri_bt,s0inv,E_R,E_Y)
  {
    times <- which(R[j,]<Cens[j])
    out=rep(0, 2)
    if(length(times)>0)
    {
      for(k in 1:length(times))
      {
        r1=which(R_u==R[j,times[k]])
        out=out+ (Ri_bt[[index]][,r1]-Y_mat[r1,index]*E_R[[r1]]/E_Y[r1])*s0inv[r1]
      }
    }
    return(out)
  }
  
  create_v <- function(index,zi,Delta,R,Cens,R_u,Y_mat,Ri_bt,s0inv,E_R,E_Y)
  {
    #output = matrix(0, ncol = 4, nrow = 4)
    if(zi[index]==0 | Delta[index]==1)
    {
      output = matrix(0, ncol = 2, nrow = 2)
    }else{
      p1<-lapply(1:length(Delta),create_v1,index,R,Cens,R_u,Y_mat,Ri_bt,s0inv,E_R,E_Y)
      output = Reduce('+', p1)/ncol(Y_mat)
      }
    return(output%*%t(output))
  }
  
  exp_R <-function(index,alphatilde,zi,Delta,Ri_bt)
  {
    r = rep(0, nrow(Ri_bt[[1]]))
    for(i in 1:length(zi)){
      r = r + zi[i]/alphatilde*(1-Delta[i])*unlist(Ri_bt[[i]])[,index]
    }
    r = r/length(zi)
    return(r)
  }
  
  exp_Y<-function(index,Delta,Y_mat)
  {
    mean((1-Delta)*Y_mat[index,])
    
  }
  
 
  beta_estimation <- function(R,R_u,Cens,Y_mat,z_obsd,Delta,wt_mat,initial_beta,tol=0.0001,max.iter=1000)
  {
    beta0<- initial_beta
    flag=0
    m=0
    while(flag==0 & m<max.iter)
    {
      s1 <- lapply(1:length(R_u),create_s1,w_mat,Y_mat,z_obsd,beta0)
      s0inv <- sapply(1:length(R_u),create_s0inv,w_mat,Y_mat,z_obsd,beta0)
      s2 <- lapply(1:length(R_u),create_s2,w_mat,Y_mat,z_obsd,beta0)
      u_beta <- create_u(beta0,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv)
      uprime_beta <- create_uprime(beta0,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv,s2)
      beta1 <- beta0 - solve(uprime_beta)%*%u_beta
      if(sum((beta1-beta0)^2)<tol) flag=1
      beta0=beta1
      m=m+1
     }
    return(beta0) 
  }
  
  
  var_estimation <- function(R,R_u,Cens,Y_mat,z_obsd,Delta,wt_mat,zi,betaest)
  {
    N=nrow(R)
    alphatilde<- mean(zi)
    s1 <- lapply(1:length(R_u),create_s1,w_mat,Y_mat,z_obsd,betaest)
    s0inv <- sapply(1:length(R_u),create_s0inv,w_mat,Y_mat,z_obsd,betaest)
    s0inv <- ifelse(is.na (s0inv),0,s0inv)
    s2 <- lapply(1:length(R_u),create_s2,w_mat,Y_mat,z_obsd,betaest)
    Ahat <- -create_uprime(betaest,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv,s2)/N
    Mhat <- lapply(1:N,create_M,betaest,R,R_u,Cens,Y_mat,w_mat,z_obsd,s1,s0inv)
    Mhat1 <- mapply("*", Mhat, zi, SIMPLIFY = FALSE) 
    Qhat <- Reduce('+', Mhat1)/alphatilde/N
    
    Ri_bt<-create_Ribt(z_obsd,Y_mat,s1,s0inv,betaest)
    E_R <-lapply(1:ncol(Ri_bt[[1]]),exp_R,alphatilde,zi,Delta,Ri_bt)
    E_Y <- rep(0, nrow(Y_mat))
    for(i in 1: nrow(Y_mat)){
      E_Y[i] =  mean((1-Delta)*Y_mat[i,])
    }
    Vhat <- lapply(1:ncol(Y_mat),create_v,zi,Delta,R,Cens,R_u,Y_mat,Ri_bt,s0inv,E_R,E_Y)
    Vhat1 <- Reduce('+', Vhat)/alphatilde/N
    
    output <- solve(Ahat)%*%(Qhat+(1-alphatilde)/alphatilde*Vhat1)%*%solve(Ahat)/N
    return(output)
  }

  
  create_dn1 <- function(i,R_u,R,Cens){
    out = ifelse(R[i,1] == R_u, 1, 0)
    for(j in 2:length(R[i,])){
      out = out +ifelse(R[i,j] == R_u, 1, 0)
    }
    return(out*ifelse(Cens[i] >= R_u, 1, 0))
  }
  
  create_dn <- function(R_u,R,Cens){
    out <- lapply(1:length(Cens),create_dn1,R_u,R,Cens)
    return(matrix(unlist(out),nrow = length(R_u), ncol = length(Cens), byrow = F)) 
  }
  
  create_mu1 <- function(i,w_mat,Y_mat,z_obsd,betaest,dN)
  {
    numer = sum(dN[i,])
    dinom=create_s0inv(i,w_mat,Y_mat,z_obsd,betaest)
    out=ifelse(dinom==0,0,numer*dinom/nrow(z_obsd))
    return(out)
    
  }
  
  create_mu <- function(R_u,w_mat,Y_mat,z_obsd,betaest,dN){
    out <- sapply(1:length(R_u),create_mu1,w_mat,Y_mat,z_obsd,betaest,dN)
    return(out) 
  }
  
  library(tidyverse)
  set.seed(12345)
  ############################################################################################################################################
      t1 <- proc.time()
      library(foreign)
      #Assume we have already read the data set and saved it as brazil.
      brazil1 = brazil[!duplicated(brazil$numcri),]
      n = nrow(brazil1)
    R_u = sort(unique(brazil$tstop[which(brazil$status == 1)]))
    r = max(aggregate(brazil$status, by=list(Category=brazil$numcri), FUN=sum)[,2])
    tau = max(brazil$tstop)
    Cens = rep(0, n)
    R_obsd = matrix(0, nrow = n, ncol = r)
    z_obsd = matrix(0, nrow = n, ncol = 5)
    l = 1
    for(i in unique(brazil$numcri)){
      k = which(brazil$numcri == i)
      time_i = brazil$tstop[k]
      age_i = brazil$age_group[k]
      group_i = as.numeric(brazil$grupo[k]) - 1
      gen_i = as.numeric(brazil$gender[k]) - 1
      sanit_i = brazil$sanit[k]
      measles_i = brazil$measles[k] 
      Cens[l] = time_i[which(brazil$status[k] == 0)]
      time_i1 = time_i[which(brazil$status[k] == 1)]
      z_obsd[l,1] = group_i[which(brazil$status[k] == 0)]
      z_obsd[l,2] = gen_i[which(brazil$status[k] == 0)]
      z_obsd[l,3] = age_i[which(brazil$status[k] == 0)]
      z_obsd[l,4] = sanit_i[which(brazil$status[k] == 0)]
      z_obsd[l,5] = measles_i[which(brazil$status[k] == 0)]
      if(length(time_i1) == 0){
        R_obsd[l,] = rep(tau + 1, r)
      } else {
        R_obsd[l,1:length(time_i1)] = sort(time_i1)
        if(length(time_i1) < r) R_obsd[l,(length(time_i1) + 1): r] = tau + 1 
      }
      l = l +1
    }
    
   
     
    Y_mat <- create_y(R_u,Cens)
    Delta <- as.numeric(R_obsd[,1]<Cens)
   
   set.seed(12345)
   zi <- rbinom(ncol(Y_mat),1,0.2)
   w_mat <- create_wt(R_u,Delta,Y_mat,zi)
   temp_out = rep(0, 2)
   z_obsd1 <- as.matrix(z_obsd[,3:4])
   temp_out <- try(beta_estimation(R_obsd,R_u,Cens,Y_mat,z_obsd1,Delta,w_mat,
                                   initial_beta = rep(0, 2)), silent = T)
   beta_ests=temp_out
     dN <- create_dn(R_u,R_obsd,Cens)
     dmu00 <- create_mu(R_u,w_mat,Y_mat,z_obsd1,beta_ests,dN)
     dmu01 <- create_mu(R_u,w_mat,Y_mat,z_obsd1,beta_ests,dN)*exp(beta_ests[2])
     dmu10 <- create_mu(R_u,w_mat,Y_mat,z_obsd1,beta_ests,dN)*exp(beta_ests[1])
     dmu11 <- create_mu(R_u,w_mat,Y_mat,z_obsd1,beta_ests,dN)*exp(sum(beta_ests))
     
     mu00 = cumsum(dmu00)
     mu01 = cumsum(dmu01)
     mu10 = cumsum(dmu10)
     mu11 = cumsum(dmu11)
     data = as.data.frame(cbind(R_u, mu00, mu01, mu10, mu11))
boot.num = 1000
data_ = NULL 

ids = 1:nrow(R_obsd)
ids_event = ids[which(Delta == 1)]
ids_noevent = ids[which(Delta == 0)]
  coh.size = nrow(R_obsd)
   coh.size1 = length(ids_event)
   coh.size2 = length(ids_noevent)
   
set.seed(12345)
for(loop in 1:boot.num){
  loc1 = sample.int(n=coh.size1,size=coh.size1,replace=T) 
  loc2 <- sample.int(n=coh.size2,size=coh.size2,replace=T) 
  loc = c(loc1, loc2)
  R_obsd1 = R_obsd[loc,]
  Cens1 = Cens[loc]
  z_obsd11 = z_obsd1[loc,]
  R_u1 <- create_r_u(R_obsd1, Cens1)
  Y_mat1 <- create_y(R_u1,Cens1)
  Delta1 <- as.numeric(R_obsd1[,1]<Cens1)
  zi1 <- rbinom(ncol(Y_mat1),1,0.2) #0.1
  eta1<- rep(0,ncol(Y_mat1))
  sub1<- which(Delta1*(1-zi1)==1)
  eta1[sub1]<-rbinom(length(sub1),1,0.1)
  w_mat1 <- create_wt(R_u1,Delta1,Y_mat1,zi1)
  temp_out <- try(beta_estimation(R_obsd1,R_u1,Cens1,Y_mat1,z_obsd11,Delta1,w_mat1,
                                  initial_beta = rep(0, 2)))
  if(length(temp_out)>1){
      beta_ests=temp_out
      dN <- create_dn(R_u1,R_obsd1,Cens1)
      dmu00 <- create_mu(R_u1,w_mat1,Y_mat1,z_obsd11,beta_ests,dN)
      dmu01 <- create_mu(R_u1,w_mat1,Y_mat1,z_obsd11,beta_ests,dN)*exp(beta_ests[2])
      dmu10 <- create_mu(R_u1,w_mat1,Y_mat1,z_obsd11,beta_ests,dN)*exp(beta_ests[1])
      dmu11 <- create_mu(R_u1,w_mat1,Y_mat1,z_obsd11,beta_ests,dN)*exp(sum(beta_ests))
      
      mu00 = cumsum(dmu00)
      mu01 = cumsum(dmu01)
      mu10 = cumsum(dmu10)
      mu11 = cumsum(dmu11)
      data1 = as.data.frame(cbind(R_u1, mu00, mu01, mu10, mu11))
      
    }else{
      data1 = NULL
    }

  data_ = bind_rows(data_, data1)
    
  }
 
data_2 = as.data.frame(data_) 
write.csv(data_2, "data_step1_new4.csv", row.names = F)

data_1 <- data_ %>% group_by(R_u1) %>%
      summarize_all(funs(sd),na.rm = T) %>% rename(R_u = R_u1)
  
data11 <- data %>% left_join(data_1, by = "R_u")
data11 <- data11 %>% mutate(mu00_ll = mu00.x - 1.96*mu00.y, mu00_ul = mu00.x + 1.96*mu00.y,
                          mu01_ll = mu01.x - 1.96*mu01.y, mu01_ul = mu01.x + 1.96*mu01.y,
                          mu10_ll = mu10.x - 1.96*mu10.y, mu10_ul = mu10.x + 1.96*mu10.y,
                          mu11_ll = mu11.x - 1.96*mu11.y, mu11_ul = mu11.x + 1.96*mu11.y)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
               "#CC79A7")
data_11 <- data11 %>% select(-c(mu00.y, mu01.y, mu10.y, mu11.y)) %>% 
  gather(key = key, value = value, -R_u)
data_11 <- data_11 %>% 
  mutate(key1 = ifelse(grepl("mu00_", key)== T, "mu00_cb",
                       ifelse(grepl("mu01_", key)== T, "mu01_cb", 
                              ifelse(grepl("mu10_", key)== T, "mu10_cb", 
                                     ifelse(grepl("mu11_", key)== T, "mu11_cb", key)))))
data_11 <- data_11 %>% 
  mutate(key2 = ifelse(grepl("mu00", key)== T, "group1", 
                       ifelse(grepl("mu01", key)== T, "group2",
                              ifelse(grepl("mu10", key)== T, "group3","group4"))),
                   key3 = ifelse(grepl("cb", key1)== T, "CI", "notCI"), 
                   key4 = ifelse(grepl("ul", key)== T, "CI_UL", 
                                 ifelse(grepl("ll", key)== T,"CI_LL", "notCI")))

write.csv(data_11, "data_for_plot_new4.csv", row.names = F)

if(F){
  
  val <- c("< 12 months (estimate)", "95% CI", ">= 12 months (estimate)",
           "95% CI")
  val1 <- c("< 12 months (estimate)", "95% CI", "95% Confidence Intervals",
            ">= 12 months (estimate)","95% Confidence Intervals", "95% Confidence Intervals")
  val2 <- c(">= 12 months","< 12 months")
  val3 <- c("Lower 95% CI", "Upper 95% CI", "Estimate")
  val3_1 <- c("95% CI", "Estimate")
  val2_1 <- c("Absence", "Presence")
  val22 <- c("< 12 months / Absence of Toilet", ">= 12 months / Presence of Toilet")
  
  p <- NULL
  
  for(i in 1:3){
    if(i == 1){
      data_new = data_31
      p[[i]] = ggplot(data = data_new, aes(x = R_u, y = value,color = key2, linetype = key4)) + 
        geom_line(size = 1.25) + 
        ggtitle("Dichotomized Age") +
        scale_colour_manual(name = "Age Groups", labels = val2, values=c(cbPalette[7], cbPalette[8])) + 
        scale_linetype_manual(values=c(6,6,1), guide = F)+
        scale_x_continuous(name = "Time (days)") +
        scale_y_continuous(name = "Cumulative Rates Estimates") +
        theme(axis.text=element_text(size=20), axis.title=element_text(size=15),
              plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
              legend.text=element_text(size=20),
              legend.title=element_text(size=25), legend.position = "bottom") 
    }else if(i == 2){
      data_new = data_41
      p[[i]] =   ggplot(data = data_new, aes(x = R_u, y = value,color = key2, linetype = key4)) + 
        geom_line(size = 1.25) + 
        ggtitle("Sanitation: Presence of Toilet at Home") +
        scale_colour_manual(name = "Sanitation Groups", labels = val2_1, values=c(cbPalette[7], cbPalette[8])) + 
        scale_linetype_manual(values=c(6,6,1), guide = F)+
        scale_x_continuous(name = "Time (days)") +
        scale_y_continuous(name = "Cumulative Rates Estimates") +
        theme(axis.text=element_text(size=20), axis.title=element_text(size=15),
              plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
              legend.text=element_text(size=20),
              legend.title=element_text(size=25), legend.position = "bottom")
    }else{
      data_new = data_41
      p[[i]] =   ggplot(data = data_new, aes(x = R_u, y = value,color = key2, linetype = key3)) + 
        geom_line(size = 1.25) + 
        ggtitle("Sanitation: Presence of Toilet at Home") +
        scale_colour_manual(values=c(cbPalette[7], cbPalette[8]), guide = F) + 
        scale_linetype_manual(name = "Estimate and Confidence Bands", labels = val3_1, 
                              values=c(6,1))+
        scale_x_continuous(name = "Time (days)") +
        scale_y_continuous(name = "Cumulative Rates Estimates") +
        theme(axis.text=element_text(size=20), axis.title=element_text(size=15),
              plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
              legend.text=element_text(size=20),
              legend.title=element_text(size=25), legend.position = "bottom") 
    }
    
  }
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  mylegend1<-g_legend(p[[1]])
  mylegend2<-g_legend(p[[2]])
  mylegend3<-g_legend(p[[3]])
  
  print(grid.arrange(p[[1]], p[[2]], nrow=1, 
                     top = textGrob(paste0("Risk Factor Specific Rate Function Estimates and Bootstrapped 95% CI"),
                                    gp=gpar(fontface="bold", fontsize=20)))) 
  
  p3<- grid.arrange(arrangeGrob(p[[1]] + theme(legend.position="none"),
                                p[[2]] + theme(legend.position="none"),
                                nrow=1),
                    arrangeGrob(mylegend1,nrow=2),
                    arrangeGrob(mylegend2,nrow=3),
                    arrangeGrob(mylegend3,nrow=4),heights = c(10,1,1,1),
                    top = textGrob(paste0("Risk Factor Specific Rate Function Estimates and Bootstrapped 95% CI"),
                                   gp=gpar(fontface="bold", fontsize=20)))
  
  pdf("sigh_data_example_ratefn__new4.pdf", width = 18, height = 10)
  
  print(
    grid.arrange(arrangeGrob(p[[1]],p[[2]],nrow=1),
                 arrangeGrob(mylegend3,nrow=2),heights = c(7.5,1)
                
    )
    
  )
  
  dev.off()
  
  
}
