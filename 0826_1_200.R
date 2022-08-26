
library(BB)
library(stats)
library(elliptic)
library(nleqslv)

set.seed(777)
alln <- c(200, 400)

allbe <- c(0, 0.2)

allth <- c(0, 0.1)

allgm <- c(0, 0.3)

#for Hn, take 1 (constant for all t) or  sum(delta1)/n (proportional to subjects under observation)
### you may start with I=J=K=L=1, 
### so that n=100, beta=theta=gm=0


I=1

J=1
K=1
L=1


n=alln[I];#sample sizes
p0=0.5; #Zi generated from bernoulli with p=0.5

beta=allbe[J]; 
theta = allth[K];
gm=allgm[L]; 

### create vectors to save results on point estimates
mj <- 1000 #1000 replications
hb <- hth <- hgm <- numeric(mj) #point estimation
hbth=matrix(0,mj,2)

rej=numeric(mj)
### create vectors to save results on CP and estimated standard errors on beta and theta
cpth <- sgth  <- matrix(0, 2, mj)

### create vectors to save results on CP and estimated standard errors on gm
cpgm <- sggm <- numeric(mj)

### largest followup time 
tau <- 2

iter = 0
ct = 1

### time points for nonparametric baseline estimates
T0 <- seq(0, tau, length=1000)
lmd0 = 0.1 #baseline hazard 0.1

### true Lambda0(t) (for eq.7)
Lmd0 <-  T0^1.25 * 4/5
NT0 <- length(T0)
### matrix storing esimated Lambda(t) in 1000 replicates by rows
Lmd_sig_est = Lmd_est=matrix(0, mj, NT0);

library("simsurv")
while (ct <= mj)	{
  iter = iter+1
  X=runif(n, 0, 1);
  Z=rbinom(n, 1, 0.5);
  
  ### generate C (T)
  betas = data.frame(bs=rep(beta, n), ts=rep(theta, n))
  covdat <- data.frame(X, Z)

  
  haz <- function(t, x, betas) {
    (1 +1*t^4*(x[["X"]])^2+ betas[["bs"]] * x[["X"]]) * exp(betas[["ts"]] * x[["Z"]])
  }
  
  
  
  C <- simsurv(hazard = haz,  x = covdat, betas = betas, maxt = tau)$eventtime
  
  C# First let W = Z
  lmd = lmd0 + gm * Z 
  
  ### generate Death times (U) 
  D <- numeric(n)
  for (i in 1:n)
    D[i] <-  rexp(1, lmd[i])
  
  ### constant censoring times
  C1 = tau
  
  ### observed TT=T1
  TT <- obsT <- apply(cbind(C, C1, D), 1, min)
  ### dtCD = delta1
  dtCD <- as.numeric(C<=D) * as.numeric(C<=C1)
  ### obc: observed C when dtCD==1
  obc <- sort(unique(C[dtCD==1]))
  ### pooled T1
  allt <- sort(unique(obsT))
  nt <- length(allt)
  ### DT: risk set indicator by subjects (rows) over allt (columns)
  DT = as.matrix(matrix(rep(allt, n), n, nt, byrow=TRUE)<= matrix(rep(TT, nt), n, nt))
  
  ### survival data
  # dt1: delta2
  dt1 <- as.numeric(D<=C1)
  # TT1: T2
  TT1 <- apply(cbind(C1, D), 1, min)
  # allt1: pooled T2
  allt1 <- sort(unique(TT1))
  #dtime <- sort(unique(D[dt1==1]))
  ntd <- length(allt1)
  
  ### DT1: risk set indicator by subjects (rows) over allt1 (columns)
  DT1 <- as.matrix(matrix(rep(allt1, n), n, ntd, byrow=TRUE)<= matrix(rep(TT1, ntd), n, ntd))
  
  
  ### dNd for eq.10
  dNd <- matrix(0, n, ntd)
  for (j in 1:ntd)
    dNd[,j] <- dt1 * as.numeric(allt1[j]==TT1)
  
  mzd <- matrix(rep(Z, ntd), n, ntd)
  
  wbar <- apply(DT1 * mzd, 2, sum) / apply(DT1, 2, sum)
  
  #dallt and dallt1
  
  dallt=rep(0,nt)
  dallt[1]=allt[1]
  for (j in 2:nt){
    dallt[j]=allt[j]-allt[j-1]
  }
  
  
  dallt1=rep(0,ntd)
  dallt1[1]=allt1[1]
  for (j in 2:ntd){
    dallt1[j]=allt1[j]-allt1[j-1]
  }
  
  
  #estimate gamma
  leftgm=matrix(0,n,ntd)
  for(i in 1:n){
    for (j in 1:ntd){
      leftgm[i,j]=DT1[i,j]*(Z[i]-wbar[j])*dNd[i,j]
    }
  }
  leftgm1=sum(leftgm)
  
  rightgm=matrix(0,n,ntd)
  fgamma=function(gamma){
    for(i in 1:n){
      rightgm[i,1]=1/2*DT1[i,1]*(Z[i]-wbar[1])*gamma*Z[i]*dallt1[1]
      for (j in 2:ntd){
        
        rightgm[i,j]=1/2*(DT1[i,j-1]*(Z[i]-wbar[j-1])*gamma*Z[i]+DT1[i,j]*(Z[i]-wbar[j])*gamma*Z[i])*dallt1[j]
        
      }
    }
    rightgm1=sum(rightgm)
    return(leftgm1-rightgm1)
  }
  
  
  #hgm[ct]=nleqslv(0, fgamma,control=list(trace=1,btol=.01,delta="newton"))$x
  hgm[ct]=dfsane(0, fgamma,control=list(maxit=2500))$par
  
  ##estimate shat estimate Deltahat with allt,length=nt
  
  Ndd <- matrix(0, n, nt)
  for (i in 1:n){
    for (j in 1:nt)
      Ndd[i,j] <- dt1[i] * as.numeric(allt[j]>=TT1[i])
  }
  
  
  
  dNdd <- matrix(0, n, nt)
  
  for(i in 1:n){
    dNdd[i,1] <- dt1[i] * as.numeric(Ndd[i,1]==1)
    
    for (j in 2:nt)
      dNdd[i,j] <- dt1[i]* as.numeric(Ndd[i,j]==1)*as.numeric(Ndd[i,j-1]==0)
    
  }
  
  mz <- matrix(rep(Z, nt), n, nt)
  
  
  DT3 <- as.matrix(matrix(rep(allt, n), n, nt, byrow=TRUE)<= matrix(rep(TT1, nt), n, nt))
  
  
  
  
  
  shat=matrix(0,n,nt)
  
  part0=apply(DT3*dNdd,2,sum)/ apply(DT3, 2, sum)
  part1=rep(0,nt)
  part1[1]=part0[1]
  for(j in 2:nt){
    part1[j]=part0[j]+part1[j-1]
  }
  
  
  
  component6=apply(DT3*hgm[ct]*mz,2,sum)/ apply(DT3, 2, sum)
  
  part2=rep(0,nt)
  
  part2[1]=1/2*component6[1]*dallt[1]
  for(j in 2:nt){
    part2[j]=1/2*(component6[j]+component6[j-1])*dallt[j]+part2[j-1]
  }
  
  
  
  
  for(i in 1:n){
    for(j in 1:nt){
      shat[i,j]=exp(-part1[j]+part2[j]-hgm[ct]*Z[i]*allt[j])
    }
  }
  
  
  
  #estimate phihat
  
  m1 <- matrix(rep(1, nt), n, nt)
  
  phihat <- DT * m1/shat
  
  mx <- matrix(rep(X, nt), n, nt)
  xbar <- apply(phihat * mx, 2, sum) / apply(phihat, 2, sum)
  zbar <- apply(phihat * mz, 2, sum) / apply(phihat, 2, sum)
  
  V=rbind(X,Z)
  vbar=rbind(xbar,zbar)
  
  
  
  #estimate beta and theta from eq.8
  
  ### dN for eq.8
  
  N <- matrix(0, n, nt)
  
  for (i in 1:n){
    for (j in 1:nt)
      N[i,j] <- dtCD[i] * as.numeric(allt[j]>=TT[i])
  }
  
  
  
  dN <- matrix(0, n, nt)
  
  for (i in 1:n){
    dN[i,1] <- dtCD[i] * as.numeric(N[i,1]==1)
    
    for (j in 2:nt)
      dN[i,j] <- dtCD[i] * as.numeric(N[i,j]==1)*as.numeric(N[i,j-1]==0)
  }
  
  
  
  
  #take H=1 here.
  #or use
  
  #dN <- matrix(0, n, nt)
  #for (j in 1:nt)
  #dN[,j] <- dtCD * as.numeric(allt[j]==TT)
  
  
  
  
  fth1=function(beth){
    
    leftth1=matrix(0,n,nt)
    for(i in 1:n){
      for (j in 1:nt){
        leftth1[i,j]=(V[1,i]-vbar[1,j])*phihat[i,j]*exp(-beth[2]*Z[i])*dN[i,j]
      }
    }
    leftth11=sum(leftth1)
    
    rightth1=matrix(0,n,nt)
    for(i in 1:n){
      rightth1[i,1]=1/2*(V[1,i]-vbar[1,1])*phihat[i,1]*beth[1]*X[i]*dallt[1]
      for (j in 2:nt){
        rightth1[i,j]=1/2*((V[1,i]-vbar[1,j])*phihat[i,j]*beth[1]*X[i]+(V[1,i]-vbar[1,j-1])*phihat[i,j-1]*beth[1]*X[i])*dallt[j]
      }
    }
    rightth11=sum(rightth1)
    
    return(leftth11-rightth11)
  }
  
  fth2=function(beth){
    
    leftth2=matrix(0,n,nt)
    for(i in 1:n){
      for (j in 1:nt){
        leftth2[i,j]=(V[2,i]-vbar[2,j])*phihat[i,j]*exp(-beth[2]*Z[i])*dN[i,j]
      }
    }
    leftth22=sum(leftth2)
    
    rightth2=matrix(0,n,nt)
    for(i in 1:n){
      rightth2[i,1]=1/2*(V[2,i]-vbar[2,1])*phihat[i,1]*beth[1]*X[i]*dallt[1]
      for (j in 2:nt){
        rightth2[i,j]=1/2*((V[2,i]-vbar[2,j])*phihat[i,j]*beth[1]*X[i]+(V[2,i]-vbar[2,j-1])*phihat[i,j-1]*beth[1]*X[i])*dallt[j]
      }
    }
    rightth22=sum(rightth2)
    
    return(leftth22-rightth22)
  }
  
  
  fgm_0=function(beth){
    y=numeric(2)
    y[1]=fth1(beth)
    y[2]=fth2(beth)
    return(y)
  }
  
  
  hbth[ct,]=dfsane(par=c(0,0), fn=fgm_0, control=list(maxit=3000))$par
  
  

  #or use dfsane(par=c(0,0), fn=fgm_0, control=list(trace=FALSE))
  
  
  #############evaluate the LMD and lmd (eq. 7) when N_i(t) has a jump. ###############
  ###use NT0 as the time point and plug the estimated LMD value(time point is nt)
  
  
  
  partlmd1=rep(0,nt)
  partlmd11=matrix(0,n,nt)
  
  for (i in 1:n){
    for (j in 1:nt){
      partlmd11[i,j]=phihat[i,j]*exp(-hbth[ct,2]*Z[i])*dN[i,j]
    }
    partlmd1=apply(partlmd11,2,sum)/apply(phihat, 2, sum)
  }
  
  partlmd2=rep(0,nt)
  
  partlmd22=matrix(0,n,nt)
  
  
  for (i in 1:n){
    partlmd22[i,1]=1/2*phihat[i,1]*hbth[ct,1]*X[i]*dallt[1]
    
    for (j in 2:nt)
      partlmd22[i,j]=1/2*(phihat[i,j]*hbth[ct,1]*X[i]+phihat[i,j-1]*hbth[ct,1]*X[i])*dallt[j]
  }
  
  partlmd2=apply(partlmd22, 2, sum)/apply(phihat,2,sum)
  
  
  
  
  ###########estimate SEE based on theorems##################
  
  xbar1 <- apply(phihat * mx, 2, sum) / apply(phihat, 2, sum)
  zbar1 <- apply(phihat * mz, 2, sum) / apply(phihat, 2, sum)
  vbar1 <- rbind(xbar1,zbar1)
  wbar1 <-apply(DT3 * mz, 2, sum) / apply(DT3, 2, sum)
  
  
  
  #estimate q1
  
  q1part1=matrix(0,n,nt)
  
  for (i in 1:n){
    for (j in 1:nt){
      q1part1[i,j]=(V[1,i]-vbar1[1,j])*phihat[i,j]*exp(-hbth[ct,2]*Z[i])*dN[i,j]
    }
  }
  
  q1part2=matrix(0,n,nt)
  for (i in 1:n){
    q1part2[i,1]=1/2*(V[1,i]-vbar1[1,1])*phihat[i,1]*(hbth[ct,1]*X[i])*dallt[1]
    
    for (j in 2:nt){
      q1part2[i,j]=1/2*((V[1,i]-vbar1[1,j])*phihat[i,j]*(hbth[ct,1]*X[i])+(V[1,i]-vbar1[1,j-1])*phihat[i,j-1]*(hbth[ct,1]*X[i]))*dallt[j]
    }
  }
  
  q1part3=matrix(0,n,nt)
  for (i in 1:n){
    
    for (j in 1:nt){
      q1part3[i,j]=(V[1,i]-vbar1[1,j])*phihat[i,j]*partlmd1[j]
    }
  }
  
  
  
  
  
  q1part4=matrix(0,n,nt)
  for (i in 1:n){
    
    for (j in 1:nt){
      q1part4[i,j]=(V[1,i]-vbar1[1,j])*phihat[i,j]*partlmd2[j]
    }
  }
  
  q11=apply(q1part1-q1part2-q1part3+q1part4,1,sum)
  
  q1part11=matrix(0,n,nt)
  
  for (i in 1:n){
    for (j in 1:nt){
      q1part11[i,j]=(V[2,i]-vbar1[2,j])*phihat[i,j]*exp(-hbth[ct,2]*Z[i])*dN[i,j]
    }
  }
  
  q1part22=matrix(0,n,nt)
  for (i in 1:n){
    q1part22[i,1]=1/2*(V[2,i]-vbar1[2,1])*phihat[i,1]*(hbth[ct,1]*X[i])*dallt[1]
    
    for (j in 2:nt){
      q1part22[i,j]=1/2*((V[2,i]-vbar1[2,j])*phihat[i,j]*(hbth[ct,1]*X[i])+(V[2,i]-vbar1[2,j-1])*phihat[i,j-1]*(hbth[ct,1]*X[i]))*dallt[j]
    }
  }
  
  q1part33=matrix(0,n,nt)
  for (i in 1:n){
    
    for (j in 1:nt){
      q1part33[i,j]=(V[2,i]-vbar1[2,j])*phihat[i,j]*partlmd1[j]
    }
  }
  
  
  
  
  q1part44=matrix(0,n,nt)
  for (i in 1:n){
    
    for (j in 1:nt){
      q1part44[i,j]=(V[2,i]-vbar1[2,j])*phihat[i,j]*partlmd2[j]
    }
  }
  
  q22=apply(q1part11-q1part22-q1part33+q1part44,1,sum)
  
  q1=rbind(q11,q22)
  #estimate q2
  #estimate PHI
  
  
  
  
  
  sumPHIpart1=apply(q1part1,1,sum)
  
  PHIpart1=matrix(0,n,nt)
  minus1=matrix(0,n,nt)
  for (i in 1:n){
    minus1[i,1]=q1part1[i,1]
    PHIpart1[i,1]=sumPHIpart1[i]-minus1[i,1]
    for (j in 2:nt){
      minus1[i,j]=minus1[i,j-1]+q1part1[i,j]
      PHIpart1[i,j]=sumPHIpart1[i]-minus1[i,j]
    }
  }
  
  sumPHIpart2=apply(q1part2,1,sum)
  
  PHIpart2=matrix(0,n,nt)
  minus2=matrix(0,n,nt)
  for (i in 1:n){
    minus2[i,1]=q1part2[i,1]
    PHIpart2[i,1]=sumPHIpart2[i]-minus2[i,1]
    for (j in 2:nt){
      minus2[i,j]=minus2[i,j-1]+q1part2[i,j]
      PHIpart2[i,j]=sumPHIpart2[i]-minus2[i,j]
    }
  }
  
  sumPHIpart3=apply(q1part3,1,sum)
  
  PHIpart3=matrix(0,n,nt)
  minus3=matrix(0,n,nt)
  for (i in 1:n){
    minus3[i,1]=q1part3[i,1]
    PHIpart3[i,1]=sumPHIpart3[i]-minus3[i,1]
    for (j in 2:nt){
      minus3[i,j]=minus3[i,j-1]+q1part3[i,j]
      PHIpart3[i,j]=sumPHIpart3[i]-minus3[i,j]
    }
  }
  
  sumPHIpart4=apply(q1part4,1,sum)
  
  PHIpart4=matrix(0,n,nt)
  minus4=matrix(0,n,nt)
  for (i in 1:n){
    minus4[i,1]=q1part4[i,1]
    PHIpart4[i,1]=sumPHIpart4[i]-minus4[i,1]
    for (j in 2:nt){
      minus4[i,j]=minus4[i,j-1]+q1part4[i,j]
      PHIpart4[i,j]=sumPHIpart4[i]-minus4[i,j]
    }
  }
  #PHihat also should have two dimension (x,z)
  PHIhat1=apply(PHIpart1-PHIpart2-PHIpart3+PHIpart4,2,sum)/n
  
  sumPHIpart11=apply(q1part11,1,sum)
  
  PHIpart11=matrix(0,n,nt)
  minus11=matrix(0,n,nt)
  for (i in 1:n){
    minus11[i,1]=q1part11[i,1]
    PHIpart11[i,1]=sumPHIpart11[i]-minus11[i,1]
    for (j in 2:nt){
      minus11[i,j]=minus11[i,j-1]+q1part11[i,j]
      PHIpart11[i,j]=sumPHIpart11[i]-minus11[i,j]
    }
  }
  
  sumPHIpart22=apply(q1part22,1,sum)
  
  PHIpart22=matrix(0,n,nt)
  minus22=matrix(0,n,nt)
  for (i in 1:n){
    minus22[i,1]=q1part22[i,1]
    PHIpart22[i,1]=sumPHIpart22[i]-minus22[i,1]
    for (j in 2:nt){
      minus22[i,j]=minus22[i,j-1]+q1part22[i,j]
      PHIpart22[i,j]=sumPHIpart22[i]-minus22[i,j]
    }
  }
  
  sumPHIpart33=apply(q1part33,1,sum)
  
  PHIpart33=matrix(0,n,nt)
  minus33=matrix(0,n,nt)
  for (i in 1:n){
    minus33[i,1]=q1part33[i,1]
    PHIpart33[i,1]=sumPHIpart33[i]-minus33[i,1]
    for (j in 2:nt){
      minus33[i,j]=minus33[i,j-1]+q1part33[i,j]
      PHIpart33[i,j]=sumPHIpart33[i]-minus33[i,j]
    }
  }
  
  sumPHIpart44=apply(q1part44,1,sum)
  
  PHIpart44=matrix(0,n,nt)
  minus44=matrix(0,n,nt)
  for (i in 1:n){
    minus44[i,1]=q1part44[i,1]
    PHIpart44[i,1]=sumPHIpart44[i]-minus44[i,1]
    for (j in 2:nt){
      minus44[i,j]=minus44[i,j-1]+q1part44[i,j]
      PHIpart44[i,j]=sumPHIpart44[i]-minus44[i,j]
    }
  }
  
  PHIhat2=apply(PHIpart11-PHIpart22-PHIpart33+PHIpart44,2,sum)/n
  PHIhat=rbind(PHIhat1,PHIhat2)
  
  #estimate q2
  #estimate R0
  R0 <- apply(DT3 * m1, 2, sum) / n
  
  component5=matrix(0,2,nt)
  
  for(i in 1:2){
    for(j in 1:nt){
      component5[i,j]=PHIhat[i,j]/R0[j]
    }
  }
  
  q2part1=matrix(0,n,nt)
  
  for (i in 1:n){
    for (j in 1:nt){
      q2part1[i,j]=component5[1,j]*dNdd[i,j]
    }
  }
  
  q2part2=matrix(0,n,nt)
  for (i in 1:n){
    q2part2[i,1]=1/2*component5[1,1]*DT3[i,1]*hgm[ct]*Z[i]*dallt[1]
    
    for (j in 2:nt){
      q2part2[i,j]=1/2*(component5[1,j]*DT3[i,j]*hgm[ct]*Z[i]+component5[1,j-1]*DT3[i,j-1]*hgm[ct]*Z[i])*dallt[j]
    }
  }
  
  
  q2part3=matrix(0,n,nt)


  dDelta0hatpart1=apply(dNdd,2,sum)/apply(DT3,2,sum)
  
  for (i in 1:n){
    for (j in 1:nt){
      q2part3[i,j]=component5[1,j]*DT3[i,j]*dDelta0hatpart1[j]
    }
  }
  
  q2part4=matrix(0,n,nt)
  
  dDelta0hatpart22=matrix(0,n,nt)
  for(i in 1:n){
    dDelta0hatpart22[i,1]=1/2*DT3[i,j]*hgm[ct]*Z[i]*dallt[1]
    
    for (j in 2:nt){
      dDelta0hatpart22[i,j]=1/2*(DT3[i,j]*hgm[ct]*Z[i]+DT3[i,j-1]*hgm[ct]*Z[i])*dallt[j]
    }
  }
  
  dDelta0hatpart2=apply(dDelta0hatpart22,2,sum)/apply(DT3,2,sum)
  
  for (i in 1:n){
    
    for (j in 1:nt){
      q2part4[i,j]=component5[1,j]*DT3[i,j]*dDelta0hatpart2[j]
    }
  }
  
  q21=apply(q2part1-q2part2-q2part3+q2part4,1,sum)
  
  
  q2part11=matrix(0,n,nt)
  
  for (i in 1:n){
    for (j in 1:nt){
      q2part11[i,j]=component5[2,j]*dNdd[i,j]
    }
  }
  
  q2part22=matrix(0,n,nt)
  for (i in 1:n){
    q2part22[i,1]=1/2*component5[2,1]*DT3[i,1]*hgm[ct]*Z[i]*dallt[1]
    
    for (j in 2:nt){
      q2part22[i,j]=1/2*(component5[2,j]*DT3[i,j]*hgm[ct]*Z[i]+component5[2,j-1]*DT3[i,j-1]*hgm[ct]*Z[i])*dallt[j]
    }
  }
  
  
  q2part33=matrix(0,n,nt)
  
  
  for (i in 1:n){
    for (j in 1:nt){
      q2part33[i,j]=component5[2,j]*DT3[i,j]*dDelta0hatpart1[j]
    }
  }
  
  q2part44=matrix(0,n,nt)
  
  
  for (i in 1:n){
    
    for (j in 1:nt){
      q2part44[i,j]=component5[2,j]*DT3[i,j]*dDelta0hatpart2[j]
    }
  }
  
  q22=apply(q2part11-q2part22-q2part33+q2part44,1,sum)
  
  q2=rbind(q21,q22)
  #estimate q3
  
  
  #estimate Eta
  
  Eta=rep(0,nt)
  Eta[1]=1/2*wbar1[1]*dallt[1]
  for (j in 2:nt){
    Eta[j]=1/2*(wbar1[j]+wbar1[j-1])*dallt[j]+Eta[j-1]
  }
  
  ETA=-Eta
  
  
  #estimate Omega
  
  
  
  product=matrix(0,n,nt)  
  for(i in 1:n){
    product[i,1]=1/2*(Z[i]-wbar1[1])%*%t(Z[i]-wbar1[1])*DT3[i,1]*dallt[1]
    
    for(j in 2:nt){
      product[i,j]=1/2*((Z[i]-wbar1[j])%*%t(Z[i]-wbar1[j])*DT3[i,j]+(Z[i]-wbar1[j-1])%*%t(Z[i]-wbar1[j-1])*DT3[i,j-1])*dallt[j]
    }
    
  }
  
  
  
  Omega=sum(product)/n
  
  
  #estimate B
  Bpart1=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      Bpart1[i,j]=(V[1,i]-vbar1[1,j])*(ETA[j]+allt[j]*Z[i])*((Omega)^(-1))*phihat[i,j]*exp(-hbth[ct,2]*Z[i])*dN[i,j]
    }
  }
  
  Bpart2=matrix(0,n,nt)
  for (i in 1:n){
    Bpart2[i,1]=1/2*(V[1,i]-vbar1[1,1])*(ETA[1]+allt[1]*Z[i])*((Omega)^(-1))*phihat[i,1]*(hbth[ct,1]*X[i])*dallt[1]
    
    for (j in 2:nt){
      Bpart2[i,j]=1/2*((V[1,i]-vbar1[1,j])*(ETA[j]+allt[j]*Z[i])*((Omega)^(-1))*phihat[i,j]*(hbth[ct,1]*X[i])
                       +(V[1,i]-vbar1[1,j-1])*(ETA[j-1]+allt[j-1]*Z[i])*((Omega)^(-1))*phihat[i,j-1]*(hbth[ct,1]*X[i]))*dallt[j]
    }
  }
  
  
  Bpart3=matrix(0,n,nt)
  
  
  for (i in 1:n){
    for (j in 1:nt){
      Bpart3[i,j]=(V[1,i]-vbar1[1,j])*(ETA[j]+allt[j]*Z[i])*((Omega)^(-1))*phihat[i,j]*partlmd1[j]
    }
  }
  
  Bpart4=matrix(0,n,nt)
  
  for (i in 1:n){
    
    for (j in 1:nt){
      Bpart4[i,j]=(V[1,i]-vbar1[1,j])*(ETA[j]+allt[j]*Z[i])*((Omega)^(-1))*phihat[i,j]*partlmd2[j]
    }
    
  }
  
  
  
  B1=sum(Bpart1-Bpart2-Bpart3+Bpart4)/n
  
  
  Bpart11=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      Bpart11[i,j]=(V[2,i]-vbar1[2,j])*(ETA[j]+allt[j]*Z[i])*((Omega)^(-1))*phihat[i,j]*exp(-hbth[ct,2]*Z[i])*dN[i,j]
    }
  }
  
  Bpart22=matrix(0,n,nt)
  for (i in 1:n){
    Bpart22[i,1]=1/2*(V[2,i]-vbar1[2,1])*(ETA[1]+allt[1]*Z[i])*((Omega)^(-1))*phihat[i,1]*(hbth[ct,1]*X[i])*dallt[1]
    
    for (j in 2:nt){
      Bpart22[i,j]=1/2*((V[2,i]-vbar1[2,j])*(ETA[j]+allt[j]*Z[i])*((Omega)^(-1))*phihat[i,j]*(hbth[ct,1]*X[i])
                        +(V[2,i]-vbar1[2,j-1])*(ETA[j-1]+allt[j-1]*Z[i])*((Omega)^(-1))*phihat[i,j-1]*(hbth[ct,1]*X[i]) )*dallt[j]
    }
  }
  
  
  Bpart33=matrix(0,n,nt)
  
  
  for (i in 1:n){
    for (j in 1:nt){
      Bpart33[i,j]=(V[2,i]-vbar[2,j])*(ETA[j]+allt[j]*Z[i])*((Omega)^(-1))*phihat[i,j]*partlmd1[j]
    }
  }
  
  Bpart44=matrix(0,n,nt)
  
  for (i in 1:n){
    
    for (j in 1:nt){
      Bpart44[i,j]=(V[2,i]-vbar[2,j])*(ETA[j]+allt[j]*Z[i])*((Omega)^(-1))*phihat[i,j]*partlmd2[j]
      
    }
    
  }
  
  B2=sum(Bpart11-Bpart22-Bpart33+Bpart44)/n
  
  Bhat=rbind(B1,B2)
  
  #then estimate q3
  
  q3part1=matrix(0,n,nt)
  
  for (i in 1:n){
    for (j in 1:nt){
      q3part1[i,j]=Bhat[1]*(Z[i]-wbar1[j])*dNdd[i,j]
    }
  }
  
  q3part2=matrix(0,n,nt)
  for (i in 1:n){
    q3part2[i,1]=1/2*Bhat[1]*(Z[i]-wbar1[1])*DT3[i,1]*hgm[ct]*Z[i]*dallt[1]
    
    for (j in 2:nt){
      q3part2[i,j]=1/2*(Bhat[1]*(Z[i]-wbar1[j])*DT3[i,j]*hgm[ct]*Z[i]+Bhat[1]*(Z[i]-wbar1[j-1])*DT3[i,j-1]*hgm[ct]*Z[i])*dallt[j]
    }
  }
  
  q3part3=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      q3part3[i,j]=Bhat[1]*(Z[i]-wbar1[j])*DT3[i,j]*dDelta0hatpart1[j]
    }
  }
  
  q3part4=matrix(0,n,nt)
  for (i in 1:n){
    
    for (j in 1:nt){
      q3part4[i,j]=Bhat[1]*(Z[i]-wbar1[j])*DT3[i,j]*dDelta0hatpart2[j]
    }
  }
  
  
  q31=apply(q3part1-q3part2-q3part3+q3part4,1,sum)
  
  
  q3part11=matrix(0,n,nt)
  
  for (i in 1:n){
    for (j in 1:nt){
      q3part11[i,j]=Bhat[2]*(Z[i]-wbar1[j])*dNdd[i,j]
    }
  }
  
  q3part22=matrix(0,n,nt)
  for (i in 1:n){
    q3part22[i,1]=Bhat[2]*(Z[i]-wbar1[1])*DT3[i,1]*hgm[ct]*Z[i]*dallt[1]
    
    for (j in 2:nt){
      q3part22[i,j]=1/2*(Bhat[2]*(Z[i]-wbar1[j])*DT3[i,j]*hgm[ct]*Z[i]
                         +Bhat[2]*(Z[i]-wbar1[j-1])*DT3[i,j-1]*hgm[ct]*Z[i])*dallt[j]
    }
  }
  
  q3part33=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      q3part33[i,j]=Bhat[2]*(Z[i]-wbar1[j])*DT3[i,j]*dDelta0hatpart1[j]
    }
  }
  
  q3part44=matrix(0,n,nt)
  for (i in 1:n){
    
    for (j in 1:nt){
      q3part44[i,j]=Bhat[2]*(Z[i]-wbar1[j])*DT3[i,j]*dDelta0hatpart2[j]
    }
  }
  
  
  q32=apply(q3part11-q3part22-q3part33+q3part44,1,sum)
  
  q3=rbind(q31,q32)
  
  
  #then estimate Sigmahat
  
  
  
  Sigma=matrix(0,n,4)
  for (i in 1:n){
    Sigma[i,1]=as.matrix((q1[,i]+q2[,i]+q3[,i])%*%t(q1[,i]+q2[,i]+q3[,i]))[1,1]
    Sigma[i,2]=as.matrix((q1[,i]+q2[,i]+q3[,i])%*%t(q1[,i]+q2[,i]+q3[,i]))[1,2]
    Sigma[i,3]=as.matrix((q1[,i]+q2[,i]+q3[,i])%*%t(q1[,i]+q2[,i]+q3[,i]))[2,1]
    Sigma[i,4]=as.matrix((q1[,i]+q2[,i]+q3[,i])%*%t(q1[,i]+q2[,i]+q3[,i]))[2,2]
  }
  
  Sigmahat=matrix(0,2,2)
  Sigmahat[1,1]=sum(Sigma[,1])/n
  Sigmahat[1,2]=sum(Sigma[,2])/n
  Sigmahat[2,1]=sum(Sigma[,3])/n
  Sigmahat[2,2]=sum(Sigma[,4])/n
  
  #estimate A
  
  #estimate A11
  A1=matrix(0,n,nt)
  for (i in 1:n){
    A1[i,1]=1/2*(X[i]-xbar1[1])*phihat[i,1]*X[i]*dallt[1]
    
    for (j in 2:nt){
      A1[i,j]=1/2*((X[i]-xbar1[j])*phihat[i,j]*X[i]+(X[i]-xbar1[j-1])*phihat[i,j-1]*X[i])*dallt[j]
    }
  }
  
  A11=sum(A1)/n
  
  #estimate A12
  
  A2=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A2[i,j]=(X[i]-xbar1[j])*phihat[i,j]*exp(-hbth[ct,2]*Z[i])*Z[i]*dN[i,j]
    }
  }
  
  A12=sum(A2)/n
  
  
  #estimate A21
  A3=matrix(0,n,nt)
  for (i in 1:n){
    A3[i,1]=1/2*(Z[i]-zbar1[1])*phihat[i,1]*X[i]*dallt[1]
    
    for (j in 2:nt){
      A3[i,j]=1/2*((Z[i]-zbar1[j])*phihat[i,j]*X[i]+(Z[i]-zbar1[j-1])*phihat[i,j-1]*X[i])*dallt[j]
    }
  }
  
  A21=sum(A3)/n
  
  
  #estimate A22
  A4=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4[i,j]=(Z[i]-zbar1[j])*phihat[i,j]*exp(-hbth[ct,2]*Z[i])*Z[i]*dN[i,j]
    }
  }
  
  A22=sum(A4)/n
  
  
  #estimate Ahat
  Ahat=matrix(c(A11,A12,A21,A22),2,2,byrow=T)
  
  Ainv=solve(Ahat)
  
  #estimate Sigmastar
  
  Sigmastar=Ainv%*%Sigmahat%*%Ainv
  
  sgth[1,ct]=sqrt(Sigmastar[1,1])/sqrt(n)
  
  sgth[2,ct]=sqrt(Sigmastar[2,2])/sqrt(n)
  
  ##counting CP
  
  
  
  
  
  cpth[1,ct]=ifelse((hbth[ct,1]-qnorm(1-0.05/2)*sgth[1,ct] < beta) & (beta < hbth[ct,1]+qnorm(1-0.05/2)*sgth[1,ct]),1,0)
  
  cpth[2,ct]=ifelse((hbth[ct,2]-qnorm(1-0.05/2)*sgth[2,ct] < theta) & (theta < hbth[ct,2]+qnorm(1-0.05/2)*sgth[2,ct]),1,0)
  
  
  Lambda=rep(0,nt)
  
  Lambda[1]=partlmd1[1]-partlmd2[1]
  
  for(i in 2:nt){
    Lambda[i]=partlmd1[i]-partlmd2[i]+Lambda[i-1]
  }
  
 
  
  #estimate LMD-sig-est
  
  k=apply(phihat,2,sum)/n
  
  #estimate u1
  
  u1part1=matrix(0,n,nt)
  
  
  
  for (i in 1:n){
    for (j in 1:nt){
      u1part1[i,j]=1/k[j]*phihat[i,j]*exp(-hbth[ct,2]*Z[i])*dN[i,j]
    }
  }
  
  u1part2=matrix(0,n,nt)
  for (i in 1:n){
    u1part2[i,1]=1/2*1/k[1]*phihat[i,1]*(hbth[ct,1]*X[i])*dallt[1]
    
    for (j in 2:nt){
      u1part2[i,j]=1/2*(1/k[j]*phihat[i,j]*(hbth[ct,1]*X[i])+1/k[j-1]*phihat[i,j-1]*(hbth[ct,1]*X[i]))*dallt[j]
    }
  }
  
  u1part3=matrix(0,n,nt)
  for (i in 1:n){
    
    for (j in 1:nt){
      u1part3[i,j]=1/k[j]*phihat[i,j]*partlmd1[j]
    }
  }
  
  
  
  
  
  u1part4=matrix(0,n,nt)
  for (i in 1:n){
    
    for (j in 1:nt){
      u1part4[i,j]=1/k[j]*phihat[i,j]*partlmd2[j]
    }
  }
  
  u1_0=matrix(0,n,nt)
  u1=matrix(0,n,nt)
  
  for(i in 1:n){
    for (j in 1:nt){
      u1_0[i,j]=u1part1[i,j]-u1part2[i,j]-u1part3[i,j]+u1part4[i,j]
    }
  }
  
  for(i in 1:n){
    u1[i,1]=u1_0[i,1]
       for (j in 2:nt){
         u1[i,j]=u1[i,j-1]+u1_0[i,j]
    }
  }
  
  #estimate Ghat
  
  sumG=apply(u1,2,sum)/n
  
  Ghat=matrix(0,nt,nt)
  for (i in 1:nt){
   for(j in i:nt){
     Ghat[i,j]=(sumG[j]-sumG[i])
   }
    }
  
  
  
  #estimate u2
  
  component8=matrix(0,nt,nt)
  
  
  for(j in 1:nt){
    for(i in 1:nt){
    component8[i,j]=Ghat[i,j]/R0[i]
  }
  }
  
  u2part1=matrix(0,n,nt)
  u2part1_0=array(0,c(n,nt,nt))
  u2part1_0sum=array(0,c(n,nt,nt))
  
  u2part3=matrix(0,n,nt)
  u2part3_0=array(0,c(n,nt,nt))
  u2part3_0sum=array(0,c(n,nt,nt))
  
  u2part4=matrix(0,n,nt)
  u2part4_0=array(0,c(n,nt,nt))
  u2part4_0sum=array(0,c(n,nt,nt))
  
  for (i in 1:n){
    
    for(s in 1:nt){
        for(j in 1:s){
      u2part1_0[i,j,s]=component8[j,s]*dNdd[i,j]
      u2part3_0[i,j,s]=component8[j,s]*DT3[i,j]*dDelta0hatpart1[j]
      u2part4_0[i,j,s]=component8[j,s]*DT3[i,j]*dDelta0hatpart2[j]
      
        }
      u2part1_0sum[i,1,s]=u2part1_0[i,1,s]
      u2part3_0sum[i,1,s]=u2part3_0[i,1,s]
      u2part4_0sum[i,1,s]=u2part4_0[i,1,s]
      
        if(s>1) {
          for(j in 2:s){
        u2part1_0sum[i,j,s]=u2part1_0sum[i,j-1,s]+u2part1_0[i,j,s]
        u2part3_0sum[i,j,s]=u2part3_0sum[i,j-1,s]+u2part3_0[i,j,s]
        u2part4_0sum[i,j,s]=u2part4_0sum[i,j-1,s]+u2part4_0[i,j,s]
        
         }
        }
      u2part1[i,s]=u2part1_0sum[i,s,s]
      u2part3[i,s]=u2part3_0sum[i,s,s]
      u2part4[i,s]=u2part4_0sum[i,s,s]
      
      }
     }
  

  
  u2part2=matrix(0,n,nt)
  u2part2_0=array(0,c(n,nt,nt))
  u2part2_0sum=array(0,c(n,nt,nt))
  
  for (i in 1:n){
    for(s in 1:nt){
    u2part2_0[i,1,s]=1/2*component8[1,s]*DT3[i,1]*hgm[ct]*Z[i]*dallt[1]
    for (j in 2:nt){
      u2part2_0[i,j,s]=1/2*(component8[j,s]*DT3[i,j]*hgm[ct]*Z[i]+component8[j-1,s]*DT3[i,j-1]*hgm[ct]*Z[i])*dallt[j]
    }
    }
  }
  
  for (i in 1:n){
    
    for(s in 1:nt){
      u2part2_0sum[i,1,s]=u2part2_0[i,1,s]
      
      if(s>1) {
        for(j in 2:s){
          u2part2_0sum[i,j,s]=u2part2_0sum[i,j-1,s]+u2part2_0[i,j,s]
        }
      }
      u2part2[i,s]=u2part2_0sum[i,s,s]
    }
  }
  
  

  
  
  
 
  u2=matrix(0,n,nt)
  

  for(i in 1:n){
    for (j in 1:nt){
      u2[i,j]=u2part1[i,j]-u2part2[i,j]-u2part3[i,j]+u2part4[i,j]
    }
  }
  
 
  
  #estimate Ohat
  
  Ohatpart1=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      Ohatpart1[i,j]=1/k[j]*(ETA[j]+allt[j]*Z[i])*((Omega)^(-1))*phihat[i,j]*exp(-hbth[ct,2]*Z[i])*dN[i,j]
    }
  }
  
  Ohatpart2=matrix(0,n,nt)
  for (i in 1:n){
    Ohatpart2[i,1]=1/2*1/k[1]*(ETA[1]+allt[1]*Z[i])*((Omega)^(-1))*phihat[i,1]*(hbth[ct,1]*X[i])*dallt[1]
    
    for (j in 2:nt){
      Ohatpart2[i,j]=1/2*(1/k[j]*(ETA[j]+allt[j]*Z[i])*((Omega)^(-1))*phihat[i,j]*(hbth[ct,1]*X[i])
                          +1/k[j-1]*(ETA[j-1]+allt[j-1]*Z[i])*((Omega)^(-1))*phihat[i,j-1]*(hbth[ct,1]*X[i]))*dallt[j]
    }
  }
  
  
  Ohatpart3=matrix(0,n,nt)
  
  
  for (i in 1:n){
    for (j in 1:nt){
      Ohatpart3[i,j]=1/k[j]*(ETA[j]+allt[j]*Z[i])*((Omega)^(-1))*phihat[i,j]*partlmd1[j]
    }
  }
  
  Ohatpart4=matrix(0,n,nt)
  
  for (i in 1:n){
    
    for (j in 1:nt){
      Ohatpart4[i,j]=1/k[j]*(ETA[j]+allt[j]*Z[i])*((Omega)^(-1))*phihat[i,j]*partlmd2[j]
    }
    
  }
  
  
  
  Ohat=rep(0,nt)
  
  Ohatpart11=apply(Ohatpart1, 2, sum)
  Ohatpart22=apply(Ohatpart2, 2, sum)
  Ohatpart33=apply(Ohatpart3, 2, sum)
  Ohatpart44=apply(Ohatpart4, 2, sum)
  
  
  
  Ohat[1]=Ohatpart11[1]-Ohatpart22[1]-Ohatpart33[1]+Ohatpart44[1]
  
  for (j in 2:nt){
    Ohat[j]=Ohat[j-1]+Ohatpart11[j]-Ohatpart22[j]-Ohatpart33[j]+Ohatpart44[j]
  }
  
  Ohat=Ohat/n
  
  
  #estimate u3
  
  u3part1=matrix(0,n,nt)
  
  for (i in 1:n){
    for (j in 1:nt){
      u3part1[i,j]=(Z[i]-wbar1[j])*dNdd[i,j]
    }
  }
  
  u3part2=matrix(0,n,nt)
  for (i in 1:n){
    u3part2[i,1]=1/2*(Z[i]-wbar1[1])*DT3[i,1]*hgm[ct]*Z[i]*dallt[1]
    
    for (j in 2:nt){
      u3part2[i,j]=1/2*((Z[i]-wbar1[j])*DT3[i,j]*hgm[ct]*Z[i]+(Z[i]-wbar1[j-1])*DT3[i,j-1]*hgm[ct]*Z[i])*dallt[j]
    }
  }
  
  u3part3=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      u3part3[i,j]=(Z[i]-wbar1[j])*DT3[i,j]*dDelta0hatpart1[j]
    }
  }
  
  u3part4=matrix(0,n,nt)
  for (i in 1:n){
    
    for (j in 1:nt){
      u3part4[i,j]=(Z[i]-wbar1[j])*DT3[i,j]*dDelta0hatpart2[j]
    }
  }
  
  
  u33=apply(u3part1-u3part2-u3part3+u3part4,1,sum)
  
  u3=matrix(0,n,nt)
  
  for( i in 1:n){
    for (j in 1:nt){
      u3[i,j]=Ohat[j]*u33[i]
    }
  }
  
  #estimate Gamma
  
  Gammapart1=matrix(0,n,nt)
  
  for (i in 1:n){
    
    Gammapart1[i,1]=1/2*phihat[i,1]/k[1]*X[i]*dallt[1]
    
    for (j in 2:nt){
      Gammapart1[i,j]=1/2*(phihat[i,j]/k[j]*X[i]+phihat[i,j-1]/k[j-1]*X[i])*dallt[j]
    }
  }
  
  Gammapart2=matrix(0,n,nt)
  
  for (i in 1:n){
    for (j in 1:nt){
      Gammapart2[i,j]=phihat[i,j]/k[j]*Z[i]*exp(-hbth[ct,2]*Z[i])*dN[i,j]
    }
  }
  
  Gammapart11_0=apply(Gammapart1,2,sum)/n
  

  Gammapart22_0=apply(Gammapart2,2,sum)/n
  
  Gammapart11=rep(0,nt)
  Gammapart11[1]=Gammapart11_0[1]
    for (j in 2:nt){
      Gammapart11[j]=Gammapart11[j-1]+Gammapart11_0[j]
    }
  
  Gammapart22=rep(0,nt)
  Gammapart22[1]=Gammapart22_0[1]
  for (j in 2:nt){
    Gammapart22[j]=Gammapart22[j-1]+Gammapart22_0[j]
  }
  
  
  Gamma=rbind(-Gammapart11,-Gammapart22)
  
  #estimate g_i
  
  sgs=matrix(0,n,nt)
  
  for (i in 1:n){
    for (j in 1:nt){
      sgs[i,j]=Gamma[,j]%*%Ainv%*%as.matrix(q1[,i]+q2[,i]+q3[,i])+u1[i,j]+u2[i,j]+u3[i,j]
    }
  }
  
  
  ## sglm, estlmds: defined on common time segments
  sglm <- numeric(NT0)
  estlmds <- numeric(NT0)
  
  ## sglm and estlmds are nonzero only when t>=min(allt)
  stt <- min(which(T0>=min(allt))) 
  for (j in stt:NT0) {
    b <- max(which(allt<=T0[j])) 
    sglm[j] <- sqrt((t(sgs[,b]) %*% sgs[,b])/n^2)
    estlmds[j] <- Lambda[b]
  }
  
  ######################### change time scale to allt
  Lmd_est[ct,] <- estlmds
  Lmd_sig_est[ct,] <- sglm
  
 

  
  mv=matrix(c(0.2,0.4,1,0,1,1),2,byrow=T)
  nv=ncol(mv)
  
  colns=c(round(nt/5),2*round(nt/5),3*round(nt/5),4*round(nt/5))
  ns=length(colns)
  
  DTv=matrix(0,n,nv)
  for(i in 1:n){
    for(j in 1:nv){
    DTv[i,j]=as.numeric(V[1,i]<=mv[1,j])*as.numeric(V[2,i]<=mv[2,j])
    }
  }
    

  
  Dhat=array(0,c(nt,nv,n))
  Dhatcomp=array(0,c(nt,nv,n))
  
  for(s in 1:nt){
    for (j in 1:nv){
      for (i in 1:n){
        Dhatcomp[s,j,i]=DTv[i,j]*phihat[i,s]
      }      }
  }
  
  Dhatcompsum=array(0,c(nt,nv,n))
 frac1=array(0,c(nt,nv))
  for(j in 1:nv){
    for (s in 1:nt){
      Dhatcompsum[s,j,1]=Dhatcomp[s,j,1]
        for (i in 2:n){
        Dhatcompsum[s,j,i]=Dhatcompsum[s,j,i-1]+Dhatcomp[s,j,i]
        }  
      frac1[s,j]=Dhatcompsum[s,j,n]
    }
  }
  
  frac2=apply(phihat,2,sum)
  
  for(s in 1:nt){
    for (j in 1:nv){
      for (i in 1:n){
        Dhat[s,j,i]=DTv[i,j]-frac1[s,j]/frac2[s]
      }      }
  }

  
  
  etahat1=array(0,c(nt,nv,n))
  etahat1part1=array(0,c(nt,nv,n))
  etahat1part2=array(0,c(nt,nv,n))
  etahat1part3=array(0,c(nt,nv,n))
  etahat1part4=array(0,c(nt,nv,n))
  

  
  for (i in 1:n)
  {
    for (j in 1:nv){
      for (s in 1:nt){
        etahat1part1[s,j,i]=Dhat[s,j,i]*phihat[i,s]*exp(-hbth[ct,2]*Z[i])*dN[i,s]
        etahat1part3[s,j,i]=Dhat[s,j,i]*phihat[i,s]*partlmd1[s]
        etahat1part4[s,j,i]=Dhat[s,j,i]*phihat[i,s]*partlmd2[s]
        
      }
      etahat1part2[1,j,i]=1/2*Dhat[s,j,i]*phihat[i,1]*(hbth[ct,1]*X[i])*dallt[1]
      for (s in 2:nt){
        etahat1part2[s,j,i]=1/2*(Dhat[s,j,i]*phihat[i,s]*(hbth[ct,1]*X[i])+Dhat[s-1,j,i]*phihat[i,s-1]*(hbth[ct,1]*X[i]))*dallt[s]
      }
    }
    }
    
   
etahat1_0=array(0,c(nt,nv,n))
    
    for(i in 1:n){
      for (j in 1:nv){
        for (s in 1:nt){
        etahat1_0[s,j,i]=etahat1part1[s,j,i]-etahat1part2[s,j,i]-etahat1part3[s,j,i]+etahat1part4[s,j,i]
        }
      }
    }
    
for(i in 1:n){
    for(j in 1:nv){
      etahat1[1,j,i]=etahat1_0[1,j,i]
      for (s in 2:nt){
        etahat1[s,j,i]=etahat1[s-1,j,i]+etahat1_0[s,j,i]
      }
    }
}

etahat11=array(0,c(ns,nv,n))  
for(i in 1:n){
  for (j in 1:nv){
   for(s in 1:ns){
     etahat11[s,j,i]=etahat1[colns[s],j,i]
   }
      }
  
}
  
#estimate etahat2

#estimate Phiuvs
  
sumetahat1_0=array(0,c(nt,nv,n))
sumetahat1=array(0,c(nt,nv))

for (j in 1:nv){
  for (s in 1:nt){
    sumetahat1_0[s,j,1]=etahat1[s,j,1]
    for (i in 2:n){
      sumetahat1_0[s,j,i]=sumetahat1_0[s,j,i-1]+etahat1[s,j,i]
    }
  }
}

for(s in 1:nt){
  for (j in 1:nv){
    sumetahat1[s,j]=sumetahat1_0[s,j,n]
  }
}


Phiuvs=array(0,c(nt,nv,nt))


for(j in 1:nv){
for (u in 1:nt){
  for(s in u:nt)
    {
    Phiuvs[u,j,s]=(sumetahat1[s,j]-sumetahat1[u,j])/n
  }
}
}


  




component9=array(0,c(nt,nv,nt))
for(u in 1:nt){
  for (j in 1:nv){
    for (s in 1:nt){
      component9[u,j,s]=Phiuvs[u,j,s]/R0[u]
    }
  }
}



etahat2part1=array(0,c(nt,nt,nv,n))
etahat2part2=array(0,c(nt,nt,nv,n))
etahat2part3=array(0,c(nt,nt,nv,n))
etahat2part4=array(0,c(nt,nt,nv,n))

for (i in 1:n){ 
  for(j in 1:nv){
  for (s in 1:nt){
    for(u in 1:nt){
    etahat2part1[u,s,j,i]=component9[u,j,s]*dNdd[i,u]
    etahat2part3[u,s,j,i]=component9[u,j,s]*DT3[i,u]*dDelta0hatpart1[u]
    etahat2part4[u,s,j,i]=component9[u,j,s]*DT3[i,u]*dDelta0hatpart2[u]
    }
  
    etahat2part2[1,s,j,i]=1/2*component9[1,j,s]*DT3[i,1]*hgm[ct]*Z[i]*dallt[1]
    for(u in 2:nt){
      etahat2part2[u,s,j,i]=1/2*(component9[u,j,s]*DT3[i,u]*hgm[ct]*Z[i]+component9[u-1,j,s]*DT3[i,u-1]*hgm[ct]*Z[i])*dallt[u]
    }
  }
    }
}


etahat2part1_0=array(0,c(nt,nv,n))
etahat2part2_0=array(0,c(nt,nv,n))
etahat2part3_0=array(0,c(nt,nv,n))
etahat2part4_0=array(0,c(nt,nv,n))

etahat2part11=array(0,c(nt,nt,nv,n))
etahat2part22=array(0,c(nt,nt,nv,n))
etahat2part33=array(0,c(nt,nt,nv,n))
etahat2part44=array(0,c(nt,nt,nv,n))


for(i in 1:n){
  for (j in 1:nv){
    for (s in 1:nt){
      etahat2part11[1,s,j,i]=etahat2part1[1,s,j,i]
      etahat2part22[1,s,j,i]=etahat2part2[1,s,j,i]
      etahat2part33[1,s,j,i]=etahat2part3[1,s,j,i]
      etahat2part44[1,s,j,i]=etahat2part4[1,s,j,i]
      
      if(s>1){
      for (u in 2:s){
        etahat2part11[u,s,j,i]=etahat2part1[u,s,j,i]+etahat2part11[u-1,s,j,i]
        etahat2part22[u,s,j,i]=etahat2part2[u,s,j,i]+etahat2part22[u-1,s,j,i]
        etahat2part33[u,s,j,i]=etahat2part3[u,s,j,i]+etahat2part33[u-1,s,j,i]
        etahat2part44[u,s,j,i]=etahat2part4[u,s,j,i]+etahat2part44[u-1,s,j,i]
        
        }
      }
      
      etahat2part1_0[s,j,i]=etahat2part11[s,s,j,i]
      etahat2part2_0[s,j,i]=etahat2part22[s,s,j,i]
      etahat2part3_0[s,j,i]=etahat2part33[s,s,j,i]
      etahat2part4_0[s,j,i]=etahat2part44[s,s,j,i]
      
      }
  }
}




etahat2=array(0,c(nt,nv,n))

for(i in 1:n){
  for (j in 1:nv){
    for (s in 1:nt){
      etahat2[s,j,i]=etahat2part1_0[s,j,i]-etahat2part2_0[s,j,i]-etahat2part3_0[s,j,i]+etahat2part4_0[s,j,i]
    }
  }
}

etahat22=array(0,c(ns,nv,n))  
for(i in 1:n){
  for (j in 1:nv){
    for(s in 1:ns){
      etahat22[s,j,i]=etahat2[colns[s],j,i]
    }
  }
  
}


#estimate Bsv

Bsvpart1=array(0,c(nt,nv,n))
Bsvpart2=array(0,c(nt,nv,n))
Bsvpart3=array(0,c(nt,nv,n))
Bsvpart4=array(0,c(nt,nv,n))

for (i in 1:n){
  for(j in 1:nv){
  for (s in 1:nt){
    Bsvpart1[s,j,i]=Dhat[s,j,i]*(ETA[s]+allt[s]*Z[i])*((Omega)^(-1))*phihat[i,s]*exp(-hbth[ct,2]*Z[i])*dN[i,s]
  
    Bsvpart3[s,j,i]=Dhat[s,j,i]*(ETA[s]+allt[s]*Z[i])*((Omega)^(-1))*phihat[i,s]*partlmd1[s]
    Bsvpart4[s,j,i]=Dhat[s,j,i]*(ETA[s]+allt[s]*Z[i])*((Omega)^(-1))*phihat[i,s]*partlmd2[s]
    
    }
    Bsvpart2[1,j,i]=1/2*Dhat[1,j,i]*(ETA[1]+allt[1]*Z[i])*((Omega)^(-1))*phihat[i,1]*(hbth[ct,1]*X[i])*dallt[1]
    
    for(s in 2:nt){
      Bsvpart2[s,j,i]=1/2*(Dhat[s,j,i]*(ETA[s]+allt[s]*Z[i])*((Omega)^(-1))*phihat[i,s]*(hbth[ct,1]*X[i])
                           +Dhat[s-1,j,i]*(ETA[s-1]+allt[s-1]*Z[i])*((Omega)^(-1))*phihat[i,s-1]*(hbth[ct,1]*X[i]))*dallt[s]
    }
  }
}



Bsv=array(0,c(nt,nv))
Bsv_0=array(0,c(nt,nv))

sumBsvpart1_0=array(0,c(nt,nv,n))
sumBsvpart1=array(0,c(nt,nv))
sumBsvpart2_0=array(0,c(nt,nv,n))
sumBsvpart2=array(0,c(nt,nv))
sumBsvpart3_0=array(0,c(nt,nv,n))
sumBsvpart3=array(0,c(nt,nv))
sumBsvpart4_0=array(0,c(nt,nv,n))
sumBsvpart4=array(0,c(nt,nv))

for(j in 1:nv){

  for (s in 1:nt){
    sumBsvpart1_0[s,j,1]=Bsvpart1[s,j,1]
    sumBsvpart2_0[s,j,1]=Bsvpart2[s,j,1]
    sumBsvpart3_0[s,j,1]=Bsvpart3[s,j,1]
    sumBsvpart4_0[s,j,1]=Bsvpart4[s,j,1]
  
    for (i in 2:n){
      sumBsvpart1_0[s,j,i]=sumBsvpart1_0[s,j,i-1]+Bsvpart1[s,j,i]
      sumBsvpart2_0[s,j,i]=sumBsvpart2_0[s,j,i-1]+Bsvpart2[s,j,i]
      sumBsvpart3_0[s,j,i]=sumBsvpart3_0[s,j,i-1]+Bsvpart3[s,j,i]
      sumBsvpart4_0[s,j,i]=sumBsvpart4_0[s,j,i-1]+Bsvpart4[s,j,i]
      
        }
  }

}

for(j in 1:nv){
  for (s in 1:nt){
    sumBsvpart1[s,j]=sumBsvpart1_0[s,j,n]
    sumBsvpart2[s,j]=sumBsvpart2_0[s,j,n]
    sumBsvpart3[s,j]=sumBsvpart3_0[s,j,n]
    sumBsvpart4[s,j]=sumBsvpart4_0[s,j,n]
    
    
  }
}

for (j in 1:nv){

for (s in 1:nt){
  Bsv_0[s,j]=sumBsvpart1[s,j]-sumBsvpart2[s,j]-sumBsvpart3[s,j]+sumBsvpart4[s,j]
}
}

Bsv_0=Bsv_0/n


  for(j in 1:nv){
    Bsv[1,j]=Bsv_0[1,j]
    for (s in 2:nt){
      Bsv[s,j]=Bsv[s-1,j]+Bsv_0[s,j]
    }
  }

Bsvv=array(0,c(ns,nv))  

  for (j in 1:nv){
    for(s in 1:ns){
     Bsvv[s,j]=Bsv[colns[s],j]
    }
  }
  


#estimate etahat3

etahat33=array(0,c(ns,nv,n))

for( i in 1:n){
  for (j in 1:nv){
  for (s in 1:ns){
    etahat33[s,j,i]=Bsvv[s,j]*u33[i]
  }
  }
}

#estimate Ksv


Ksvpart1=array(0,c(nt,nv,n))
Ksvpart2=array(0,c(nt,nv,n))
sumKsvpart1=array(0,c(nt,nv))
sumKsvpart2=array(0,c(nt,nv))
sumKsvpart1_0=array(0,c(nt,nv,n))
sumKsvpart2_0=array(0,c(nt,nv,n))


for (i in 1:n){
  for(j in 1:nv){
  Ksvpart1[1,j,i]=1/2*phihat[i,1]*Dhat[1,j,i]*X[i]*dallt[1]
  
  for (s in 2:nt){
    Ksvpart1[s,j,i]=1/2*(phihat[i,s]*Dhat[s,j,i]*X[i]+phihat[i,s-1]*Dhat[s-1,j,i]*X[i])*dallt[s]
  }
  for(s in 1:nt){
    Ksvpart2[s,j,i]=phihat[i,s]*Dhat[s,j,i]*Z[i]*exp(-hbth[ct,2]*Z[i])*dN[i,s]
    
  }
  }
}

for(j in 1:nv){

  for (s in 1:nt){
    sumKsvpart1_0[s,j,1]=Ksvpart1[s,j,1]
    sumKsvpart2_0[s,j,1]=Ksvpart2[s,j,1]
    
    for (i in 2:n){
      sumKsvpart1_0[s,j,i]=sumKsvpart1_0[s,j,i-1]+Ksvpart1[s,j,i]
      sumKsvpart2_0[s,j,i]=sumKsvpart2_0[s,j,i-1]+Ksvpart2[s,j,i]
      
    }
  }
}


for(j in 1:nv){
  for (s in 1:nt){
    sumKsvpart1[s,j]=sumKsvpart1_0[s,j,n]/n
    sumKsvpart2[s,j]=sumKsvpart2_0[s,j,n]/n
    }
}


sumKsvpart11=array(0,c(nt,nv))
sumKsvpart22=array(0,c(nt,nv))

for(j in 1:nv){
  sumKsvpart11[1,j]=sumKsvpart1[1,j]
  sumKsvpart22[1,j]=sumKsvpart2[1,j]
  
for(s in 2:nt){
  sumKsvpart11[s,j]=sumKsvpart11[s-1,j]+sumKsvpart1[s,j]
  sumKsvpart22[s,j]=sumKsvpart22[s-1,j]+sumKsvpart2[s,j]
  
  }}




Ksvv=array(0,c(2,ns,nv))

for(j in 1:nv){
  for (s in 1:ns){
    Ksvv[1,s,j]=-sumKsvpart11[colns[s],j]
    Ksvv[2,s,j]=-sumKsvpart22[colns[s],j]
    
  }
}

  


#estimating obsF



Fsvpart1=array(0,c(nt,nv,n))
Fsvpart2=array(0,c(nt,nv,n))
Fsvpart3=array(0,c(nt,nv,n))
Fsvpart4=array(0,c(nt,nv,n))

for (i in 1:n)
{
  for (j in 1:nv){
    for (s in 1:nt){
      Fsvpart1[s,j,i]=DTv[i,j]*phihat[i,s]*exp(-hbth[ct,2]*Z[i])*dN[i,s]
      Fsvpart3[s,j,i]=DTv[i,j]*phihat[i,s]*partlmd1[s]
      Fsvpart4[s,j,i]=DTv[i,j]*phihat[i,s]*partlmd2[s]
      
    }
    Fsvpart2[1,j,i]=1/2*DTv[i,j]*phihat[i,1]*(hbth[ct,1]*X[i])*dallt[1]
    for (s in 2:nt){
      Fsvpart2[s,j,i]=1/2*(DTv[i,j]*phihat[i,s]*(hbth[ct,1]*X[i])+DTv[i,j]*phihat[i,s-1]*(hbth[ct,1]*X[i]))*dallt[s]
    }
  }
  }


Fsv_0=array(0,c(nt,nv,n))
Fsv_1=array(0,c(nt,nv,n))
Fsv_2=array(0,c(nt,nv,n))
Fsv_3=array(0,c(nt,nv))
obsF=array(0,c(ns,nv))

for(i in 1:n){
  for (j in 1:nv){
    for (s in 1:nt){
      Fsv_0[s,j,i]=Fsvpart1[s,j,i]-Fsvpart2[s,j,i]-Fsvpart3[s,j,i]+Fsvpart4[s,j,i]
    }
  }
}

for(i in 1:n){
  for(j in 1:nv){
    Fsv_1[1,j,i]=Fsv_0[1,j,i]
    for (s in 2:nt){
      Fsv_1[s,j,i]=Fsv_1[s-1,j,i]+Fsv_0[s,j,i]
    }
  }
}

for(j in 1:nv){
  for(s in 1:nt){
    Fsv_2[s,j,1]=Fsv_1[s,j,1]
    for (i in 2:n){
      Fsv_2[s,j,i]=Fsv_2[s,j,i-1]+Fsv_1[s,j,i]
    }
    Fsv_3[s,j]=Fsv_2[s,j,n]/sqrt(n)
  }
}

for (j in 1:nv){
  for (s in 1:ns){
    obsF[s,j]=Fsv_3[colns[s],j]
  }
}



#within each replication, calculate max(abs(obsF)) and max(abs(calF))


  
  mf <- numeric(1000)
  for (ck in 1:1000) {
     gs=rnorm(n)
      # calculate calF[i] 
      Fsvhat_0=array(0,c(ns,nv,n))
    
    for (i in 1:n){
      
      for (j in 1:nv){
        for(s in 1:ns){
          Fsvhat_0[s,j,i]=(etahat11[s,j,i]+etahat22[s,j,i]+etahat33[s,j,i]+t(as.matrix(Ksvv[,s,j]))%*%Ainv%*%as.matrix(q1[,i]+q2[,i]+q3[,i]))*gs[i]
        }
      }
    }
    
    calF=array(0,c(ns,nv))
    Fsvhatsum=array(0,c(ns,nv,n))
    
    for(j in 1:nv){
      for (s in 1:ns){
        Fsvhatsum[s,j,1]=Fsvhat_0[s,j,1]
        for (i in 2:n){
          Fsvhatsum[s,j,i]=Fsvhat_0[s,j,i]+Fsvhatsum[s,j,i-1]
        }
        calF[s,j]=Fsvhatsum[s,j,n]/sqrt(n)
      }
    }
    
      mf[ck] <- max(abs(calF))
  }
  
  
  
  rej[ct] <- as.numeric(mean(mf>= max(abs(obsF)))<=0.05)
  
 
  



ct=ct+1
  
  
}


# if model is correct, mean(rej) is the empirical test size and should be close to 0.05
mean(rej)

save.image(file="0826_1_200.RData")

