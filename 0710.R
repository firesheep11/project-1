
library(BB)
library(stats)
library(elliptic)
library(nleqslv)

task=read.csv(file.choose())

task<-subset(task, task$RATRIAL!='' )
attach(task)
n=nrow(task)

try2<-subset(task, task$RVISINF=="Y" )
attach(try2)
n=nrow(try2)

try3<-subset(task, task$AGE>85)
attach(try3)
n=nrow(try3)

try4<-subset(task, task$RCONSC=="U")
attach(try4)
n=nrow(try4)

try5<-subset(task, task$FDENNIS=="Y")
attach(try5)
n=nrow(try5)

try6<-subset(task, task$STYPE=="POCS")
attach(try6)
n=nrow(try6)

try7<-subset(task, task$RDELAY>29)
attach(try7)
n=nrow(try7)
#data size

#additive effects X:age
X1=rep(0,n)
X1[AGE<50]=1
X1[50<=AGE&AGE<59]=2
X1[59<=AGE&AGE<69]=3
X1[69<=AGE&AGE<79]=4
X1[AGE>=79]=5

X1=AGE

#or
X2=RDELAY

#or
X3=RSBP

X=X1

#multiplicative effects Z:treatment group
group1=rep(0,n)
group1[RXASP=="Y"]=1

group1[RXASP=="N"]=0


group2=rep(0,n)

group2[(RXHEP=="M"|RXHEP=="H"|RXHEP=="L")]=1
group2[RXHEP=="N"]=0


Z=group
#or
Z1=group


RATRIAL[RATRIAL=="N"]=0
RATRIAL[RATRIAL=="Y"]=1

RATRIAL=as.numeric(RATRIAL)

RVISINF[RVISINF=="N"]=0
RVISINF[RVISINF=="Y"]=1

RVISINF=as.numeric(RVISINF)

X=X1
Z1=group1
Z2=group2
Z3=RATRIAL
Z4=RVISINF
Z5=X3

Z=cbind(Z1,Z2,Z3,Z4,Z5)




#event of interest: the first stroke recurrence (no matter which type)
T0=cbind(DRSISCD,DRSHD,DRSUNKD)
T0[is.na(T0)]=1400
C=apply(T0,1,min)

#C=DALIVED
#C[is.na(C)]=580

D=FDEADD
D[is.na(D)]=580


C1=rep(0,n)
C1[DIED==0]=TD[DIED==0]
C1[C1==0]=580
C1[is.na(C1)]=580

#DALIVE indicates patient discharge alive within 14 days.

FDEADD[DALIVE== "Y"] 
#still have death time for patients who discharged alive. 
#so the discharge is not work as censoring here.





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

fit=cox.aalen(Surv(TT,dtCD)~prop(X)+prop(Z1)+prop(Z2)+prop(Z3)+prop(Z4)+prop(Z5),n.sim=400,max.time=10)


summary(fit)
### survival data
# dt1: delta2
dt1 <- as.numeric(D<=C1)
# TT1: T2
TT1 <- apply(cbind(C1,D), 1, min)
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



W=cbind(X,Z)
mm=ncol(W)
mws <- array(0, dim=c(n, ntd, mm))
mwbar=array(0,dim=c(ntd,mm))
for(k in 1:mm){
  mws[,,k]=matrix(rep(W[,k],ntd),n,ntd)
  mwbar[,k]=apply(DT1 * mws[,,k], 2, sum) / apply(DT1*1, 2,sum)
  
}



#estimate gamma
leftgm=matrix(0,n,ntd)
rightgm=matrix(0,n,ntd)
fgamma=function(gamma){
  scores=numeric(mm)
  for(k in 1:mm){
    for(i in 1:n){
      for (j in 1:ntd){
        leftgm[i,j]=DT1[i,j]*(W[i,k]-mwbar[j,k])*dNd[i,j]
      }
    }
    leftgm1=sum(leftgm)
    for(i in 1:n){
      rightgm[i,1]=1/2*DT1[i,1]*(W[i,k]-mwbar[1,k])*gamma[k]*W[i,k]*dallt1[1]
      for (j in 2:ntd){
        
        rightgm[i,j]=1/2*(DT1[i,j-1]*(W[i,k]-mwbar[j-1,k])*gamma[k]*Z[i]+DT1[i,j]*(W[i,k]-mwbar[j,k])*gamma[k]*W[i,k])*dallt1[j]
        
      }
    }
    rightgm1=sum(rightgm)
    scores[k]=rightgm1-leftgm1
  }
  
  scores
}


#hgm[ct]=nleqslv(0, fgamma,control=list(trace=1,btol=.01,delta="newton"))$x

hgm=nleqslv(numeric(mm), fgamma,control=list(trace=1,btol=.01,delta="newton"))$x




##estimate shat estimate Deltahat with allt,length=nt
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



DT3 <- as.matrix(matrix(rep(allt, n), n, nt, byrow=TRUE)<= matrix(rep(TT1, nt), n, nt))





shat=matrix(0,n,nt)

part0=apply(DT3*dNdd,2,sum)/ apply(DT3, 2, sum)
part1=rep(0,nt)
part1[1]=part0[1]
for(j in 2:nt){
  part1[j]=part0[j]+part1[j-1]
}



component5=matrix(0,n,nt)
for(i in 1:n){
  for(j in 1:nt){
    component5[i,j]=DT3[i,j]*((hgm)%*%W[i,])
  }
}
component6=apply(component5,2,sum)/ apply(DT3, 2, sum)


part2=rep(0,nt)

part2[1]=1/2*component6[1]*dallt[1]
for(j in 2:nt){
  part2[j]=1/2*(component6[j]+component6[j-1])*dallt[j]+part2[j-1]
}




for(i in 1:n){
  for(j in 1:nt){
    shat[i,j]=exp(-part1[j]+part2[j]-((hgm)%*%W[i,])*allt[j])
  }
}



#estimate phihat

m1 <- matrix(rep(1, nt), n, nt)

phihat <- DT * m1/shat

mx <- matrix(rep(X, nt), n, nt)


xbar <- apply(phihat * mx, 2, sum) / apply(phihat, 2, sum)



pro1=matrix(0,n,nt)
pro2=rep(0,nt)
for(j in 1:nt){
  pro1[1,j]=Z[1,1]*phihat[1,j]
  for(i in 2:n){
    pro1[i,j]=pro1[i-1,j]+Z[i,1]*phihat[i,j]
  }
  pro2[j]=pro1[n,j]
}


pro3=matrix(0,n,nt)
pro4=rep(0,nt)
for(j in 1:nt){
  pro3[1,j]=Z[1,2]*phihat[1,j]
  for(i in 2:n){
    pro3[i,j]=pro3[i-1,j]+Z[i,2]*phihat[i,j]
  }
  pro4[j]=pro3[n,j]
}


pro5=matrix(0,n,nt)
pro6=rep(0,nt)
for(j in 1:nt){
  pro5[1,j]=Z[1,3]*phihat[1,j]
  for(i in 2:n){
    pro5[i,j]=pro5[i-1,j]+Z[i,3]*phihat[i,j]
  }
  pro6[j]=pro5[n,j]
}

pro7=matrix(0,n,nt)
pro8=rep(0,nt)
for(j in 1:nt){
  pro7[1,j]=Z[1,4]*phihat[1,j]
  for(i in 2:n){
    pro7[i,j]=pro7[i-1,j]+Z[i,4]*phihat[i,j]
  }
  pro8[j]=pro7[n,j]
}


pro9=matrix(0,n,nt)
pro10=rep(0,nt)
for(j in 1:nt){
  pro9[1,j]=Z[1,5]*phihat[1,j]
  for(i in 2:n){
    pro9[i,j]=pro9[i-1,j]+Z[i,5]*phihat[i,j]
  }
  pro10[j]=pro9[n,j]
}

pro11=matrix(0,n,nt)
pro12=rep(0,nt)
for(j in 1:nt){
  pro11[1,j]=Z[1,6]*phihat[1,j]
  for(i in 2:n){
    pro11[i,j]=pro11[i-1,j]+Z[i,6]*phihat[i,j]
  }
  pro12[j]=pro11[n,j]
}

zbar <- cbind(pro2,pro4,pro6,pro8,pro10)/apply(phihat, 2, sum)




V=t(cbind(X,Z))
vbar=t(cbind(xbar,zbar))



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
      leftth1[i,j]=(V[1,i]-vbar[1,j])*phihat[i,j]*exp(-(beth[2]*Z[i,1]+beth[3]*Z[i,2]+beth[4]*Z[i,3]+beth[5]*Z[i,4]+beth[6]*Z[i,5]))*dN[i,j]
    }
  }
  leftth11=sum(leftth1)
  
  rightth1=matrix(0,n,nt)
  for(i in 1:n){
    rightth1[i,1]=1/2*(V[1,i]-vbar[1,1])*phihat[i,1]*(beth[1]*X[i])*dallt[1]
    for (j in 2:nt){
      rightth1[i,j]=1/2*((V[1,i]-vbar[1,j])*phihat[i,j]*(beth[1]*X[i])+(V[1,i]-vbar[1,j-1])*phihat[i,j-1]*(beth[1]*X[i]))*dallt[j]
    }
  }
  rightth11=sum(rightth1)
  
  return(leftth11-rightth11)
}

fth2=function(beth){
  
  leftth2=matrix(0,n,nt)
  for(i in 1:n){
    for (j in 1:nt){
      leftth2[i,j]=(V[2,i]-vbar[2,j])*phihat[i,j]*exp(-(beth[2]*Z[i,1]+beth[3]*Z[i,2]+beth[4]*Z[i,3]+beth[5]*Z[i,4]+beth[6]*Z[i,5]))*dN[i,j]
    }
  }
  leftth22=sum(leftth2)
  
  rightth2=matrix(0,n,nt)
  for(i in 1:n){
    rightth2[i,1]=1/2*(V[2,i]-vbar[2,1])*phihat[i,1]*(beth[1]*X[i])*dallt[1]
    for (j in 2:nt){
      rightth2[i,j]=1/2*((V[2,i]-vbar[2,j])*phihat[i,j]*(beth[1]*X[i])+(V[2,i]-vbar[2,j-1])*phihat[i,j-1]*(beth[1]*X[i]))*dallt[j]
    }
  }
  rightth22=sum(rightth2)
  
  return(leftth22-rightth22)
}


fth3=function(beth){
  
  leftth3=matrix(0,n,nt)
  for(i in 1:n){
    for (j in 1:nt){
      leftth3[i,j]=(V[3,i]-vbar[3,j])*phihat[i,j]*exp(-(beth[2]*Z[i,1]+beth[3]*Z[i,2]+beth[4]*Z[i,3]+beth[5]*Z[i,4]+beth[6]*Z[i,5]))*dN[i,j]
    }
  }
  leftth33=sum(leftth3)
  
  rightth3=matrix(0,n,nt)
  for(i in 1:n){
    rightth3[i,1]=1/2*(V[3,i]-vbar[3,1])*phihat[i,1]*(beth[1]*X[i])*dallt[1]
    for (j in 2:nt){
      rightth3[i,j]=1/2*((V[3,i]-vbar[3,j])*phihat[i,j]*(beth[1]*X[i])+(V[3,i]-vbar[3,j-1])*phihat[i,j-1]*(beth[1]*X[i]))*dallt[j]
    }
  }
  rightth33=sum(rightth3)
  
  return(leftth33-rightth33)
}

fth4=function(beth){
  
  leftth4=matrix(0,n,nt)
  for(i in 1:n){
    for (j in 1:nt){
      leftth4[i,j]=(V[4,i]-vbar[4,j])*phihat[i,j]*exp(-(beth[2]*Z[i,1]+beth[3]*Z[i,2]+beth[4]*Z[i,3]+beth[5]*Z[i,4]+beth[6]*Z[i,5]))*dN[i,j]
    }
  }
  leftth44=sum(leftth4)
  
  rightth4=matrix(0,n,nt)
  for(i in 1:n){
    rightth4[i,1]=1/2*(V[4,i]-vbar[4,1])*phihat[i,1]*(beth[1]*X[i])*dallt[1]
    for (j in 2:nt){
      rightth4[i,j]=1/2*((V[4,i]-vbar[4,j])*phihat[i,j]*(beth[1]*X[i])+(V[4,i]-vbar[4,j-1])*phihat[i,j-1]*(beth[1]*X[i]))*dallt[j]
    }
  }
  rightth44=sum(rightth4)
  
  return(leftth44-rightth44)
}


fth5=function(beth){
  
  leftth5=matrix(0,n,nt)
  for(i in 1:n){
    for (j in 1:nt){
      leftth5[i,j]=(V[5,i]-vbar[5,j])*phihat[i,j]*exp(-(beth[2]*Z[i,1]+beth[3]*Z[i,2]+beth[4]*Z[i,3]+beth[5]*Z[i,4]+beth[6]*Z[i,5]))*dN[i,j]
    }
  }
  leftth55=sum(leftth5)
  
  rightth5=matrix(0,n,nt)
  for(i in 1:n){
    rightth5[i,1]=1/2*(V[5,i]-vbar[5,1])*phihat[i,1]*(beth[1]*X[i])*dallt[1]
    for (j in 2:nt){
      rightth5[i,j]=1/2*((V[5,i]-vbar[5,j])*phihat[i,j]*(beth[1]*X[i])+(V[5,i]-vbar[5,j-1])*phihat[i,j-1]*(beth[1]*X[i]))*dallt[j]
    }
  }
  rightth55=sum(rightth5)
  
  return(leftth55-rightth55)
}


fth6=function(beth){
  
  leftth6=matrix(0,n,nt)
  for(i in 1:n){
    for (j in 1:nt){
      leftth6[i,j]=(V[6,i]-vbar[6,j])*phihat[i,j]*exp(-(beth[2]*Z[i,1]+beth[3]*Z[i,2]+beth[4]*Z[i,3]+beth[5]*Z[i,4]+beth[6]*Z[i,5]))*dN[i,j]
    }
  }
  leftth66=sum(leftth6)
  
  rightth6=matrix(0,n,nt)
  for(i in 1:n){
    rightth6[i,1]=1/2*(V[6,i]-vbar[6,1])*phihat[i,1]*(beth[1]*X[i])*dallt[1]
    for (j in 2:nt){
      rightth6[i,j]=1/2*((V[6,i]-vbar[6,j])*phihat[i,j]*(beth[1]*X[i])+(V[6,i]-vbar[6,j-1])*phihat[i,j-1]*(beth[1]*X[i]))*dallt[j]
    }
  }
  rightth66=sum(rightth6)
  
  return(leftth66-rightth66)
}





fgm_0=function(beth){
  y=numeric(6)
  y[1]=fth1(beth)
  y[2]=fth2(beth)
  y[3]=fth3(beth)
  y[4]=fth4(beth)
  y[5]=fth5(beth)
  y[6]=fth6(beth)
  
  return(y)
}



hbth=nleqslv(c(0,0,0,0,0,0), fgm_0,control=list(trace=1,btol=.01,delta="newton"))$x



iter=1


while(iter<=1){
  partlmd1=rep(0,nt)
  partlmd11=matrix(0,n,nt)
  
  for (i in 1:n){
    for (j in 1:nt){
      partlmd11[i,j]=phihat[i,j]*exp(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5])*dN[i,j]
    }
    partlmd1=apply(partlmd11,2,sum)/apply(phihat, 2, sum)
  }
  
  partlmd2=rep(0,nt)
  
  partlmd22=matrix(0,n,nt)
  
  
  for (i in 1:n){
    partlmd22[i,1]=1/2*phihat[i,1]*hbth[1]*X[i]*dallt[1]
    
    for (j in 2:nt)
      partlmd22[i,j]=1/2*(phihat[i,j]*hbth[1]*X[i]+phihat[i,j-1]*hbth[1]*X[i])*dallt[j]
  }
  
  partlmd2=apply(partlmd22, 2, sum)/apply(phihat,2,sum)
  
  
  
  
  ###########estimate SEE based on theorems##################
  
  xbar1 <- xbar
  zbar1 <- zbar
  vbar1 <- vbar
  
  
  mws=array(0,dim=c(n,nt,mm))
  mwbar1=array(0,dim=c(nt,mm))
  for(k in 1:mm){
    mws[,,k]=matrix(rep(W[,k],nt),n,nt)
    mwbar1[,k]=apply(DT3 * mws[,,k], 2, sum) / apply(DT3*1, 2,sum)
    
  }
  
  
  
  
  
  #estimate q1
  
  q1part1=array(0,dim=c(n,nt,mm))
  q1part2=array(0,dim=c(n,nt,mm))
  q1part3=array(0,dim=c(n,nt,mm))
  q1part4=array(0,dim=c(n,nt,mm))
  
  for(k in 1:mm){
    for (i in 1:n){
      q1part2[i,1,k]=1/2*(V[k,i]-vbar1[k,1])*phihat[i,1]*(hbth[1]*X[i])*dallt[1]
      
      for (j in 1:nt){
        q1part1[i,j,k]=(V[k,i]-vbar1[k,j])*phihat[i,j]*exp(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5])*dN[i,j]
        
        q1part3[i,j,k]=(V[k,i]-vbar1[k,j])*phihat[i,j]*partlmd1[j]
        q1part4[i,j,k]=(V[k,i]-vbar1[k,j])*phihat[i,j]*partlmd2[j]
        
      }
      
      for (j in 2:nt){
        q1part2[i,j,k]=1/2*((V[k,i]-vbar1[k,j])*phihat[i,j]*(hbth[1]*X[i])+(V[k,i]-vbar1[k,j-1])*phihat[i,j-1]*(hbth[1]*X[i]))*dallt[j]
      }
    }
  }
  
  
  q11=matrix(0,n,mm)
  for(k in 1:mm){
    q11[,k]=apply(as.matrix(q1part1[,,k])-as.matrix(q1part2[,,k])-as.matrix(q1part3[,,k])+as.matrix(q1part1[,,k]),1,sum)
  }
  
  q1=t(q11)
  #estimate q2
  #estimate PHI
  
  
  
  PHIpart1=array(0,c(n,nt,mm))
  PHIpart2=array(0,c(n,nt,mm))
  PHIpart3=array(0,c(n,nt,mm))
  PHIpart4=array(0,c(n,nt,mm))
  
  for(k in 1:mm){
    sumPHIpart1=apply(as.matrix(q1part1[,,k]),1,sum)
    minus1=matrix(0,n,nt)
    for (i in 1:n){
      minus1[i,1]=q1part1[i,1,k]
      PHIpart1[i,1,k]=sumPHIpart1[i]-minus1[i,1]
      for (j in 2:nt){
        minus1[i,j]=minus1[i,j-1]+q1part1[i,j,k]
        PHIpart1[i,j,k]=sumPHIpart1[i]-minus1[i,j]
      }
    }
    
    sumPHIpart2=apply(as.matrix(q1part2[,,k]),1,sum)
    
    minus2=matrix(0,n,nt)
    for (i in 1:n){
      minus2[i,1]=q1part2[i,1,k]
      PHIpart2[i,1,k]=sumPHIpart2[i]-minus2[i,1]
      for (j in 2:nt){
        minus2[i,j]=minus2[i,j-1]+q1part2[i,j,k]
        PHIpart2[i,j,k]=sumPHIpart2[i]-minus2[i,j]
      }
    }
    
    sumPHIpart3=apply(as.matrix(q1part3[,,k]),1,sum)
    
    minus3=matrix(0,n,nt)
    for (i in 1:n){
      minus3[i,1]=q1part3[i,1,k]
      PHIpart3[i,1,k]=sumPHIpart3[i]-minus3[i,1]
      for (j in 2:nt){
        minus3[i,j]=minus3[i,j-1]+q1part3[i,j,k]
        PHIpart3[i,j,k]=sumPHIpart3[i]-minus3[i,j]
      }
    }
    
    sumPHIpart4=apply(as.matrix(q1part4[,,k]),1,sum)
    
    minus4=matrix(0,n,nt)
    for (i in 1:n){
      minus4[i,1]=q1part4[i,1,k]
      PHIpart4[i,1,k]=sumPHIpart4[i]-minus4[i,1]
      for (j in 2:nt){
        minus4[i,j]=minus4[i,j-1]+q1part4[i,j,k]
        PHIpart4[i,j,k]=sumPHIpart4[i]-minus4[i,j]
      }
    }
  }
  
  
  
  #PHihat also should have two dimension (x,z)
  PHIhat=matrix(0,nt,mm)
  for(k in 1:mm){
    PHIhat[,k]=apply(as.matrix(PHIpart1[,,k])-as.matrix(PHIpart2[,,k])-as.matrix(PHIpart3[,,k])+as.matrix(PHIpart4[,,k]),2,sum)/n
  }
  
  
  #estimate q2
  #estimate R0
  R0 <- apply(DT3 * m1, 2, sum) / n
  
  component5=matrix(0,mm,nt)
  
  for(k in 1:mm){
    for(j in 1:nt){
      component5[k,j]=PHIhat[j,k]/R0[j]
    }
  }
  
  q2part1=array(0,c(n,nt,mm))
  q2part2=array(0,c(n,nt,mm))
  q2part3=array(0,c(n,nt,mm))
  q2part4=array(0,c(n,nt,mm))
  
  dDelta0hatpart1=apply(dNdd,2,sum)/apply(DT3,2,sum)
  
  dDelta0hatpart22=matrix(0,n,nt)
  for(i in 1:n){
    dDelta0hatpart22[i,1]=1/2*DT3[i,j]*(hgm%*%W[i,])*dallt[1]
    
    for (j in 2:nt){
      dDelta0hatpart22[i,j]=1/2*(DT3[i,j]*(hgm%*%W[i,])+DT3[i,j-1]*(hgm%*%W[i,]))*dallt[j]
    }
  }
  
  dDelta0hatpart2=apply(dDelta0hatpart22,2,sum)/apply(DT3,2,sum)
  
  
  for(k in 1:mm){
    
    for (i in 1:n){
      q2part2[i,1,k]=1/2*component5[k,1]*DT3[i,1]*(hgm%*%W[i,])*dallt[1]
      
      for (j in 1:nt){
        q2part1[i,j,k]=component5[k,j]*dNdd[i,j]
        q2part3[i,j,k]=component5[k,j]*DT3[i,j]*dDelta0hatpart1[j]
        q2part4[i,j,k]=component5[k,j]*DT3[i,j]*dDelta0hatpart2[j]
        
      }
      
      
      for (j in 2:nt){
        q2part2[i,j,k]=1/2*(component5[k,j]*DT3[i,j]*(hgm%*%W[i,])+component5[k,j-1]*DT3[i,j-1]*(hgm%*%W[i,]))*dallt[j]
      }
      
    }
  }
  
  
  
  q2=matrix(0,n,mm)
  for(k in 1:mm){
    q2[,k]=apply(as.matrix(q2part1[,,k])-as.matrix(q2part2[,,k])-as.matrix(q2part3[,,k])+as.matrix(q2part4[,,k]),1,sum)
  }
  
  q2=t(q2)
  
  
  #estimate Eta
  
  
  
  mEta=array(0,dim=c(nt,mm))
  
  for(k in 1:mm){
    mEta[1,k]=1/2*mwbar1[1,k]*dallt[1]
    for (j in 2:nt){
      mEta[j,k]=1/2*(mwbar1[j,k]+mwbar1[j-1,k])*dallt[j]+mEta[j-1,k]
    }
  }
  
  mETA=-mEta
  
  
  #estimate Omega
  

  
  product=array(0,dim=c(n,nt,mm^2))  
  
  for(i in 1:n){
    
    for(k in 1:mm^2){
      
    if(k%%mm!=0){
    product[i,1,k]=as.matrix(1/2*(W[i,]-mwbar1[1,])%*%t(W[i,]-mwbar1[1,])*DT3[i,1]*dallt[1])[ceiling(k/mm),k%%mm]
    }
    else{    product[i,1,k]=as.matrix(1/2*(W[i,]-mwbar1[1,])%*%t(W[i,]-mwbar1[1,])*DT3[i,1]*dallt[1])[ceiling(k/mm),k%%mm+mm]
}
    
    for(j in 2:nt){
      if(k%%mm!=0){
      product[i,j,k]=as.matrix(1/2*((W[i,]-mwbar1[j,])%*%t(W[i,]-mwbar1[j,])*DT3[i,j]+(W[i,]-mwbar1[j-1,])%*%t(W[i,]-mwbar1[j-1,])*DT3[i,j-1])*dallt[j])[ceiling(k/mm),k%%mm]
     
      }
      else{product[i,j,k]=as.matrix(1/2*((W[i,]-mwbar1[j,])%*%t(W[i,]-mwbar1[j,])*DT3[i,j]+(W[i,]-mwbar1[j-1,])%*%t(W[i,]-mwbar1[j-1,])*DT3[i,j-1])*dallt[j])[ceiling(k/mm),k%%mm+mm]
}
    
  }
    }
    }
  
  
  Omega=matrix(0,mm,mm)
  for(k in 1:mm^2){
    
  if(k%%mm!=0){
  Omega[ceiling(k/mm),k%%mm]=sum(product[,,k])/n
  }
  else{Omega[ceiling(k/mm),k%%mm+mm]=sum(product[,,k])/n}

}
  
  
  Omeinv=solve(Omega)
  
  
  #estimate B
  Bpart1=array(0,dim=c(n,nt,mm^2))
  for (i in 1:n){
    for (j in 1:nt){
      for(k in 1:mm^2){
        if(k%%mm!=0){
          Bpart1[i,j,k]=as.matrix((V[,i]-vbar1[,j])%*%t(mETA[j,]+allt[j]*W[i,])%*%Omeinv*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*dN[i,j])[ceiling(k/mm),k%%mm]
          
        }
        else{ Bpart1[i,j,k]=as.matrix((V[,i]-vbar1[,j])%*%t(mETA[j,]+allt[j]*W[i,])%*%Omeinv*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*dN[i,j])[ceiling(k/mm),k%%mm+mm]
}
      }

    }
  }
  
  Bpart2=array(0,dim=c(n,nt,mm^2))
  for (i in 1:n){
    for(k in 1:mm^2){
      if(k%%mm!=0){
        Bpart2[i,1,k]=as.matrix(1/2*(V[,i]-vbar1[,1])%*%t(mETA[1,]+allt[1]*W[i,])%*%Omeinv*phihat[i,1]*(hbth[1]*X[i])*dallt[1])[ceiling(k/mm),k%%mm]
        }
      else{        Bpart2[i,1,k]=as.matrix(1/2*(V[,i]-vbar1[,1])%*%t(mETA[1,]+allt[1]*W[i,])%*%Omeinv*phihat[i,1]*(hbth[1]*X[i])*dallt[1])[ceiling(k/mm),k%%mm+mm]
}
    }

    for (j in 2:nt){
      for(k in 1:mm^2){
        if(k%%mm!=0){
          Bpart2[i,j,1]=as.matrix(1/2*((V[,i]-vbar1[,j])%*%t(mETA[j,]+allt[j]*W[i,])%*%Omeinv*phihat[i,j]*(hbth[1]*X[i])
                                       +(V[,i]-vbar1[,j-1])%*%t(mETA[j,]+allt[j]*W[i,])%*%Omeinv*phihat[i,j-1]*(hbth[1]*X[i]))*dallt[j])[ceiling(k/mm),k%%mm]
        }
        else{Bpart2[i,j,1]=as.matrix(1/2*((V[,i]-vbar1[,j])%*%t(mETA[j,]+allt[j]*W[i,])%*%Omeinv*phihat[i,j]*(hbth[1]*X[i])
                                          +(V[,i]-vbar1[,j-1])%*%t(mETA[j,]+allt[j]*W[i,])%*%Omeinv*phihat[i,j-1]*(hbth[1]*X[i]))*dallt[j])[ceiling(k/mm),k%%mm+mm]
        }
      }
     
    }
  }
  
  
  Bpart3=array(0,dim=c(n,nt,mm^2))
  
  
  for (i in 1:n){
    for (j in 1:nt){
      for(k in 1:mm^2){
        if(k%%mm!=0){
      Bpart3[i,j,k]=as.matrix((V[,i]-vbar1[,j])%*%t(mETA[j,]+allt[j]*W[i,])%*%Omeinv*phihat[i,j]*partlmd1[j])[ceiling(k/mm),k%%mm]
        }
        else{      Bpart3[i,j,k]=as.matrix((V[,i]-vbar1[,j])%*%t(mETA[j,]+allt[j]*W[i,])%*%Omeinv*phihat[i,j]*partlmd1[j])[ceiling(k/mm),k%%mm+mm]
}
    }
  }}
  
  Bpart4=array(0,dim=c(n,nt,mm^2))
  
  for (i in 1:n){
    for (j in 1:nt){
      for(k in 1:mm^2){
        if(k%%mm!=0){
          Bpart4[i,j,k]=as.matrix((V[,i]-vbar1[,j])%*%t(mETA[j,]+allt[j]*W[i,])%*%Omeinv*phihat[i,j]*partlmd2[j])[ceiling(k/mm),k%%mm]
        }
        else{      Bpart4[i,j,k]=as.matrix((V[,i]-vbar1[,j])%*%t(mETA[j,]+allt[j]*W[i,])%*%Omeinv*phihat[i,j]*partlmd2[j])[ceiling(k/mm),k%%mm+mm]
        }
      }
    }}
  
  
  
  B1=Bpart1-Bpart2-Bpart3+Bpart4
  Bhat=matrix(0,mm,mm)

  for(k in 1:mm^2){
    
    if(k%%mm!=0){
      Bhat[ceiling(k/mm),k%%mm]=sum(B1[,,k])/n
    }
    else{    Bhat[ceiling(k/mm),k%%mm+mm]=sum(B1[,,k])/n
}
    
  }
  
  
  #then estimate q3
  
  q3part1=array(0,dim=c(n,nt,mm))
  q3part2=array(0,dim=c(n,nt,mm))
  q3part3=array(0,dim=c(n,nt,mm))
  q3part4=array(0,dim=c(n,nt,mm))
  
  for(k in 1:mm){
    for (i in 1:n){
      q3part2[i,1,k]=as.matrix(1/2*Bhat%*%(W[i,]-mwbar1[1,])*DT3[i,1]*as.numeric(hgm%*%W[i,])*dallt[1])[k,1]
      
      for (j in 1:nt){
        q3part1[i,j,k]=as.matrix(Bhat%*%(W[i,]-mwbar1[j,])*dNdd[i,j])[k,1]
        q3part3[i,j,k]=as.matrix(Bhat%*%(W[i,]-mwbar1[j,])*DT3[i,j]*dDelta0hatpart1[j])[k,1]
        q3part4[i,j,k]=as.matrix(Bhat%*%(W[i,]-mwbar1[j,])*DT3[i,j]*dDelta0hatpart2[j])[k,1]
        
      }
      
      for (j in 2:nt){
        q3part2[i,j,k]=as.matrix(1/2*(Bhat%*%(W[i,]-mwbar1[j,])*DT3[i,j]*as.numeric(hgm%*%W[i,])+Bhat%*%(W[i,]-mwbar1[j-1,])*DT3[i,j-1]*as.numeric(hgm%*%W[i,]))*dallt[j])[k,1]
        
      }
      
    }
  }
  
  
  
  q31=array(0,dim=c(n,nt))
  q32=array(0,dim=c(n,nt))
  q33=array(0,dim=c(n,nt))
  q34=array(0,dim=c(n,nt))
  q35=array(0,dim=c(n,nt))
  q36=array(0,dim=c(n,nt))
  
  for(i in 1:n){
    for(j in 1:nt){
     
      q31[i,j]=q3part1[i,j,1]-q3part2[i,j,1]-q3part3[i,j,1]+q3part4[i,j,1]
      q32[i,j]=q3part1[i,j,2]-q3part2[i,j,2]-q3part3[i,j,2]+q3part4[i,j,2]
      q33[i,j]=q3part1[i,j,3]-q3part2[i,j,3]-q3part3[i,j,3]+q3part4[i,j,3]
      q34[i,j]=q3part1[i,j,4]-q3part2[i,j,4]-q3part3[i,j,4]+q3part4[i,j,4]
      q35[i,j]=q3part1[i,j,5]-q3part2[i,j,5]-q3part3[i,j,5]+q3part4[i,j,5]
      q36[i,j]=q3part1[i,j,6]-q3part2[i,j,6]-q3part3[i,j,6]+q3part4[i,j,6]
      
    }
  }
  
  q3_1=apply(q31,1,sum)
  q3_2=apply(q32,1,sum)
  q3_3=apply(q33,1,sum)
  q3_4=apply(q34,1,sum)
  q3_5=apply(q35,1,sum)
  q3_6=apply(q36,1,sum)
  
  q3=rbind(q3_1,q3_2,q3_3,q3_4,q3_5,q3_6)
  
  
  
  #then estimate Sigmahat
  
  
  
  Sigma=matrix(0,n,mm^2)
  for (i in 1:n){
    for(k in 1:mm^2){
      if(k%%mm!=0){
        Sigma[i,k]=as.matrix((q1[,i]+q2[,i]+q3[,i])%*%t(q1[,i]+q2[,i]+q3[,i]))[ceiling(k/mm),k%%mm]
      }
      else{ Sigma[i,k]=as.matrix((q1[,i]+q2[,i]+q3[,i])%*%t(q1[,i]+q2[,i]+q3[,i]))[ceiling(k/mm),k%%mm+mm]
}
    }
  }
  
  Sigmahat=matrix(0,mm,mm)
  for(k in 1:mm^2){
    
    if(k%%mm!=0){
      Sigmahat[ceiling(k/mm),k%%mm]=sum(Sigma[,k])/n
    }
    else{    Sigmahat[ceiling(k/mm),k%%mm+mm]=sum(Sigma[,k])/n
    }
    
  }
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
  
  A2_1=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A2_1[i,j]=(X[i]-xbar1[j])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,1]*dN[i,j]
    }
  }
  
  A12_1=sum(A2_1)/n
  
  A2_2=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A2_2[i,j]=(X[i]-xbar1[j])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,2]*dN[i,j]
    }
  }
  
  A12_2=sum(A2_2)/n
  
  A2_3=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A2_3[i,j]=(X[i]-xbar1[j])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,3]*dN[i,j]
    }
  }
  
   
  A12_3=sum(A2_3)/n
  
  A2_4=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A2_4[i,j]=(X[i]-xbar1[j])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,4]*dN[i,j]
    }
  }
  
  
  A12_4=sum(A2_4)/n
  
  A2_5=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A2_5[i,j]=(X[i]-xbar1[j])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,5]*dN[i,j]
    }
  }
  
  
  A12_5=sum(A2_5)/n
  
  
  #estimate A21
  A3_1=matrix(0,n,nt)
  for (i in 1:n){
    A3_1[i,1]=1/2*(Z[i,1]-zbar1[1,1])*phihat[i,1]*X[i]*dallt[1]
    
    for (j in 2:nt){
      A3_1[i,j]=1/2*((Z[i,1]-zbar1[j,1])*phihat[i,j]*X[i]+(Z[i,1]-zbar1[j-1,1])*phihat[i,j-1]*X[i])*dallt[j]
    }
  }
  
  A21_1=sum(A3_1)/n
  
  A3_2=matrix(0,n,nt)
  for (i in 1:n){
    A3_2[i,1]=1/2*(Z[i,2]-zbar1[1,2])*phihat[i,1]*X[i]*dallt[1]
    
    for (j in 2:nt){
      A3_2[i,j]=1/2*((Z[i,2]-zbar1[j,2])*phihat[i,j]*X[i]+(Z[i,2]-zbar1[j-1,2])*phihat[i,j-1]*X[i])*dallt[j]
    }
  }
  
  A21_2=sum(A3_2)/n
  
  
  A3_3=matrix(0,n,nt)
  for (i in 1:n){
    A3_3[i,1]=1/2*(Z[i,3]-zbar1[1,3])*phihat[i,1]*X[i]*dallt[1]
    
    for (j in 2:nt){
      A3_3[i,j]=1/2*((Z[i,3]-zbar1[j,3])*phihat[i,j]*X[i]+(Z[i,3]-zbar1[j-1,3])*phihat[i,j-1]*X[i])*dallt[j]
    }
  }
  
  A21_3=sum(A3_3)/n
  
  
  A3_4=matrix(0,n,nt)
  for (i in 1:n){
    A3_4[i,1]=1/2*(Z[i,4]-zbar1[1,4])*phihat[i,1]*X[i]*dallt[1]
    
    for (j in 2:nt){
      A3_4[i,j]=1/2*((Z[i,4]-zbar1[j,4])*phihat[i,j]*X[i]+(Z[i,4]-zbar1[j-1,4])*phihat[i,j-1]*X[i])*dallt[j]
    }
  }
  
  A21_4=sum(A3_4)/n
  
  A3_5=matrix(0,n,nt)
  for (i in 1:n){
    A3_5[i,1]=1/2*(Z[i,5]-zbar1[1,5])*phihat[i,1]*X[i]*dallt[1]
    
    for (j in 2:nt){
      A3_5[i,j]=1/2*((Z[i,5]-zbar1[j,5])*phihat[i,j]*X[i]+(Z[i,5]-zbar1[j-1,5])*phihat[i,j-1]*X[i])*dallt[j]
    }
  }
  
  A21_5=sum(A3_5)/n
  
  #estimate A22
  
  A4_11=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_11[i,j]=(Z[i,1]-zbar1[j,1])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,1]*dN[i,j]
    }
  }
  
  A22_11=sum(A4_11)/n
  
  A4_12=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_12[i,j]=(Z[i,1]-zbar1[j,1])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,2]*dN[i,j]
    }
  }
  
  A22_12=sum(A4_12)/n
  
  A4_13=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_13[i,j]=(Z[i,1]-zbar1[j,1])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,3]*dN[i,j]
    }
  }
  
  A22_13=sum(A4_13)/n
  
  A4_14=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_14[i,j]=(Z[i,1]-zbar1[j,1])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,4]*dN[i,j]
    }
  }
  
  A22_14=sum(A4_14)/n
  
  A4_15=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_15[i,j]=(Z[i,1]-zbar1[j,1])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,5]*dN[i,j]
    }
  }
  
  A22_15=sum(A4_15)/n
  
  
  A4_21=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_21[i,j]=(Z[i,2]-zbar1[j,2])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,1]*dN[i,j]
    }
  }
  
  A22_21=sum(A4_21)/n
  
  
  A4_22=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_22[i,j]=(Z[i,2]-zbar1[j,2])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,2]*dN[i,j]
    }
  }
  
  A22_22=sum(A4_22)/n
  
  A4_23=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_23[i,j]=(Z[i,2]-zbar1[j,2])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,3]*dN[i,j]
    }
  }
  
  A22_23=sum(A4_23)/n
  
  A4_24=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_24[i,j]=(Z[i,2]-zbar1[j,2])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,4]*dN[i,j]
    }
  }
  
  A22_24=sum(A4_24)/n
  
  A4_25=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_25[i,j]=(Z[i,2]-zbar1[j,2])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,5]*dN[i,j]
    }
  }
  
  A22_25=sum(A4_25)/n
  
  
  A4_31=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_31[i,j]=(Z[i,3]-zbar1[j,3])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,1]*dN[i,j]
    }
  }
  
  A22_31=sum(A4_31)/n
  
  A4_32=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_32[i,j]=(Z[i,3]-zbar1[j,3])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,2]*dN[i,j]
    }
  }
  
  A22_32=sum(A4_32)/n
  
  
  
  
  A4_33=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_33[i,j]=(Z[i,3]-zbar1[j,3])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,3]*dN[i,j]
    }
  }
  
  A22_33=sum(A4_33)/n
  
  A4_34=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_34[i,j]=(Z[i,3]-zbar1[j,3])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,4]*dN[i,j]
    }
  }
  
  A22_34=sum(A4_34)/n
  
  A4_35=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_35[i,j]=(Z[i,3]-zbar1[j,3])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,5]*dN[i,j]
    }
  }
  
  A22_35=sum(A4_35)/n
  
  A4_41=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_41[i,j]=(Z[i,4]-zbar1[j,4])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,1]*dN[i,j]
    }
  }
  
  A22_41=sum(A4_41)/n
  
  
  A4_42=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_42[i,j]=(Z[i,4]-zbar1[j,4])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,2]*dN[i,j]
    }
  }
  
  A22_42=sum(A4_42)/n
  
  
  A4_43=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_43[i,j]=(Z[i,4]-zbar1[j,4])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,3]*dN[i,j]
    }
  }
  
  A22_43=sum(A4_43)/n
  
  
  A4_44=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_44[i,j]=(Z[i,4]-zbar1[j,4])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,4]*dN[i,j]
    }
  }
  
  A22_44=sum(A4_44)/n
  
  
  A4_45=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_45[i,j]=(Z[i,4]-zbar1[j,4])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,5]*dN[i,j]
    }
  }
  
  A22_45=sum(A4_45)/n
  
  
  A4_51=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_51[i,j]=(Z[i,5]-zbar1[j,5])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,1]*dN[i,j]
    }
  }
  
  A22_51=sum(A4_51)/n
  
  A4_52=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_52[i,j]=(Z[i,5]-zbar1[j,5])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,2]*dN[i,j]
    }
  }
  
  A22_52=sum(A4_52)/n
  
  A4_53=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_53[i,j]=(Z[i,5]-zbar1[j,5])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,3]*dN[i,j]
    }
  }
  
  A22_53=sum(A4_53)/n
  
  A4_54=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_54[i,j]=(Z[i,5]-zbar1[j,5])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,4]*dN[i,j]
    }
  }
  
  A22_54=sum(A4_54)/n
  
  A4_55=matrix(0,n,nt)
  for (i in 1:n){
    for (j in 1:nt){
      A4_55[i,j]=(Z[i,5]-zbar1[j,5])*phihat[i,j]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*Z[i,5]*dN[i,j]
    }
  }
  
  A22_55=sum(A4_55)/n
  
  
  
  #estimate Ahat
  Ahat=matrix(c(A11,A12_1,A12_2,A12_3,A12_4,A12_5,
                A21_1,A22_11,A22_12,A22_13,A22_14,A22_15,
                A21_2,A22_21,A22_22,A22_23,A22_24,A22_25,
                A21_3,A22_31,A22_32,A22_33,A22_34,A22_35,
                A21_4,A22_41,A22_42,A22_43,A22_44,A22_45,
                A21_5,A22_51,A22_52,A22_53,A22_54,A22_55),6,6,byrow=T)
  
  Ainv=solve(Ahat)
  
  #estimate Sigmastar
  
  Sigmastar=Ainv%*%Sigmahat%*%Ainv
  
  sgth=rep(0,6)
  sgth[1]=sqrt(Sigmastar[1,1])/sqrt(n)
  
  sgth[2]=sqrt(Sigmastar[2,2])/sqrt(n)
  
  sgth[3]=sqrt(Sigmastar[3,3])/sqrt(n)
  
  sgth[4]=sqrt(Sigmastar[4,4])/sqrt(n)
  sgth[5]=sqrt(Sigmastar[5,5])/sqrt(n)
  sgth[6]=sqrt(Sigmastar[6,6])/sqrt(n)
  
  iter=iter+1
}




u3part1=array(0,dim=c(n,nt,6))

for (i in 1:n){
  for (j in 1:nt){
    u3part1[i,j,1]=as.matrix((W[i,]-mwbar1[j,])*dNdd[i,j])[1,1]
    u3part1[i,j,2]=as.matrix((W[i,]-mwbar1[j,])*dNdd[i,j])[2,1]
    u3part1[i,j,3]=as.matrix((W[i,]-mwbar1[j,])*dNdd[i,j])[3,1]
    u3part1[i,j,4]=as.matrix((W[i,]-mwbar1[j,])*dNdd[i,j])[4,1]
    u3part1[i,j,5]=as.matrix((W[i,]-mwbar1[j,])*dNdd[i,j])[5,1]
    u3part1[i,j,6]=as.matrix((W[i,]-mwbar1[j,])*dNdd[i,j])[6,1]
    
  }
}

u3part2=array(0,dim=c(n,nt,6))
for (i in 1:n){
  u3part2[i,1,1]=as.matrix(1/2*(W[i,]-mwbar1[1,])*DT3[i,1]*as.numeric(hgm%*%W[i,])*dallt[1])[1,1]
  u3part2[i,1,2]=as.matrix(1/2*(W[i,]-mwbar1[1,])*DT3[i,1]*as.numeric(hgm%*%W[i,])*dallt[1])[2,1]
  u3part2[i,1,3]=as.matrix(1/2*(W[i,]-mwbar1[1,])*DT3[i,1]*as.numeric(hgm%*%W[i,])*dallt[1])[3,1]
  u3part2[i,1,4]=as.matrix(1/2*(W[i,]-mwbar1[1,])*DT3[i,1]*as.numeric(hgm%*%W[i,])*dallt[1])[4,1]
  u3part2[i,1,5]=as.matrix(1/2*(W[i,]-mwbar1[1,])*DT3[i,1]*as.numeric(hgm%*%W[i,])*dallt[1])[5,1]
  u3part2[i,1,6]=as.matrix(1/2*(W[i,]-mwbar1[1,])*DT3[i,1]*as.numeric(hgm%*%W[i,])*dallt[1])[6,1]
  
  for (j in 2:nt){
    u3part2[i,j,1]=as.matrix(1/2*((W[i,]-mwbar1[j,])*DT3[i,j]*as.numeric(hgm%*%W[i,])+(W[i,]-mwbar1[j-1,])*DT3[i,j-1]*as.numeric(hgm%*%W[i,]))*dallt[j])[1,1]
    u3part2[i,j,2]=as.matrix(1/2*((W[i,]-mwbar1[j,])*DT3[i,j]*as.numeric(hgm%*%W[i,])+(W[i,]-mwbar1[j-1,])*DT3[i,j-1]*as.numeric(hgm%*%W[i,]))*dallt[j])[2,1]
    u3part2[i,j,3]=as.matrix(1/2*((W[i,]-mwbar1[j,])*DT3[i,j]*as.numeric(hgm%*%W[i,])+(W[i,]-mwbar1[j-1,])*DT3[i,j-1]*as.numeric(hgm%*%W[i,]))*dallt[j])[3,1]
    u3part2[i,j,4]=as.matrix(1/2*((W[i,]-mwbar1[j,])*DT3[i,j]*as.numeric(hgm%*%W[i,])+(W[i,]-mwbar1[j-1,])*DT3[i,j-1]*as.numeric(hgm%*%W[i,]))*dallt[j])[4,1]
    u3part2[i,j,5]=as.matrix(1/2*((W[i,]-mwbar1[j,])*DT3[i,j]*as.numeric(hgm%*%W[i,])+(W[i,]-mwbar1[j-1,])*DT3[i,j-1]*as.numeric(hgm%*%W[i,]))*dallt[j])[5,1]
    u3part2[i,j,6]=as.matrix(1/2*((W[i,]-mwbar1[j,])*DT3[i,j]*as.numeric(hgm%*%W[i,])+(W[i,]-mwbar1[j-1,])*DT3[i,j-1]*as.numeric(hgm%*%W[i,]))*dallt[j])[6,1]
    
  }
}

u3part3=array(0,dim=c(n,nt,6))
for (i in 1:n){
  for (j in 1:nt){
    u3part3[i,j,1]=as.matrix((W[i,]-mwbar1[j,])*DT3[i,j]*dDelta0hatpart1[j])[1,1]
    u3part3[i,j,2]=as.matrix((W[i,]-mwbar1[j,])*DT3[i,j]*dDelta0hatpart1[j])[2,1]
    u3part3[i,j,3]=as.matrix((W[i,]-mwbar1[j,])*DT3[i,j]*dDelta0hatpart1[j])[3,1]
    u3part3[i,j,4]=as.matrix((W[i,]-mwbar1[j,])*DT3[i,j]*dDelta0hatpart1[j])[4,1]
    u3part3[i,j,5]=as.matrix((W[i,]-mwbar1[j,])*DT3[i,j]*dDelta0hatpart1[j])[5,1]
    u3part3[i,j,6]=as.matrix((W[i,]-mwbar1[j,])*DT3[i,j]*dDelta0hatpart1[j])[6,1]
    
  }
}

u3part4=array(0,dim=c(n,nt,6))
for (i in 1:n){
  for (j in 1:nt){
    u3part4[i,j,1]=as.matrix((W[i,]-mwbar1[j,])*DT3[i,j]*dDelta0hatpart2[j])[1,1]
    u3part4[i,j,2]=as.matrix((W[i,]-mwbar1[j,])*DT3[i,j]*dDelta0hatpart2[j])[2,1]
    u3part4[i,j,3]=as.matrix((W[i,]-mwbar1[j,])*DT3[i,j]*dDelta0hatpart2[j])[3,1]
    u3part4[i,j,4]=as.matrix((W[i,]-mwbar1[j,])*DT3[i,j]*dDelta0hatpart2[j])[4,1]
    u3part4[i,j,5]=as.matrix((W[i,]-mwbar1[j,])*DT3[i,j]*dDelta0hatpart2[j])[5,1]
    u3part4[i,j,6]=as.matrix((W[i,]-mwbar1[j,])*DT3[i,j]*dDelta0hatpart2[j])[6,1]
    
  }
}


u31=array(0,dim=c(n,nt))
u32=array(0,dim=c(n,nt))
u33=array(0,dim=c(n,nt))
u34=array(0,dim=c(n,nt))
u35=array(0,dim=c(n,nt))
u36=array(0,dim=c(n,nt))

for(i in 1:n){
  for(j in 1:nt){
    u31[i,j]=u3part1[i,j,1]-u3part2[i,j,1]-u3part3[i,j,1]+u3part4[i,j,1]
    u32[i,j]=u3part1[i,j,2]-u3part2[i,j,2]-u3part3[i,j,2]+u3part4[i,j,2]
    u33[i,j]=u3part1[i,j,3]-u3part2[i,j,3]-u3part3[i,j,3]+u3part4[i,j,3]
    u34[i,j]=u3part1[i,j,4]-u3part2[i,j,4]-u3part3[i,j,4]+u3part4[i,j,4]
    u35[i,j]=u3part1[i,j,5]-u3part2[i,j,5]-u3part3[i,j,5]+u3part4[i,j,5]
    u36[i,j]=u3part1[i,j,6]-u3part2[i,j,6]-u3part3[i,j,6]+u3part4[i,j,6]
    
  }
}

u3_1=apply(u31,1,sum)
u3_2=apply(u32,1,sum)
u3_3=apply(u33,1,sum)
u3_4=apply(u34,1,sum)
u3_5=apply(u35,1,sum)
u3_6=apply(u36,1,sum)

u3=rbind(u3_1,u3_2,u3_3,u3_4,u3_5,u3_6)

mv=matrix(c(30,50,10,0,0,1,0,0,1,0,1,1,0,1,1,0,1,2),6,byrow=T)

nv=ncol(mv)

colns=c(round(nt/5),2*round(nt/5),3*round(nt/5),4*round(nt/5))
ns=length(colns)


example_df <- data.frame(phihat,hbth,Z,X,dallt,dN,dNd,dallt1,nt,n,mv,nv,colns,u3,ns,
                         dNdd,partlmd1,partlmd2,R0,DT3,dDelta0hatpart1,dDelta0hatpart2,
                         hgm,W,mETA,allt,Omeinv,Ainv,q1,q2,q3)



save(phihat,file = "phihat.RData")
save(X,file = "X.RData")
save(Z,file = "Z.RData")
save(dallt,file = "dallt.RData")
save(dN,file = "dN.RData")
save(dNd,file = "dNd.RData")
save(dallt1,file = "dallt1.RData")
save(nt,file = "nt.RData")
save(n,file = "n.RData")
save(mv,file = "mv.RData")
save(nv,file = "nv.RData")
save(colns,file = "colns.RData")
save(u3,file = "u3.RData")
save(ns,file = "ns.RData")
save(dNdd,file = "dNdd.RData")
save(partlmd1,file = "partlmd1.RData")
save(partlmd2,file = "partlmd2.RData")
save(R0,file = "R0.RData")
save(DT3,file = "DT3.RData")
save(dDelta0hatpart1,file = "dDelta0hatpart1.RData")
save(dDelta0hatpart2,file = "dDelta0hatpart2.RData")
save(hgm,file = "hgm.RData")
save(W,file = "W.RData")
save(mETA,file = "mETA.RData")
save(allt,file = "allt.RData")
save(Omeinv,file = "Omeinv.RData")
save(Ainv,file = "Ainv.RData")
save(q1,file = "q1.RData")
save(q2,file = "q2.RData")
save(q3,file = "q3.RData")

save(hbth,file = "hbth.RData")
save(V,file = "V.RData")

save(etahat2_1,file = "etahat2_1_0814.RData")
save(component9,file = "component9.RData")

save(etahat2_2,file = "etahat2_2_0814.RData")
save(etahat2_3,file = "etahat2_3_0814.RData")

save(etahat11,file = "etahat11.RData")



























iter=1
while(iter<=1){
  DTv=matrix(0,n,nv)
  for(i in 1:n){
    for(j in 1:nv){
      DTv[i,j]=as.numeric(V[1,i]<=mv[1,j])*as.numeric(V[2,i]<=mv[2,j])*as.numeric(V[3,i]<=mv[3,j])*as.numeric(V[4,i]<=mv[4,j])*as.numeric(V[5,i]<=mv[5,j])*as.numeric(V[6,i]<=mv[6,j])
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
        etahat1part1[s,j,i]=Dhat[s,j,i]*phihat[i,s]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*dN[i,s]
        etahat1part3[s,j,i]=Dhat[s,j,i]*phihat[i,s]*partlmd1[s]
        etahat1part4[s,j,i]=Dhat[s,j,i]*phihat[i,s]*partlmd2[s]
        
      }
      etahat1part2[1,j,i]=1/2*Dhat[s,j,i]*phihat[i,1]*(hbth[1]*X[i])*dallt[1]
      for (s in 2:nt){
        etahat1part2[s,j,i]=1/2*(Dhat[s,j,i]*phihat[i,s]*(hbth[1]*X[i])+Dhat[s-1,j,i]*phihat[i,s-1]*(hbth[1]*X[i]))*dallt[s]
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
        
        etahat2part2[1,s,j,i]=1/2*component9[1,j,s]*DT3[i,1]*(hgm%*%W[i,])*dallt[1]
        for(u in 2:nt){
          etahat2part2[u,s,j,i]=1/2*(component9[u,j,s]*DT3[i,u]*(hgm%*%W[i,])+component9[u-1,j,s]*DT3[i,u-1]*(hgm%*%W[i,]))*dallt[u]
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
  
  
  
  Bsvpart1=array(0,c(nt,nv,n,6))
  Bsvpart2=array(0,c(nt,nv,n,6))
  Bsvpart3=array(0,c(nt,nv,n,6))
  Bsvpart4=array(0,c(nt,nv,n,6))
 
  
  for (i in 1:n){
    for(j in 1:nv){
      for (s in 1:nt){
        Bsvpart1[s,j,i,1]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*dN[i,s])[1,1]
        Bsvpart1[s,j,i,2]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*dN[i,s])[1,2]
        Bsvpart1[s,j,i,3]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*dN[i,s])[1,3]
        Bsvpart1[s,j,i,4]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*dN[i,s])[1,4]
        Bsvpart1[s,j,i,5]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*dN[i,s])[1,5]
        Bsvpart1[s,j,i,6]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*dN[i,s])[1,6]
        
        Bsvpart3[s,j,i,1]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*partlmd1[s])[1,1]
        Bsvpart3[s,j,i,2]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*partlmd1[s])[1,2]
        Bsvpart3[s,j,i,3]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*partlmd1[s])[1,3]
        Bsvpart3[s,j,i,4]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*partlmd1[s])[1,4]
        Bsvpart3[s,j,i,5]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*partlmd1[s])[1,5]
        Bsvpart3[s,j,i,6]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*partlmd1[s])[1,6]
        
        Bsvpart4[s,j,i,1]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*partlmd2[s])[1,1]
        Bsvpart4[s,j,i,2]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*partlmd2[s])[1,2]
        Bsvpart4[s,j,i,3]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*partlmd2[s])[1,3]
        Bsvpart4[s,j,i,4]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*partlmd2[s])[1,4]
        Bsvpart4[s,j,i,5]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*partlmd2[s])[1,5]
        Bsvpart4[s,j,i,6]=as.matrix(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*partlmd2[s])[1,6]
        
      }
      Bsvpart2[1,j,i,1]=as.matrix(1/2*Dhat[s,j,i]*(mETA[1,]+allt[1]*W[i,])%*%Omeinv*phihat[i,1]*(hbth[1]*X[i])*dallt[1])[1,1]
      Bsvpart2[1,j,i,2]=as.matrix(1/2*Dhat[s,j,i]*(mETA[1,]+allt[1]*W[i,])%*%Omeinv*phihat[i,1]*(hbth[1]*X[i])*dallt[1])[1,2]
      Bsvpart2[1,j,i,3]=as.matrix(1/2*Dhat[s,j,i]*(mETA[1,]+allt[1]*W[i,])%*%Omeinv*phihat[i,1]*(hbth[1]*X[i])*dallt[1])[1,3]
      Bsvpart2[1,j,i,4]=as.matrix(1/2*Dhat[s,j,i]*(mETA[1,]+allt[1]*W[i,])%*%Omeinv*phihat[i,1]*(hbth[1]*X[i])*dallt[1])[1,4]
      Bsvpart2[1,j,i,5]=as.matrix(1/2*Dhat[s,j,i]*(mETA[1,]+allt[1]*W[i,])%*%Omeinv*phihat[i,1]*(hbth[1]*X[i])*dallt[1])[1,5]
      Bsvpart2[1,j,i,6]=as.matrix(1/2*Dhat[s,j,i]*(mETA[1,]+allt[1]*W[i,])%*%Omeinv*phihat[i,1]*(hbth[1]*X[i])*dallt[1])[1,6]
      
      for(s in 2:nt){
        Bsvpart2[s,j,i,1]=as.matrix(1/2*(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*(hbth[1]*X[i])
                                         +Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s-1]*(hbth[1]*X[i]))*dallt[s])[1,1]
        Bsvpart2[s,j,i,2]=as.matrix(1/2*(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*(hbth[1]*X[i])
                                         +Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s-1]*(hbth[1]*X[i]))*dallt[s])[1,2]
        Bsvpart2[s,j,i,3]=as.matrix(1/2*(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*(hbth[1]*X[i])
                                         +Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s-1]*(hbth[1]*X[i]))*dallt[s])[1,3]
        Bsvpart2[s,j,i,4]=as.matrix(1/2*(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*(hbth[1]*X[i])
                                         +Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s-1]*(hbth[1]*X[i]))*dallt[s])[1,4]
        Bsvpart2[s,j,i,5]=as.matrix(1/2*(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*(hbth[1]*X[i])
                                         +Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s-1]*(hbth[1]*X[i]))*dallt[s])[1,5]
        Bsvpart2[s,j,i,6]=as.matrix(1/2*(Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s]*(hbth[1]*X[i])
                                         +Dhat[s,j,i]*(mETA[s,]+allt[s]*W[i,])%*%Omeinv*phihat[i,s-1]*(hbth[1]*X[i]))*dallt[s])[1,6]
        
      }
    }
  }
  
  
  Bsv=array(0,c(nt,nv,6))
  Bsv_0=array(0,c(nt,nv,6))
  
  sumBsvpart1_0=array(0,c(nt,nv,n,6))
  sumBsvpart1=array(0,c(nt,nv,6))
  sumBsvpart2_0=array(0,c(nt,nv,n,6))
  sumBsvpart2=array(0,c(nt,nv,6))
  sumBsvpart3_0=array(0,c(nt,nv,n,6))
  sumBsvpart3=array(0,c(nt,nv,6))
  sumBsvpart4_0=array(0,c(nt,nv,n,6))
  sumBsvpart4=array(0,c(nt,nv,6))
  
  
  
  
  for(k in 1:6){
    for(j in 1:nv){
      
      for (s in 1:nt){
        sumBsvpart1_0[s,j,1,k]=Bsvpart1[s,j,1,k]
        sumBsvpart2_0[s,j,1,k]=Bsvpart2[s,j,1,k]
        sumBsvpart3_0[s,j,1,k]=Bsvpart3[s,j,1,k]
        sumBsvpart4_0[s,j,1,k]=Bsvpart4[s,j,1,k]
        
        for (i in 2:n){
          sumBsvpart1_0[s,j,i,k]=sumBsvpart1_0[s,j,i-1,k]+Bsvpart1[s,j,i,k]
          sumBsvpart2_0[s,j,i,k]=sumBsvpart2_0[s,j,i-1,k]+Bsvpart2[s,j,i,k]
          sumBsvpart3_0[s,j,i,k]=sumBsvpart3_0[s,j,i-1,k]+Bsvpart3[s,j,i,k]
          sumBsvpart4_0[s,j,i,k]=sumBsvpart4_0[s,j,i-1,k]+Bsvpart4[s,j,i,k]
          
        }
      }
      
    }
  }
  
  
  for(k in 1:6){
    for(j in 1:nv){
      for (s in 1:nt){
        sumBsvpart1[s,j,k]=sumBsvpart1_0[s,j,n,k]
        sumBsvpart2[s,j,k]=sumBsvpart2_0[s,j,n,k]
        sumBsvpart3[s,j,k]=sumBsvpart3_0[s,j,n,k]
        sumBsvpart4[s,j,k]=sumBsvpart4_0[s,j,n,k]
        
        
      }
    }
  }
  
  for(k in 1:6){
    for (j in 1:nv){
      
      for (s in 1:nt){
        Bsv_0[s,j,k]=sumBsvpart1[s,j,k]-sumBsvpart2[s,j,k]-sumBsvpart3[s,j,k]+sumBsvpart4[s,j,k]
      }
    }
    
  }
  
  
  Bsv_0=Bsv_0/n
  
  
  for(k in 1:6){
    for(j in 1:nv){
      Bsv[1,j,k]=Bsv_0[1,j,k]
      for (s in 2:nt){
        Bsv[s,j,k]=Bsv[s-1,j,k]+Bsv_0[s,j,k]
      }
    }
  }
  
  Bsvv=array(0,c(ns,nv,6))  
  
  for(k in 1:6){
    for (j in 1:nv){
      for(s in 1:ns){
        Bsvv[s,j,k]=Bsv[colns[s],j,k]
      }
    }
  }
  
  
  
  #estimate etahat3
  
  etahat33=array(0,c(ns,nv,n))
  
  for( i in 1:n){
    for (j in 1:nv){
      for (s in 1:ns){
        etahat33[s,j,i]=(Bsvv[s,j,])%*%(u3[,i])
      }
    }
  }
  
  #estimate Ksv
  
  
  Ksvpart1=array(0,c(nt,nv,n))
  Ksvpart2=array(0,c(nt,nv,n))
  Ksvpart3=array(0,c(nt,nv,n))
  Ksvpart4=array(0,c(nt,nv,n))
  Ksvpart5=array(0,c(nt,nv,n))
  Ksvpart6=array(0,c(nt,nv,n))
  
  sumKsvpart1=array(0,c(nt,nv))
  sumKsvpart2=array(0,c(nt,nv))
  sumKsvpart3=array(0,c(nt,nv))
  sumKsvpart4=array(0,c(nt,nv))
  sumKsvpart5=array(0,c(nt,nv))
  sumKsvpart6=array(0,c(nt,nv))
  
  sumKsvpart1_0=array(0,c(nt,nv,n))
  sumKsvpart2_0=array(0,c(nt,nv,n))
  sumKsvpart3_0=array(0,c(nt,nv,n))
  sumKsvpart4_0=array(0,c(nt,nv,n))
  sumKsvpart5_0=array(0,c(nt,nv,n))
  sumKsvpart6_0=array(0,c(nt,nv,n))
  
  
  for (i in 1:n){
    for(j in 1:nv){
      Ksvpart1[1,j,i]=1/2*phihat[i,1]*Dhat[1,j,i]*X[i]*dallt[1]
      
      for (s in 2:nt){
        Ksvpart1[s,j,i]=1/2*(phihat[i,s]*Dhat[s,j,i]*X[i]+phihat[i,s-1]*Dhat[s-1,j,i]*X[i])*dallt[s]
        
      }
      for(s in 1:nt){
        Ksvpart2[s,j,i]=phihat[i,s]*Dhat[s,j,i]*Z[i,1]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*dN[i,s]
        Ksvpart3[s,j,i]=phihat[i,s]*Dhat[s,j,i]*Z[i,2]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*dN[i,s]
        Ksvpart4[s,j,i]=phihat[i,s]*Dhat[s,j,i]*Z[i,3]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,3]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*dN[i,s]
        Ksvpart5[s,j,i]=phihat[i,s]*Dhat[s,j,i]*Z[i,4]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,3]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*dN[i,s]
        Ksvpart6[s,j,i]=phihat[i,s]*Dhat[s,j,i]*Z[i,5]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,3]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*dN[i,s]
        
      }
    }
  }
  
  for(j in 1:nv){
    
    for (s in 1:nt){
      sumKsvpart1_0[s,j,1]=Ksvpart1[s,j,1]
      sumKsvpart2_0[s,j,1]=Ksvpart2[s,j,1]
      sumKsvpart3_0[s,j,1]=Ksvpart3[s,j,1]
      sumKsvpart4_0[s,j,1]=Ksvpart4[s,j,1]
      sumKsvpart5_0[s,j,1]=Ksvpart5[s,j,1]
      sumKsvpart6_0[s,j,1]=Ksvpart6[s,j,1]
      
      for (i in 2:n){
        sumKsvpart1_0[s,j,i]=sumKsvpart1_0[s,j,i-1]+Ksvpart1[s,j,i]
        sumKsvpart2_0[s,j,i]=sumKsvpart2_0[s,j,i-1]+Ksvpart2[s,j,i]
        sumKsvpart3_0[s,j,i]=sumKsvpart3_0[s,j,i-1]+Ksvpart3[s,j,i]
        sumKsvpart4_0[s,j,i]=sumKsvpart4_0[s,j,i-1]+Ksvpart4[s,j,i]
        sumKsvpart5_0[s,j,i]=sumKsvpart5_0[s,j,i-1]+Ksvpart5[s,j,i]
        sumKsvpart6_0[s,j,i]=sumKsvpart6_0[s,j,i-1]+Ksvpart6[s,j,i]
        
      }
    }
  }
  
  
  for(j in 1:nv){
    for (s in 1:nt){
      sumKsvpart1[s,j]=sumKsvpart1_0[s,j,n]/n
      sumKsvpart2[s,j]=sumKsvpart2_0[s,j,n]/n
      sumKsvpart3[s,j]=sumKsvpart3_0[s,j,n]/n
      sumKsvpart4[s,j]=sumKsvpart4_0[s,j,n]/n
      sumKsvpart5[s,j]=sumKsvpart5_0[s,j,n]/n
      sumKsvpart6[s,j]=sumKsvpart6_0[s,j,n]/n
      
    }
  }
  
  
  sumKsvpart11=array(0,c(nt,nv))
  sumKsvpart22=array(0,c(nt,nv))
  sumKsvpart33=array(0,c(nt,nv))
  sumKsvpart44=array(0,c(nt,nv))
  sumKsvpart55=array(0,c(nt,nv))
  sumKsvpart66=array(0,c(nt,nv))
  
  for(j in 1:nv){
    sumKsvpart11[1,j]=sumKsvpart1[1,j]
    sumKsvpart22[1,j]=sumKsvpart2[1,j]
    sumKsvpart33[1,j]=sumKsvpart3[1,j]
    sumKsvpart44[1,j]=sumKsvpart4[1,j]
    sumKsvpart55[1,j]=sumKsvpart5[1,j]
    sumKsvpart66[1,j]=sumKsvpart6[1,j]
    
    for(s in 2:nt){
      sumKsvpart11[s,j]=sumKsvpart11[s-1,j]+sumKsvpart1[s,j]
      sumKsvpart22[s,j]=sumKsvpart22[s-1,j]+sumKsvpart2[s,j]
      sumKsvpart33[s,j]=sumKsvpart33[s-1,j]+sumKsvpart3[s,j]
      sumKsvpart44[s,j]=sumKsvpart44[s-1,j]+sumKsvpart4[s,j]
      sumKsvpart55[s,j]=sumKsvpart55[s-1,j]+sumKsvpart5[s,j]
      sumKsvpart66[s,j]=sumKsvpart66[s-1,j]+sumKsvpart6[s,j]
      
    }}
  
  
  
  
  Ksvv=array(0,c(6,ns,nv))
  
  for(j in 1:nv){
    for (s in 1:ns){
      Ksvv[1,s,j]=-sumKsvpart11[colns[s],j]
      Ksvv[2,s,j]=-sumKsvpart22[colns[s],j]
      Ksvv[3,s,j]=-sumKsvpart33[colns[s],j]
      Ksvv[4,s,j]=-sumKsvpart44[colns[s],j]
      Ksvv[5,s,j]=-sumKsvpart55[colns[s],j]
      Ksvv[6,s,j]=-sumKsvpart66[colns[s],j]
      
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
        Fsvpart1[s,j,i]=DTv[i,j]*phihat[i,s]*exp(-(hbth[2]*Z[i,1]+hbth[3]*Z[i,2]+hbth[4]*Z[i,3]+hbth[5]*Z[i,4]+hbth[6]*Z[i,5]))*dN[i,s]
        Fsvpart3[s,j,i]=DTv[i,j]*phihat[i,s]*partlmd1[s]
        Fsvpart4[s,j,i]=DTv[i,j]*phihat[i,s]*partlmd2[s]
        
      }
      Fsvpart2[1,j,i]=1/2*DTv[i,j]*phihat[i,1]*(hbth[1]*X[i])*dallt[1]
      for (s in 2:nt){
        Fsvpart2[s,j,i]=1/2*(DTv[i,j]*phihat[i,s]*(hbth[1]*X[i])+DTv[i,j]*phihat[i,s-1]*(hbth[1]*X[i]))*dallt[s]
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
  
  
  
  as.numeric(mean(mf>= max(abs(obsF)))<=0.05)
  
  
  iter=iter+1
}



