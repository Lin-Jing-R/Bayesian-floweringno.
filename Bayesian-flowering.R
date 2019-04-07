#####normal linear regression
lm.f18<-glm(flower18~treat*group, family=binomial(log))
d<-c(0:200)
anova.18<-anova(lm.f18,test="Chisq")
#recover seperate intercapts and slopes using contrastsmatrix
matrix<-cbind(1,constraste(group))%*%matrix(coef(lm.f),2)
#prediction for mean response with s.e.
pre.eu<-predict.glm(lm.f18,data.frame(treatment=d,group=factor(rep("EU",length(d)),
                                                               levels=levels(group))),
                    type="response",se.fit=T)
pre.am<-predict.glm(lm.f18,data.frame(treatment=d,group=factor(rep("NAint",length(d)),
                                                               levels=levels(group))),
                    type="response",se.fit=T)
par(mfrow=c(3,1))
plot(diam,flower18,ylim=c(0,90))
plot(d,pred.eu$fit,type="l",ylim=(0,80))
lines(d,(pred.eu$fit+pred.eu$se),lty=2)
lines(d,(pred.eu$fit-pred.eu$se),lty=2)
lines(d,pred.am$fit)
summary(lm.f18)

flo<-c(0,40)
val<-seq(5,20,by=5)
plot(d,exp(matrix[1,1]+matrix[1,2]*d),type='l',
     ylim=c(0,80))
for(i in 1:length(val)){
  pv<-exp(matrix[1,1]+matrix[1,2]*val[i])
  dist<-dpois(flo,val)
  lines((20*dist+val[i]),flo,type='s')
  abline(v=val[i])
}

flo<-flower18; flo[flo>0]<-1
lm.f<-glm(flo~treat*group, family=binomial(logit))
anova.f<-anova(lm.f,test="Chisq")

###transfer
linkit<-function(lvec)exp(lvec)/(1+exp(lvec))
newmatrix<-cbind(1,constraste(group))%*%matrix(coef(lm.f),2)
dseq<-c(0:200)
#prediction
pre.eu<-predict.glm(lm.f,data.frame(treatment=dseq,group=factor(rep("EU",length(dseq)),
                                                                levels=levels(group))),
                    type="response",se.fit=T)
pre.am<-predict.glm(lm.f,data.frame(treatment=dseq,group=factor(rep("NAint",length(dseq)),
                                                                levels=levels(group))),
                    type="response",se.fit=T)
plot(dseq,pred.eu$fit,type="l",ylim=(0,80))
lines(d,(pred.eu$fit+pred.eu$se),lty=2)
lines(d,(pred.eu$fit-pred.eu$se),lty=2)
lines(d,pred.am$fit)

####Bayeaian linear regression
nutri<-treatment
c1<-flo
c2<-flower18
g<-group
b0<-20;b1<-2;a0<-.05
param<-c(b0,b1,a0)
likfit<-function(param){
  b0<-param[1]
  b1<-param[2]
  a0<-param[3]
  ##shoot with no flower
  zd<-nutri[c1==0]
  pc<-a0*zd^2
  theta1<-pnorm(zd,b0,b1)
  pz<-log((1-theta1)+theta1*exp(-pc))
  ##with flowers
  kd<-nutri[c1>0]
  kc<-c2[c1>0]
  qc<-a0*kd^2
  theta2<-pnorm(kd,b0,b1)
  pk<-log(theta2)-qc+kc*log(qc)-gamma(kc+1)
  return(-sum(pz,pk))
}
out<-optim(param,likfit,lower=c(10,1,0.01),upper=c(30,10,1),
           ,etod="L-BFGS-B")
btot<-out$par[1]
b1tot<-out$par[2]
atot<-out$par[3]
ltot<-out$value
ntot<-length(flo)

d1<-diam[group=="EU"]
c1<-flo[group=="EU"]
c2<-flower18[group=="EU"]
out<-optim(param,likfit,lower-c(10,1,0.01),upper=c(30,10,1),
           method="L-BFGS-B")
b.eu<-out$par[1]
b1.eu<-out$par[2]
a.eu<-out$par[3]
l.eu<-out$value
n.eu<-length(c1)

d1<-diam[group=="NAint"]
c1<-flo[group=="NAint"]
c2<-flower18[group=="NAint"]
out<-optim(param,likfit,lower-c(10,1,0.001),upper=c(30,10,1),
           method="L-BFGS-B")
b.na<-out$par[1]
b1.na<-out$par[2]
a.na<-out$par[3]
l.na<-out$value
n.na<-length(c1)

nboot<-2000
mle<-matrix(NA,nrow=nboot,ncol=3)
dseq<-seq(0,30,length=100)
pred.mat<-matrix(NA,nrow=nboot,ncol=100)
pred.con<-matrix(NA,nrow=nboot,ncol=100)
flo.t<-treatment[group=="EU"]
c1.t<-flo[group=="EU"]
c2.t<-flower18[group=="EU"]
nt.t<-length(flo.t)
for(i in 1:nboot){
  cindex<-sample(nt.t,replace=T)
  d1<-flo.t[cindex]
  c1<-flo[cindex]
  c2<-flower18[cindex]
  out<-optim(param,likfit,lower=c(10,1,0.001),upper=c(30,10,1),
             method="L-BFGS-B")
  mle[i,]<-out$par
  print(c(i,mle[i,]))
  pred.mat[i,]<-pnorm(dseq,mle[i,1],mle[i,2])
  pred.flo[i,]<-mle[i,3]*dseq2
}

separ<-sqrt(diag((var(mle))))
cib0<-quantile(mle[,1],c(0.025,0.975))##95% confident level
cib1<-quantile(mle[,2],c(0.025,0.975))
cia0<-quantile(mle[,3],c(0.025,0.975))
bootout<-
  rbind(c(b.na,bi.na,a.na),separ,cbind(cib0,cib1,cia0))

###density plot
plot(density(mle[,1],width=1.5),type='l',xlim=c(10,30),ylim=c(0,.5))
abine(v=cib0,lty=2)
plot(density(mle[,2],width=1.5),type='l',xlim=c(10,20),ylim=c(0,.5))
abine(v=cib1,lty=2)
plot(density(mle[,3],width=0.05),type='l',xlim=c(0,0.05),ylim=c(0,120))
abine(v=cia0,lty=2)

#predict responses
ci.m<-matrix(NA,nrow=100,ncol=2)
ci.c<-matrix(NA,nrow=100,ncol=2)
for(i in 1:100){
  ci.m[i,]<-quantile(pred.mat[,i],c(0.025,0.975))
  ci.c[i,]<-quantile(pred.mat[,i],c(0.025,0.975))
}
plot(treatment,jitter(flo),xlim=c(0,30))
lines(dseq,predmat)
lines(dseq,ci.m[,1],lty=2)
lines(dseq,ci.m[,2],lyt=2)
plot(rtreatment,predcon)
lines(dseq,ci.c[,1],lty=2)
lines(dseq,ci.c[,2],lty=2)

