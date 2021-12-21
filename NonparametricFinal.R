n<-1000
x<-c(rep(0,n))
components.random <- sample(1:2,prob=c(0.5,0.5),size=n,replace=TRUE)
mus <- c(-1.5,1.5)
sds <- c(sqrt(3),sqrt(0.1))
x <- rnorm(n=n,mean=mus[components.random],sd=sds[components.random])
f<-function(x){0.5*dnorm(x,-1.5,sqrt(3))+0.5*dnorm(x,1.5,sqrt(0.1))}
xgrid<-seq(-10,10,0.02)
ygrid<-f(xgrid)
GK<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)}
GK2<-function(u){((1/sqrt(2*pi))*exp(-(u^2)/2))^2}
RK<-integrate(GK2,-Inf,Inf)$value
GK3<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)*u^2}
mu2K<-integrate(GK3,-Inf,Inf)$value
df2<-function(x){
  part1<-((x+1.5)^2-3)/(54*sqrt(6*pi))*exp(-(x+1.5)^2/6)
  part2<-(50*(x-1.5)^2-5)/(sqrt(0.2*pi))*exp(-(x-1.5)^2/0.2)
  y<-part1+part2
  return(y)
}
f2<-function(x){(df2(x))^2}
Rf2<-integrate(f2,-Inf,Inf)$value
h_AMISE<-(RK/(n*Rf2*mu2K^2))^(1/5)
h_AMISE
par(mfrow=c(1,2))
layout_mat <- matrix(seq(1,2), nrow = 1, ncol = 2,
                     byrow = TRUE)
layout(mat = layout_mat, 
       widths = c(1,1), respect =TRUE)
plot(xgrid,ygrid,type = "l",lty=2,xlab = "x",ylab ="f",main = "(a)")
lines(density(x,h_AMISE))
plot(xgrid,ygrid,type = "l",lty=2,xlab = "x",ylab ="f",main = "(b)")
lines(density(x,0.32))
### Fig 2
n<-100
xs<-rnorm(n)
hs<-seq(0.3,0.5,0.02)
MISE<-c(rep(0,length(hs)))
AMISE<-c(rep(0,length(hs)))
l2<-c(rep(0,length(hs)))
## RK
GK<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)}
GK2<-function(u){((1/sqrt(2*pi))*exp(-(u^2)/2))^2}
RK<-integrate(GK2,-Inf,Inf)$value
## Rf2
df2<-function(x){(dnorm(x)*(x^2-1))^2}
Rf2<-integrate(df2,-Inf,Inf)$value
##
for (i in 1:length(hs)) {
  h<-hs[i]
  # MISE[i]<-(1/(2*sqrt(pi)))*(1/(n*h)+(1-1/n)/(1+h^2)+1-(2^(3/2))/(sqrt(2+h^2)))
  MISE[i]<-(1/(n*h)+(1-1/n)/(sqrt(1+h^2))+1-(2^(3/2))/(sqrt(2+h^2)))/(2*sqrt(pi))
  ###
  AMISE[i]<-(h^4*Rf2)/4+(RK)/(n*h)
  l2[i]<-AMISE[i]-1/(3*pi*n)
}
plot(hs,MISE,type = "l",lty=3,xlab = "h",ylab ="error",ylim = range(0,MISE,AMISE,l2))
lines(hs,AMISE,lty=2)
lines(hs,l2)
### Fig3
ns<-seq(20,200,10)
hMISE<-c(rep(0,length(ns)))
for (i in 1:length(ns)) {
    hs<-seq(0.01,1,0.001)
    hMISE[i]<-hs[which.min((1/(ns[i]*hs)+(1-1/ns[i])/(sqrt(1+hs^2))+1-(2^(3/2))/(sqrt(2+hs^2)))/(2*sqrt(pi)))]
}
MISE<-c(rep(0,length(ns)))
AMISE<-c(rep(0,length(ns)))
l2<-c(rep(0,length(ns)))
for (i in 1:length(ns)) {
  n<-ns[i]
  h<-hMISE[i]
  MISE[i]<-(1/(n*h)+(1-1/n)/(sqrt(1+h^2))+1-(2^(3/2))/(sqrt(2+h^2)))/(2*sqrt(pi))
  AMISE[i]<-(h^4*Rf2)/4+(RK)/(n*h)
  l2[i]<-AMISE[i]-1/(3*pi*n)
}
plot(ns,MISE,type = "l",lty=3,xlab = "n",ylab ="error",ylim = range(0,MISE,AMISE,l2))
lines(ns,AMISE,lty=2)
lines(ns,l2)
### Fig 4
n<-100
x<-rnorm(n)
f<-function(x){dnorm(x)}
xgrid<-seq(-10,10,0.02)
ygrid<-f(xgrid)
###
par(mfrow=c(1,2))
layout_mat <- matrix(seq(1,2), nrow = 1, ncol = 2,
                     byrow = TRUE)
layout(mat = layout_mat, 
       widths = c(1,1), respect =TRUE)
###
plot(xgrid,ygrid,type = "l",lty=2,xlab = "x",ylab ="f",main = "(a)")
lines(density(x,0.26))
plot(xgrid,ygrid,type = "l",lty=2,xlab = "x",ylab ="f",main = "(b)")
lines(density(x,0.61))
### Fig 5
f1<-function(x){dnorm(x)}#a
f2<-function(x){0.6*dnorm(x,-1.5,1)+0.4*dnorm(x,1.5,1)}#b
f3<-function(x){0.5*dnorm(x,-1.5,sqrt(2))+0.5*dnorm(x,1.5,sqrt(0.2))}#c
f4<-function(x){0.5*dnorm(x,0,sqrt(3))+0.5*dnorm(x,0,sqrt(0.1))}#d
f5<-function(x){(1/3)*(dnorm(x,-1.5,1)+dnorm(x)+dnorm(x,1.5,1))}#e
f6<-function(x){(1/3)*(dnorm(x,-2.5,1)+dnorm(x)+dnorm(x,2.5,1))}#f
f7<-function(x){0.3*dnorm(x,-2,sqrt(2))+0.4*dnorm(x,0,sqrt(0.1))+0.3*dnorm(x,2,sqrt(2))}#g
f8<-function(x){(1/3)*(dnorm(x,-2,1)+dnorm(x,0,sqrt(0.3))+dnorm(x,1,sqrt(0.1)))}#h
f9<-function(x){0.1*dnorm(x,-2,sqrt(0.05))+0.8*dnorm(x,0,2)+0.1*dnorm(x,2,sqrt(0.05))}#i
xgrid<-seq(-10,10,0.02)
y1grid<-f1(xgrid)
y2grid<-f2(xgrid)
y3grid<-f3(xgrid)
y4grid<-f4(xgrid)
y5grid<-f5(xgrid)
y6grid<-f6(xgrid)
y7grid<-f7(xgrid)
y8grid<-f8(xgrid)
y9grid<-f9(xgrid)
###
par(mfrow=c(3,3))
layout_mat <- matrix(seq(1,9), nrow = 3, ncol = 3,
                     byrow = TRUE)
layout(mat = layout_mat, 
       heights = c(1,1,1),
       widths = c(1,1,1), respect =TRUE)
par(mar = c(3.75, 3.75, 2, 0.6))
###
plot(xgrid,y1grid,type = "l",lty=1,xlab = "x",ylab ="f",main = "(a)")
plot(xgrid,y2grid,type = "l",lty=1,xlab = "x",ylab ="f",main = "(b)")
plot(xgrid,y3grid,type = "l",lty=1,xlab = "x",ylab ="f",main = "(c)")
plot(xgrid,y4grid,type = "l",lty=1,xlab = "x",ylab ="f",main = "(d)")
plot(xgrid,y5grid,type = "l",lty=1,xlab = "x",ylab ="f",main = "(e)")
plot(xgrid,y6grid,type = "l",lty=1,xlab = "x",ylab ="f",main = "(f)")
plot(xgrid,y7grid,type = "l",lty=1,xlab = "x",ylab ="f",main = "(g)")
plot(xgrid,y8grid,type = "l",lty=1,xlab = "x",ylab ="f",main = "(h)")
plot(xgrid,y9grid,type = "l",lty=1,xlab = "x",ylab ="f",main = "(i)")
### Table of different h
n<-200
hs<-seq(0.0000001,1,0.0000001)
####### (a)
hMISE<-hs[which.min((1/(n*hs)+(1-1/n)/(sqrt(1+hs^2))+1-(2^(3/2))/(sqrt(2+hs^2)))/(2*sqrt(pi)))]
hMISE
###
GK<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)}
GK2<-function(u){((1/sqrt(2*pi))*exp(-(u^2)/2))^2}
RK<-integrate(GK2,-Inf,Inf)$value
GK3<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)*u^2}
mu2K<-integrate(GK3,-Inf,Inf)$value
df2<-function(x){(x^2)*dnorm(x)-dnorm(x)}
f2<-function(x){(df2(x))^2}
Rf2<-integrate(f2,-Inf,Inf)$value
h_AMISE<-(RK/(n*Rf2*mu2K^2))^(1/5)
h_AMISE
###
a<-(1/4)*(mu2K^2)*Rf2
b<-RK/n
c<-(1/4)*(h_AMISE^4)*(mu2K^2)*Rf2+(1/n)*((RK/h_AMISE)+(1/(3*pi)))
equation<-function(h){a*h^4+b/h-c}
hgrid<-seq(0.01,1,0.01)
ygrid<-equation(hgrid)
plot(hgrid,ygrid,type="l",lty=1)
abline(h=0)
h_star<-uniroot(equation,lower = 0.4,upper = 0.6)$root
h_star
####### (b)
### search for h_MISE
MISE<-function(h){
  n<-200
  w<-matrix(c(0.6,0.4),2,1)
  Omega2<-matrix(c(dnorm(0,0,sqrt(2*h^2+2)),dnorm(3,0,sqrt(2*h^2+2)),dnorm(-3,0,sqrt(2*h^2+2)),dnorm(0,0,sqrt(2*h^2+2))),2,2)
  Omega1<-matrix(c(dnorm(0,0,sqrt(h^2+2)),dnorm(3,0,sqrt(h^2+2)),dnorm(-3,0,sqrt(h^2+2)),dnorm(0,0,sqrt(h^2+2))),2,2)
  Omega0<-matrix(c(dnorm(0,0,sqrt(2)),dnorm(3,0,sqrt(2)),dnorm(-3,0,sqrt(2)),dnorm(0,0,sqrt(2))),2,2)
  y<-1/(2*n*h*sqrt(pi))+t(w)%*%((1-1/n)*Omega2-2*Omega1+Omega0)%*%w
  return(y)
}
h_MISE<-optimize(MISE,c(0,1))$minimum
h_MISE
### search for h_AMISE
n<-200
GK<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)}
GK2<-function(u){((1/sqrt(2*pi))*exp(-(u^2)/2))^2}
RK<-integrate(GK2,-Inf,Inf)$value
GK3<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)*u^2}
mu2K<-integrate(GK3,-Inf,Inf)$value
df2<-function(x){
  part1<-(0.6/sqrt(2*pi))*((x+1.5)^2-1)*exp(-(x+1.5)^2/2)
  part2<-(0.4/sqrt(2*pi))*((x-1.5)^2-1)*exp(-(x-1.5)^2/2)
  y<-part1+part2
  return(y)
}
f2<-function(x){(df2(x))^2}
Rf2<-integrate(f2,-Inf,Inf)$value
h_AMISE<-(RK/(n*Rf2*mu2K^2))^(1/5)
h_AMISE
### search for h_star
w<-matrix(c(0.6,0.4),2,1)
S<-matrix(c(1,1),2,1)
M<-matrix(c(-1.5,1.5),2,1)
V<-t(w)%*%S+t(w)%*%(M^2)-(t(w)%*%M)^2
n<-200
a<-(1/4)*(mu2K^2)*Rf2
b<-RK/n
c<-(1/4)*(h_AMISE^4)*(mu2K^2)*Rf2+(1/n)*((RK/h_AMISE)+(1/(3*pi*sqrt(c(V)))))###
equation<-function(h){a*h^4+b/h-c}
hgrid<-seq(0.01,1,0.01)
ygrid<-equation(hgrid)
plot(hgrid,ygrid,type="l",lty=1)
abline(h=0)
h_star<-uniroot(equation,lower = 0.4,upper = 0.6)$root
h_star
####### (c)
### search for h_MISE
MISE<-function(h){
  n<-200
  w<-matrix(c(0.5,0.5),2,1)
  Omega2<-matrix(c(dnorm(0,0,sqrt(2*h^2+2+2)),dnorm(3,0,sqrt(2*h^2+2.2)),dnorm(-3,0,sqrt(2*h^2+2.2)),dnorm(0,0,sqrt(2*h^2+0.4))),2,2)
  Omega1<-matrix(c(dnorm(0,0,sqrt(h^2+2+2)),dnorm(3,0,sqrt(h^2+2.2)),dnorm(-3,0,sqrt(h^2+2.2)),dnorm(0,0,sqrt(h^2+0.4))),2,2)
  Omega0<-matrix(c(dnorm(0,0,sqrt(2+2)),dnorm(3,0,sqrt(2.2)),dnorm(-3,0,sqrt(2.2)),dnorm(0,0,sqrt(0.4))),2,2)
  y<-1/(2*n*h*sqrt(pi))+t(w)%*%((1-1/n)*Omega2-2*Omega1+Omega0)%*%w
  return(y)
}
h_MISE<-optimize(MISE,c(0,1))$minimum
h_MISE
### search for h_AMISE
n<-200
GK<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)}
GK2<-function(u){((1/sqrt(2*pi))*exp(-(u^2)/2))^2}
RK<-integrate(GK2,-Inf,Inf)$value
GK3<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)*u^2}
mu2K<-integrate(GK3,-Inf,Inf)$value
df2<-function(x){
  part1<-(0.5/sqrt(2*pi*2))*(((x+1.5)^2)/4-1/2)*exp(-(x+1.5)^2/(2*2))
  part2<-(0.5/sqrt(2*pi*0.2))*(((x-1.5)^2)/(0.2^2)-1/0.2)*exp(-(x-1.5)^2/(2*0.2))
  y<-part1+part2
  return(y)
}
f2<-function(x){(df2(x))^2}
Rf2<-integrate(f2,-Inf,Inf)$value
h_AMISE<-(RK/(n*Rf2*mu2K^2))^(1/5)
h_AMISE
### search for h_star
n<-200
w<-matrix(c(0.5,0.5),2,1)
S<-matrix(c(2,0.2),2,1)
M<-matrix(c(-1.5,1.5),2,1)
V<-t(w)%*%S+t(w)%*%(M^2)-(t(w)%*%M)^2
a<-(1/4)*(mu2K^2)*Rf2
b<-RK/n
c<-(1/4)*(h_AMISE^4)*(mu2K^2)*Rf2+(1/n)*((RK/h_AMISE)+(1/(3*pi*sqrt(c(V)))))###
equation<-function(h){a*h^4+b/h-c}
hgrid<-seq(0.01,1,0.01)
ygrid<-equation(hgrid)
plot(hgrid,ygrid,type="l",lty=1)
abline(h=0)
h_star<-uniroot(equation,lower = 0.2,upper = 0.4)$root
h_star
####### (d)
### search for h_MISE
MISE<-function(h){
  n<-200
  w<-matrix(c(0.5,0.5),2,1)
  Omega2<-matrix(c(dnorm(0,0,sqrt(2*h^2+6)),dnorm(0,0,sqrt(2*h^2+3.1)),dnorm(0,0,sqrt(2*h^2+3.1)),dnorm(0,0,sqrt(2*h^2+0.2))),2,2)
  Omega1<-matrix(c(dnorm(0,0,sqrt(h^2+6)),dnorm(0,0,sqrt(h^2+3.1)),dnorm(0,0,sqrt(h^2+3.1)),dnorm(0,0,sqrt(h^2+0.2))),2,2)
  Omega0<-matrix(c(dnorm(0,0,sqrt(6)),dnorm(0,0,sqrt(3.1)),dnorm(0,0,sqrt(3.1)),dnorm(0,0,sqrt(0.2))),2,2)
  y<-1/(2*n*h*sqrt(pi))+t(w)%*%((1-1/n)*Omega2-2*Omega1+Omega0)%*%w
  return(y)
}
h_MISE<-optimize(MISE,c(0,1))$minimum
h_MISE
### search for h_AMISE
n<-200
GK<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)}
GK2<-function(u){((1/sqrt(2*pi))*exp(-(u^2)/2))^2}
RK<-integrate(GK2,-Inf,Inf)$value
GK3<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)*u^2}
mu2K<-integrate(GK3,-Inf,Inf)$value
df2<-function(x){
  part1<-(0.5/sqrt(2*pi*3))*(((x)^2)/9-1/3)*exp(-(x)^2/(2*3))
  part2<-(0.5/sqrt(2*pi*0.1))*(((x)^2)/(0.1^2)-1/0.1)*exp(-(x)^2/(2*0.1))
  y<-part1+part2
  return(y)
}
f2<-function(x){(df2(x))^2}
Rf2<-integrate(f2,-Inf,Inf)$value
h_AMISE<-(RK/(n*Rf2*mu2K^2))^(1/5)
h_AMISE
### search for h_star
n<-200
w<-matrix(c(0.5,0.5),2,1)
S<-matrix(c(3,0.1),2,1)
M<-matrix(c(0,0),2,1)
V<-t(w)%*%S+t(w)%*%(M^2)-(t(w)%*%M)^2
a<-(1/4)*(mu2K^2)*Rf2
b<-RK/n
c<-(1/4)*(h_AMISE^4)*(mu2K^2)*Rf2+(1/n)*((RK/h_AMISE)+(1/(3*pi*sqrt(c(V)))))###
equation<-function(h){a*h^4+b/h-c}
hgrid<-seq(0.01,1,0.01)
ygrid<-equation(hgrid)
plot(hgrid,ygrid,type="l",lty=1)
abline(h=0)
h_star<-uniroot(equation,lower = 0.15,upper = 0.2)$root
h_star
####### (e)
### search for h_MISE
MISE<-function(h){
  n<-200
  w<-matrix(c(1/3,1/3,1/3),3,1)
  Omega2<-matrix(c(dnorm(0,0,sqrt(2*h^2+2)),dnorm(1.5,0,sqrt(2*h^2+2)),dnorm(3,0,sqrt(2*h^2+2)),
                   dnorm(-1.5,0,sqrt(2*h^2+2)),dnorm(0,0,sqrt(2*h^2+2)),dnorm(1.5,0,sqrt(2*h^2+2)),
                   dnorm(-3,0,sqrt(2*h^2+2)),dnorm(-1.5,0,sqrt(2*h^2+2)),dnorm(0,0,sqrt(2*h^2+2))),3,3)
  Omega1<-matrix(c(dnorm(0,0,sqrt(h^2+2)),dnorm(1.5,0,sqrt(h^2+2)),dnorm(3,0,sqrt(h^2+2)),
                   dnorm(-1.5,0,sqrt(h^2+2)),dnorm(0,0,sqrt(h^2+2)),dnorm(1.5,0,sqrt(h^2+2)),
                   dnorm(-3,0,sqrt(h^2+2)),dnorm(-1.5,0,sqrt(h^2+2)),dnorm(0,0,sqrt(h^2+2))),3,3)
  Omega0<-matrix(c(dnorm(0,0,sqrt(3.1)),dnorm(1.5,0,sqrt(3.1)),dnorm(3,0,sqrt(3.1)),
                   dnorm(-1.5,0,sqrt(3.1)),dnorm(0,0,sqrt(3.1)),dnorm(1.5,0,sqrt(3.1)),
                   dnorm(-3,0,sqrt(3.1)),dnorm(-1.5,0,sqrt(3.1)),dnorm(0,0,sqrt(3.1))),3,3)
  y<-1/(2*n*h*sqrt(pi))+t(w)%*%((1-1/n)*Omega2-2*Omega1+Omega0)%*%w
  return(y)
}
h_MISE<-optimize(MISE,c(0,1))$minimum
h_MISE
### search for h_AMISE
n<-200
GK<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)}
GK2<-function(u){((1/sqrt(2*pi))*exp(-(u^2)/2))^2}
RK<-integrate(GK2,-Inf,Inf)$value
GK3<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)*u^2}
mu2K<-integrate(GK3,-Inf,Inf)$value
df2<-function(x){
  part1<-((1/3)/sqrt(2*pi))*(((x+1.5)^2)-1)*exp(-(x+1.5)^2/(2*1))
  part2<-((1/3)/sqrt(2*pi))*(((x)^2)-1)*exp(-(x)^2/(2*1))
  part3<-((1/3)/sqrt(2*pi))*(((x-1.5)^2)-1)*exp(-(x-1.5)^2/(2*1))
  y<-part1+part2+part3
  return(y)
}
f2<-function(x){(df2(x))^2}
Rf2<-integrate(f2,-Inf,Inf)$value
h_AMISE<-(RK/(n*Rf2*mu2K^2))^(1/5)
h_AMISE
### search for h_star
n<-200
w<-matrix(c(1/3,1/3,1/3),3,1)
S<-matrix(c(1,1,1),3,1)
M<-matrix(c(-1.5,0,1.5),3,1)
V<-t(w)%*%S+t(w)%*%(M^2)-(t(w)%*%M)^2
a<-(1/4)*(mu2K^2)*Rf2
b<-RK/n
c<-(1/4)*(h_AMISE^4)*(mu2K^2)*Rf2+(1/n)*((RK/h_AMISE)+(1/(3*pi*sqrt(c(V)))))###
equation<-function(h){a*h^4+b/h-c}
hgrid<-seq(0.01,1,0.01)
ygrid<-equation(hgrid)
plot(hgrid,ygrid,type="l",lty=1)
abline(h=0)
h_star<-uniroot(equation,lower = 0.6,upper = 0.9)$root
h_star
####### (f)
### search for h_MISE
MISE<-function(h){
  n<-200
  w<-matrix(c(1/3,1/3,1/3),3,1)
  Omega2<-matrix(c(dnorm(0,0,sqrt(2*h^2+2)),dnorm(2.5,0,sqrt(2*h^2+2)),dnorm(5,0,sqrt(2*h^2+2)),
                   dnorm(-2.5,0,sqrt(2*h^2+2)),dnorm(0,0,sqrt(2*h^2+2)),dnorm(2.5,0,sqrt(2*h^2+2)),
                   dnorm(-5,0,sqrt(2*h^2+2)),dnorm(-2.5,0,sqrt(2*h^2+2)),dnorm(0,0,sqrt(2*h^2+2))),3,3)
  Omega1<-matrix(c(dnorm(0,0,sqrt(h^2+2)),dnorm(2.5,0,sqrt(h^2+2)),dnorm(5,0,sqrt(h^2+2)),
                   dnorm(-2.5,0,sqrt(h^2+2)),dnorm(0,0,sqrt(h^2+2)),dnorm(2.5,0,sqrt(h^2+2)),
                   dnorm(-5,0,sqrt(h^2+2)),dnorm(-2.5,0,sqrt(h^2+2)),dnorm(0,0,sqrt(h^2+2))),3,3)
  Omega0<-matrix(c(dnorm(0,0,sqrt(3.1)),dnorm(2.5,0,sqrt(3.1)),dnorm(5,0,sqrt(3.1)),
                   dnorm(-2.5,0,sqrt(3.1)),dnorm(0,0,sqrt(3.1)),dnorm(2.5,0,sqrt(3.1)),
                   dnorm(-5,0,sqrt(3.1)),dnorm(-2.5,0,sqrt(3.1)),dnorm(0,0,sqrt(3.1))),3,3)
  y<-1/(2*n*h*sqrt(pi))+t(w)%*%((1-1/n)*Omega2-2*Omega1+Omega0)%*%w
  return(y)
}
h_MISE<-optimize(MISE,c(0,1))$minimum
h_MISE
### search for h_AMISE
n<-200
GK<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)}
GK2<-function(u){((1/sqrt(2*pi))*exp(-(u^2)/2))^2}
RK<-integrate(GK2,-Inf,Inf)$value
GK3<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)*u^2}
mu2K<-integrate(GK3,-Inf,Inf)$value
df2<-function(x){
  part1<-((1/3)/sqrt(2*pi))*(((x+2.5)^2)-1)*exp(-(x+2.5)^2/(2*1))
  part2<-((1/3)/sqrt(2*pi))*(((x)^2)-1)*exp(-(x)^2/(2*1))
  part3<-((1/3)/sqrt(2*pi))*(((x-2.5)^2)-1)*exp(-(x-2.5)^2/(2*1))
  y<-part1+part2+part3
  return(y)
}
f2<-function(x){(df2(x))^2}
Rf2<-integrate(f2,-Inf,Inf)$value
h_AMISE<-(RK/(n*Rf2*mu2K^2))^(1/5)
h_AMISE
### search for h_star
n<-200
w<-matrix(c(1/3,1/3,1/3),3,1)
S<-matrix(c(1,1,1),3,1)
M<-matrix(c(-2.5,0,2.5),3,1)
V<-t(w)%*%S+t(w)%*%(M^2)-(t(w)%*%M)^2
a<-(1/4)*(mu2K^2)*Rf2
b<-RK/n
c<-(1/4)*(h_AMISE^4)*(mu2K^2)*Rf2+(1/n)*((RK/h_AMISE)+(1/(3*pi*sqrt(c(V)))))###
equation<-function(h){a*h^4+b/h-c}
hgrid<-seq(0.01,1,0.01)
ygrid<-equation(hgrid)
plot(hgrid,ygrid,type="l",lty=1)
abline(h=0)
h_star<-uniroot(equation,lower = 0.6,upper = 0.9)$root
h_star
####### (g)
### search for h_MISE
MISE<-function(h){
  n<-200
  w<-matrix(c(0.3,0.4,0.3),3,1)
  Omega2<-matrix(c(dnorm(0,0,sqrt(2*h^2+4)),dnorm(2,0,sqrt(2*h^2+2.1)),dnorm(4,0,sqrt(2*h^2+4)),
                   dnorm(-2,0,sqrt(2*h^2+2.1)),dnorm(0,0,sqrt(2*h^2+0.2)),dnorm(2,0,sqrt(2*h^2+2.1)),
                   dnorm(-4,0,sqrt(2*h^2+4)),dnorm(-2,0,sqrt(2*h^2+2.1)),dnorm(0,0,sqrt(2*h^2+4))),3,3)
  Omega1<-matrix(c(dnorm(0,0,sqrt(h^2+4)),dnorm(2,0,sqrt(h^2+2.1)),dnorm(4,0,sqrt(h^2+4)),
                   dnorm(-2,0,sqrt(h^2+2.1)),dnorm(0,0,sqrt(h^2+0.2)),dnorm(2,0,sqrt(h^2+2.1)),
                   dnorm(-4,0,sqrt(h^2+4)),dnorm(-2,0,sqrt(h^2+2.1)),dnorm(0,0,sqrt(h^2+4))),3,3)
  Omega0<-matrix(c(dnorm(0,0,sqrt(4)),dnorm(2,0,sqrt(2.1)),dnorm(4,0,sqrt(4)),
                   dnorm(-2,0,sqrt(2.1)),dnorm(0,0,sqrt(0.2)),dnorm(2,0,sqrt(2.1)),
                   dnorm(-4,0,sqrt(4)),dnorm(-2,0,sqrt(2.1)),dnorm(0,0,sqrt(4))),3,3)
  y<-1/(2*n*h*sqrt(pi))+t(w)%*%((1-1/n)*Omega2-2*Omega1+Omega0)%*%w
  return(y)
}
h_MISE<-optimize(MISE,c(0,1))$minimum
h_MISE
### search for h_AMISE
n<-200
GK<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)}
GK2<-function(u){((1/sqrt(2*pi))*exp(-(u^2)/2))^2}
RK<-integrate(GK2,-Inf,Inf)$value
GK3<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)*u^2}
mu2K<-integrate(GK3,-Inf,Inf)$value
df2<-function(x){
  part1<-((0.3)/sqrt(2*pi*2))*(((x+2)^2)/4-1/2)*exp(-(x+2.5)^2/(2*2))
  part2<-((0.4)/sqrt(2*pi*0.1))*(((x)^2)/0.01-1/0.1)*exp(-(x)^2/(2*0.1))
  part3<-((0.3)/sqrt(2*pi*2))*(((x-2)^2)/4-1/2)*exp(-(x-2.5)^2/(2*2))
  y<-part1+part2+part3
  return(y)
}
f2<-function(x){(df2(x))^2}
Rf2<-integrate(f2,-Inf,Inf)$value
h_AMISE<-(RK/(n*Rf2*mu2K^2))^(1/5)
h_AMISE
### search for h_star
n<-200
w<-matrix(c(0.3,0.4,0.3),3,1)
S<-matrix(c(2,0.1,2),3,1)
M<-matrix(c(-2,0,2),3,1)
V<-t(w)%*%S+t(w)%*%(M^2)-(t(w)%*%M)^2
a<-(1/4)*(mu2K^2)*Rf2
b<-RK/n
c<-(1/4)*(h_AMISE^4)*(mu2K^2)*Rf2+(1/n)*((RK/h_AMISE)+(1/(3*pi*sqrt(c(V)))))###
equation<-function(h){a*h^4+b/h-c}
hgrid<-seq(0.01,1,0.01)
ygrid<-equation(hgrid)
plot(hgrid,ygrid,type="l",lty=1)
abline(h=0)
h_star<-uniroot(equation,lower = 0.15,upper = 0.25)$root
h_star
####### (h)
### search for h_MISE
MISE<-function(h){
  n<-200
  w<-matrix(c(1/3,1/3,1/3),3,1)
  Omega2<-matrix(c(dnorm(0,0,sqrt(2*h^2+2)),dnorm(2,0,sqrt(2*h^2+1.3)),dnorm(3,0,sqrt(2*h^2+1.1)),
                   dnorm(-2,0,sqrt(2*h^2+1.3)),dnorm(0,0,sqrt(2*h^2+0.6)),dnorm(1,0,sqrt(2*h^2+0.4)),
                   dnorm(-3,0,sqrt(2*h^2+1.1)),dnorm(-1,0,sqrt(2*h^2+0.4)),dnorm(0,0,sqrt(2*h^2+0.2))),3,3)
  Omega1<-matrix(c(dnorm(0,0,sqrt(h^2+2)),dnorm(2,0,sqrt(h^2+1.3)),dnorm(3,0,sqrt(h^2+1.1)),
                   dnorm(-2,0,sqrt(h^2+1.3)),dnorm(0,0,sqrt(h^2+0.6)),dnorm(1,0,sqrt(h^2+0.4)),
                   dnorm(-3,0,sqrt(h^2+1.1)),dnorm(-1,0,sqrt(h^2+0.4)),dnorm(0,0,sqrt(h^2+0.2))),3,3)
  Omega0<-matrix(c(dnorm(0,0,sqrt(2)),dnorm(2,0,sqrt(1.3)),dnorm(3,0,sqrt(1.1)),
                   dnorm(-2,0,sqrt(1.3)),dnorm(0,0,sqrt(0.6)),dnorm(1,0,sqrt(0.4)),
                   dnorm(-3,0,sqrt(1.1)),dnorm(-1,0,sqrt(0.4)),dnorm(0,0,sqrt(0.2))),3,3)
  y<-1/(2*n*h*sqrt(pi))+t(w)%*%((1-1/n)*Omega2-2*Omega1+Omega0)%*%w
  return(y)
}
h_MISE<-optimize(MISE,c(0,1))$minimum
h_MISE
### search for h_AMISE
n<-200
GK<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)}
GK2<-function(u){((1/sqrt(2*pi))*exp(-(u^2)/2))^2}
RK<-integrate(GK2,-Inf,Inf)$value
GK3<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)*u^2}
mu2K<-integrate(GK3,-Inf,Inf)$value
df2<-function(x){
  part1<-((1/3)/sqrt(2*pi))*(((x+2)^2)-1)*exp(-(x+2)^2/(2*1))
  part2<-((1/3)/sqrt(2*pi*0.3))*(((x)^2)/0.3^2-1/0.3)*exp(-(x)^2/(2*0.3))
  part3<-((1/3)/sqrt(2*pi*0.1))*(((x-1)^2)/0.1^2-1/0.1)*exp(-(x-1)^2/(2*0.1))
  y<-part1+part2+part3
  return(y)
}
f2<-function(x){(df2(x))^2}
Rf2<-integrate(f2,-Inf,Inf)$value
h_AMISE<-(RK/(n*Rf2*mu2K^2))^(1/5)
h_AMISE
### search for h_star
n<-200
w<-matrix(c(1/3,1/3,1/3),3,1)
S<-matrix(c(1,0.3,0.1),3,1)
M<-matrix(c(-2,0,1),3,1)
V<-t(w)%*%S+t(w)%*%(M^2)-(t(w)%*%M)^2
a<-(1/4)*(mu2K^2)*Rf2
b<-RK/n
c<-(1/4)*(h_AMISE^4)*(mu2K^2)*Rf2+(1/n)*((RK/h_AMISE)+(1/(3*pi*sqrt(c(V)))))###
equation<-function(h){a*h^4+b/h-c}
hgrid<-seq(0.01,1,0.01)
ygrid<-equation(hgrid)
plot(hgrid,ygrid,type="l",lty=1)
abline(h=0)
h_star<-uniroot(equation,lower = 0.2,upper = 0.3)$root
h_star
####### (i)
### search for h_MISE
MISE<-function(h){
  n<-200
  w<-matrix(c(0.1,0.8,0.1),3,1)
  Omega2<-matrix(c(dnorm(0,0,sqrt(2*h^2+0.1)),dnorm(2,0,sqrt(2*h^2+4.05)),dnorm(4,0,sqrt(2*h^2+0.1)),
                   dnorm(-2,0,sqrt(2*h^2+4.05)),dnorm(0,0,sqrt(2*h^2+8)),dnorm(2,0,sqrt(2*h^2+4.05)),
                   dnorm(-4,0,sqrt(2*h^2+0.1)),dnorm(-2,0,sqrt(2*h^2+4.05)),dnorm(0,0,sqrt(2*h^2+0.1))),3,3)
  Omega1<-matrix(c(dnorm(0,0,sqrt(h^2+0.1)),dnorm(2,0,sqrt(h^2+4.05)),dnorm(4,0,sqrt(h^2+0.1)),
                   dnorm(-2,0,sqrt(h^2+4.05)),dnorm(0,0,sqrt(h^2+8)),dnorm(2,0,sqrt(h^2+4.05)),
                   dnorm(-4,0,sqrt(h^2+0.1)),dnorm(-2,0,sqrt(h^2+4.05)),dnorm(0,0,sqrt(h^2+0.1))),3,3)
  Omega0<-matrix(c(dnorm(0,0,sqrt(0.1)),dnorm(2,0,sqrt(4.05)),dnorm(4,0,sqrt(0.1)),
                   dnorm(-2,0,sqrt(4.05)),dnorm(0,0,sqrt(8)),dnorm(2,0,sqrt(4.05)),
                   dnorm(-4,0,sqrt(0.1)),dnorm(-2,0,sqrt(4.05)),dnorm(0,0,sqrt(0.1))),3,3)
  y<-1/(2*n*h*sqrt(pi))+t(w)%*%((1-1/n)*Omega2-2*Omega1+Omega0)%*%w
  return(y)
}
h_MISE<-optimize(MISE,c(0,1))$minimum
h_MISE
### search for h_AMISE
n<-200
GK<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)}
GK2<-function(u){((1/sqrt(2*pi))*exp(-(u^2)/2))^2}
RK<-integrate(GK2,-Inf,Inf)$value
GK3<-function(u){(1/sqrt(2*pi))*exp(-(u^2)/2)*u^2}
mu2K<-integrate(GK3,-Inf,Inf)$value
df2<-function(x){
  part1<-((0.1)/sqrt(2*pi*0.05))*(((x+2)^2)/0.05^2-1/0.05)*exp(-(x+2)^2/(2*0.05))
  part2<-((0.8)/sqrt(2*pi*4))*(((x)^2)/4^2-1/4)*exp(-(x)^2/(2*4))
  part3<-((0.1)/sqrt(2*pi*0.05))*(((x-2)^2)/0.05^2-1/0.05)*exp(-(x-2)^2/(2*0.05))
  y<-part1+part2+part3
  return(y)
}
f2<-function(x){(df2(x))^2}
Rf2<-integrate(f2,-Inf,Inf)$value
h_AMISE<-(RK/(n*Rf2*mu2K^2))^(1/5)
h_AMISE
### search for h_star
n<-200
w<-matrix(c(0.1,0.8,0.1),3,1)
S<-matrix(c(0.05,4,0.05),3,1)
M<-matrix(c(-2,0,2),3,1)
V<-t(w)%*%S+t(w)%*%(M^2)-(t(w)%*%M)^2
a<-(1/4)*(mu2K^2)*Rf2
b<-RK/n
c<-(1/4)*(h_AMISE^4)*(mu2K^2)*Rf2+(1/n)*((RK/h_AMISE)+(1/(3*pi*sqrt(c(V)))))###
equation<-function(h){a*h^4+b/h-c}
hgrid<-seq(0.01,1,0.01)
ygrid<-equation(hgrid)
plot(hgrid,ygrid,type="l",lty=1)
abline(h=0)
h_star<-uniroot(equation,lower = 0.2,upper = 0.3)$root
h_star
###
result<-matrix(c(0.3670978, 0.3830395,0.4524205,
                 0.4293822, 0.4745363,0.5106065,
                 0.2167415, 0.2372072,0.2455497,
                 0.1531062, 0.1672649,0.1738521,
                 0.5886252, 0.6348775,0.7263407,
                 0.5299138, 0.6619018,0.628397,
                 0.1675756, 0.187514,0.1867968,
                 0.1853014, 0.218298,0.2111211,
                 0.1795002, 0.2453338,0.2002967),9,3,byrow = T)
result