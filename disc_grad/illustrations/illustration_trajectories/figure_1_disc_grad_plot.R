

initial.step <- function(fun,root,t0,y0){
  k4 <- fun(t0,y0)
  r4 <- root(t0,y0,k4)
  return(list(k4=k4,r4=r4,y_end=y0,t_end=t0))
}

bs32.step <- function(fun,root,eps,oldStep=NULL){
  
  y_start <- oldStep$y_end
  k1 <- oldStep$k4
  r1 <- oldStep$r4
  t_start <- oldStep$t_end
  t_end <- t_start+eps
  
  yt <- y_start + 0.5*eps*k1
  k2 <- fun(t_start+0.5*eps,yt)
  r2 <- root(t_start+0.5*eps,yt,k2)
  
  yt <- y_start + 0.75*eps*k2
  k3 <- fun(t_start+0.75*eps,yt)
  r3 <- root(t_start+0.75*eps,yt,k3)
  
  y_end <- y_start + eps*((2.0/9.0)*k1+(1.0/3.0)*k2 +(4.0/9.0)*k3)
  r_end <- eps*((2.0/9.0)*r1+(1.0/3.0)*r2+(4.0/9.0)*r3)
  
  k4 <- fun(t_end,y_end)
  r4 <- root(t_end,y_end,k4)
  
  y_low <- y_start + eps*((7.0/24.0)*k1 + (1.0/4.0)*k2 + (1.0/3.0)*k3 + (1.0/8.0)*k4)
  r_low <- eps*((7.0/24.0)*r1 + (1.0/4.0)*r2 + (1.0/3.0)*r3 + (1.0/8.0)*r4)
  
  
  
  return(list(y_start=y_start,y_end=y_end,r_end=r_end,k1=k1,r1=r1,k4=k4,r4=r4,
              t_start=t_start,t_end=t_end,eps=eps,y_low=y_low,r_low=r_low))
}


dense.y <- function(t,step){
  tau <- (t-step$t_start)
  if(any(tau>step$eps) || any(tau<0.0)) stop("bad time t in dense.y")
  t2 <- tau^2
  t3 <- tau^3
  e1 <- step$eps
  e2 <- step$eps^2
  e3 <- step$eps^3
  
  return((2*t3/e3-3*t2/e2+1)*step$y_start+
           (t3/e2-2*t2/e1+tau)*step$k1+
           (-2*t3/e3+3*t2/e2)*step$y_end+
           (t3/e2-t2/e1)*step$k4)
}



fun <- function(t,y){
  return(c(y[3:4],-(y[1:2]-y[5:6]),c(0,0)))
}
lsfun <- function(t,y,parms){
  return(list(c(y[3:4],-(y[1:2]-y[5:6]),c(0,0))))
}

root <- function(t,y,g){return(-1.0)} # not used

lsroot <- function(t,y,parms){return(y[1])}
eventfun <- function(t,y,parms){
  y.out <- y
  y.out[5] <- -y.out[5]
  print("event")
  print(t)
  return(y.out)
}


q0 <- c(-0.5,1.5)
v0 <- c(1.4,-2.5)

nstep <- 10
qs <- matrix(0.0,2,nstep+1)
qs[,1] <- q0
qsN <- matrix(0.0,2,nstep+1)
qsN[,1] <- q0
qsH <- matrix(0.0,2,nstep+1)
qsH[,1] <- q0

m0 <- c(-1,0)
y0 <- c(q0,v0,m0)
highRes <- deSolve::lsodar(y0,times=seq(from=0,to=3.0,length.out=1000),func=lsfun,rootfunc = lsroot,events = list(func=eventfun,root=TRUE),rtol=1e-10,atol=1.0e-10)
highResL <- deSolve::lsodar(y0,times=seq(from=0,to=3.0,length.out=1000),func=lsfun,rtol=1e-10,atol=1.0e-10)

oldStep <- initial.step(fun,root,0,y0)
h <- 0.3
for(i in 1:nstep){
  step <- bs32.step(fun,root,h,oldStep)
  qs[,i+1] <- step$y_end[1:2]
  oldStep <- step
}

oldStep <- initial.step(fun,root,0,y0)

for(i in 1:nstep){
  step <- bs32.step(fun,root,h,oldStep)
  qsN[,i+1] <- step$y_end[1:2]
  if(step$y_end[1]*step$y_start[1]<0.0){
    yy <- step$y_end
    yy[5] <- -yy[5]
    oldStep <- initial.step(fun,root,step$t_end,yy)
  } else {
  oldStep <- step
  }
}

oldStep <- initial.step(fun,root,0,y0)

for(i in 1:nstep){
  step <- bs32.step(fun,root,h,oldStep)
  qsH[,i+1] <- step$y_end[1:2]
  if(step$y_end[1]*step$y_start[1]< -1e-6){
    ii <- i
    fn <- function(t){dense.y(t,step)[1]}
    tt <- uniroot(fn,interval=c(step$t_start,step$t_end))$root
    print(tt)
    yy <- dense.y(tt,step)
    qsH <- cbind(qsH,yy[1:2])
    yy[5] <- -yy[5]
    oldStep <- initial.step(fun,root,tt,yy)
  } else {
    oldStep <- step
  }
}


pdf("disc_grad/illustrations/illustration_trajectories/stepIll.pdf",width=14,height = 7)
plot(highResL[,2],highResL[,3],xlim=c(-0.5,1.0),type="l",lty=2,lwd=3,col="darkgrey",
     xlab=latex2exp::TeX("$q_1$"),ylab=latex2exp::TeX("$q_2$"))

points(qsN[1,],qsN[2,],pch=19,col="green",cex=2)
points(qsH[1,],qsH[2,],pch=18,col="blue",cex=1.5)
lines(highRes[,2],highRes[,3])
arrows(highRes[10,2],highRes[10,3],highRes[11,2],highRes[11,3])
lines(c(0,0),c(-10,10),col="red")

m1 <- (qsH[1,ii+1]+qsH[1,12])/2
m2 <- (qsH[2,ii+1]+qsH[2,12])/2

r <- sqrt((qsH[1,ii+1]+qsH[1,12])^2 + (qsH[2,ii+1]+qsH[2,12])^2)

t1<-uniroot(function(theta){m1+r*sin(theta)-qsH[1,ii+1]},interval=c(1,3))$root
t2<-uniroot(function(theta){m1+r*sin(theta)-qsH[1,12]},interval=c(-1,1))$root
th <- seq(from=t1,to=t2,length.out=1000)
lines(m1+r*sin(th),m2+r*cos(th),col="blue",lwd=1)

thm <- th[990]
arrows(m1+r*sin(thm),m2+r*cos(thm),m1+r*sin(thm-0.01),m2+r*cos(thm-0.01),col="blue")


dev.off()


