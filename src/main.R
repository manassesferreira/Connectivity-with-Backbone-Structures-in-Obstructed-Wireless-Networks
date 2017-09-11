

library(menina)
require(plotrix) #draw.circle

g=15
mu=15
r=0.175
eps=0.05

D=definaDispositivo(g,mu)

C=computeComponente(D,g,r,eps)
  #print(C)
  print(length(C$x))

  par(mar=c(0,0,0,0)+0.1)
  plot( -1,-1, xlim=c(0,g-1), ylim=c(0,g-1), xaxt='n', yaxt='n', ann=FALSE)
  #points( C$x,C$y, cex=0.4, col=C$comp, pch=16)#text( C$x,C$y, labels=C$comp, cex= 0.15)
  abline(v=-1,lty=3);	abline(h=-1,lty=3); for(i in -1:g+1){ abline(v=i,lty=3); abline(h=i,lty=3);}
  #for(dx in 0:g-1){	for(dy in 0:g-1){
	#	segments(dx+0+eps,dy+0+eps,dx+0+eps,dy+1-eps); segments(dx+0+eps,dy+1-eps,dx+1-eps,dy+1-eps);
	#	segments(dx+1-eps,dy+1-eps,dx+1-eps,dy+0+eps); segments(dx+1-eps,dy+0+eps,dx+0+eps,dy+0+eps);
	#}}

  ##para observar conectividade para um dado dispositivo
  #node=round(runif(1,1,length(C$x)))
  #points(C$x[node],C$y[node],col=2)
  #draw.circle(C$x[node],C$y[node],r,border=2)
  #C$u=C$u+1;  C$v=C$v+1; #No R, os indíces começam em 1, ao invés de 0, como no Cpp.
  #lista=sort(unique(c(C$u[which(C$v[]==node)],C$v[which(C$u[]==node)])))
  #points(C$x[lista],C$y[lista],col=2,cex=0.3)

  print(length(C$u))
  print(length(unique(C$comp)))

B=busqueBase(C,D,g,r,eps)
  print(length(B$x))
  #points(B$x,B$y,pch=2,col=2)

A=assenteAcesso(B,C,D,g,r,eps)
  print(length(A$x))
  points(A$x,A$y,pch=2,col=2,cex=0.8) #pch=17


