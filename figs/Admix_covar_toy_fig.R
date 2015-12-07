library(igraph)
iArrows <- igraph:::igraph.Arrows   #http://kbroman.wordpress.com/2012/10/11/curved-arrows-in-r/

alpha=1
par(mar=c(0,3.5,0.1,1))
loc<-c(0,1,2,3)

covar.d<-function(d){
exp(-alpha*sqrt(d^2)^1.3)
}
cols<-c("red","black","blue","orange")

d<-seq(0,3,length=1000)
plot(d,covar.d(d-loc[1]),type="l",col=cols[1],lwd=2,ylim=c(-0.3,1),axes=FALSE,xlab="",ylab="",cex=1.5)
axis(side=2,at=c(0,0.5,1),cex=1.2)
mtext(side=2,line=2,"Covariance",cex=1.4)
abline(h=0)

lines(d,covar.d(d-loc[3]),col=cols[3],lwd=2)
lines(d,covar.d(d-loc[4]),col=cols[4],lwd=2)
w1<-0.4
text(x=1, y = (1-w1)*covar.d(loc[1]-loc[2]) + w1*covar.d(loc[1]-loc[4]),col=cols[1],"B-A",cex=1.5)  #points(x=1, y = (1-w1)*covar.d(0-1) + w1*covar.d(0-3),col=cols[1],pch=19)
text(x=1, y = (1-w1)*covar.d(loc[3]-loc[2]) + w1*covar.d(loc[3]-loc[4]),col=cols[3],"B-C",cex=1.5)
text(x=1, y = (1-w1)*covar.d(loc[4]-loc[1]) + w1*covar.d(0),col=cols[4],"B-D" ,cex=1.5)
text(x=1, y = (1-w1)*(1-w1)*1 + w1*w1*1 + 2*(1-w1)*w1*covar.d(loc[4]-loc[1]),col=cols[2],"B-B" ,cex=1.5)

text(c(0,1,2,3),rep(-0.1,4),c("A","B","C",'D'),col=cols,cex=1.5)
iArrows(loc[4], -0.18, loc[1], -0.18,
          h.lwd=2, sh.lwd=2, sh.col="black",
          curve=0.2 , width=1, size=0.7)
          dev.copy2pdf(file="~/Desktop/Admix_covar_toy_fig.pdf")
          
#          points(x=1, y = (1-w1)*covar.d(0-1) + w1*covar.d(0-3),col=cols[1],pch=19)
#points(x=1, y = (1-w1)*covar.d(2-1) + w1*covar.d(3-2),col=cols[3],pch=19 )
#points(x=1, y = (1-w1)*covar.d(3-1) + w1*covar.d(3-3),col=cols[4],pch=19 )
#points(x=1, y = (1-w1)*(1-w1)*1 + w1*w1*1 + 2*(1-w1)*w1*covar.d(3-1),col=cols[2],pch=19 )
