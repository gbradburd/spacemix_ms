################################################################
################################################################
#	make figures for SpaceMix paper
################################################################
################################################################

procrusteez <- function(obs.locs,target.locs,k,source.locs = NULL,option){
	require(vegan)
	proc.loc <- procrustes(obs.locs,target.locs,scale=TRUE)
	if(option==1){
		proc.pop.loc <- proc.loc$scale * target.locs %*% proc.loc$rotation + matrix(proc.loc$translation,nrow=k,ncol=2,byrow=TRUE)
	} else if(option==2){
		proc.pop.loc <- proc.loc$scale * source.locs %*% proc.loc$rotation + matrix(proc.loc$translation,nrow=k,ncol=2,byrow=TRUE)
	}
	return(proc.pop.loc)	
}

Covariance <- function(a0,aD,a2,GeoDist) {
	covariance <- (1/a0)*exp(-(aD*GeoDist)^a2)
	return(covariance)
}

admixed.Covariance <- function(covariance,admix.proportions,nugget,k,inv.mean.sample.sizes,ident.mat){
	# recover()
	if(any(admix.proportions !=0)){
		w_k <- admix.proportions/2
		admixed.Covariance <- 	tcrossprod((1-w_k),(1-w_k)) * 	covariance[1:k,1:k] + 
								tcrossprod((1-w_k),(w_k)) 	* 	covariance[1:k,(k+1):(2*k)] +
								tcrossprod(w_k,(1-w_k)) 	*	covariance[(k+1):(2*k),1:k] +
								tcrossprod(w_k,w_k)			*	covariance[(k+1):(2*k),(k+1):(2*k)]
		admixed.Covariance <- admixed.Covariance + ident.mat * nugget + ident.mat * inv.mean.sample.sizes
	} else {
		admixed.Covariance <- covariance[1:k,1:k] + ident.mat * nugget + ident.mat * inv.mean.sample.sizes
	}
	return(admixed.Covariance)
}

transformed.Covariance <- function(covariance,projection.matrix){
	transformed.covariance <- 
		crossprod(projection.matrix,covariance) %*% projection.matrix
	return(transformed.covariance)		
}

fade.admixture.source.points <- function(pop.cols,admix.proportions){
	faded.colors <- numeric(length(pop.cols))
	for(i in 1:length(pop.cols)){
		faded.colors[i] <- adjustcolor(pop.cols[i],admix.proportions[i])
	}
	return(faded.colors)
}


################################
#	WARBLER POP FIGS
################################
load("~/Desktop/Dropbox/space.mix/data/warblers/warbler_spacemix/pop/warbler_pop_spaceruns/real_prior1/warbler_pop_dataset.Robj")
	pops <- row.names(warbler.pop.coords)
	k <- length(pops)
	pop.col <- numeric(k)
	pop.col[match(c("XN"),pops)] <- "orange"
	pop.col[match(c("EM","TB","LN"),pops)] <- "gold"
	pop.col[match(c("MN","ML","PA","KL","KS","PKA","PKB"),pops)] <- "mediumseagreen"
	pop.col[match(c("AA"),pops)] <- "dodgerblue2"
	pop.col[match(c("YK","TL","AB"),pops)] <- "dodgerblue2"		#"slateblue4"
	pop.col[match(c("ST","UY","IL","AN","BK","TA","SL"),pops)] <- "red"
	pop.col[match(c("TU"),pops)] <- "slateblue4"

################
#	RealPrior1
################
load("~/Desktop/Dropbox/space.mix/data/warblers/warbler_spacemix/pop/warbler_pop_spaceruns/real_prior1/warb_pop_spaceruns_realpr1_LongRun/warb_pop_spaceruns_realpr1space_MCMC_output1.Robj")
	best <- which.max(Prob)
	target.coords <- procrusteez(warbler.pop.coords,population.coordinates[[best]][1:k,],k,option=1)
	source.coords <- procrusteez(warbler.pop.coords,population.coordinates[[best]][1:k,],k,source.locs=population.coordinates[[best]][(k+1):(2*k),],option=2)
	pop.plot.cols <- fade.admixture.source.points(pop.col,admix.proportions[,best])

	png(file="~/Desktop/Dropbox/space.mix/ms/population_warbler_map_realpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(68,116), #realpr1: c(68,115), realpr2: c(53,101), randpr: c(71,114)
					ylim=c(26,53), #realpr1: c(26,53), realpr2: c(25,54), randpr: c(71,114)
					xlab="long",
					ylab="lat")
				text(target.coords[c(1:k),],
						labels=pops,
						col=pop.col,
						font=2,cex=0.9)
				points(source.coords,
							col=pop.plot.cols,
							pch=20)
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col= pop.plot.cols,
					lwd=admix.proportions[,best],#last.params$admix.proportions,
					length=0.1)
				box(lwd=2)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/population_warbler_map_no_arrows_realpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(68,101), #realpr1: c(68,101), realpr2: c(68,101), randpr: c(71,99)
					ylim=c(25,55), #realpr1: c(25,55), realpr2: c(25,55), randpr: c(20.5,56)
					xlab="long",
					ylab="lat")
				text(target.coords[c(1:k),],
						labels=pops,
						col=pop.col,
						font=2,cex=0.9)
				box(lwd=2)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/population_warbler_admix_values_nugget_realpr1.png",res=300,width=4.5*300,height=5*300,pointsize=9)
		#quartz(width=5,height=5,pointsize=9)
		order <- c(1,3:15,2,16:22)
		par(mfrow=c(2,1),mar=c(1,4,0,1),oma=c(1,1,1,1))
			plot(admix.proportions[order,best]/2,type='n',
					xlab="",
					xaxt='n',
					ylab="estimated adxmiture proportion")
				text(admix.proportions[order,best]/2,
						labels=pops[order],
						col=pop.col[order],
						font=2,cex=0.9)
				box(lwd=2)
			plot(nugget[order,best],type='n',
					xaxt='n',
					ylab="estimated population nugget")
				text(nugget[order,best],
						labels=pops[order],
						col=pop.col[order],
						font=2,cex=0.9)
				box(lwd=2)
	dev.off()


################
#	RealPrior2
################
load("~/Desktop/Dropbox/space.mix/data/warblers/warbler_spacemix/pop/warbler_pop_spaceruns/real_prior2/warb_pop_spaceruns_realpr2_LongRun/warb_pop_spaceruns_realpr2space_MCMC_output1.Robj")
	best <- which.max(Prob)
	target.coords <- procrusteez(warbler.pop.coords,population.coordinates[[best]][1:k,],k,option=1)
	source.coords <- procrusteez(warbler.pop.coords,population.coordinates[[best]][1:k,],k,source.locs=population.coordinates[[best]][(k+1):(2*k),],option=2)
	pop.plot.cols <- fade.admixture.source.points(pop.col,admix.proportions[,best])

	png(file="~/Desktop/Dropbox/space.mix/ms/population_warbler_map_realpr2.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(53,101), #realpr1: c(68,115), realpr2: c(53,101), randpr: c(71,114)
					ylim=c(25,54), #realpr1: c(26,53), realpr2: c(25,54), randpr: c(71,114)
					xlab="long",
					ylab="lat")
				text(target.coords[c(1:k),],
						labels=pops,
						col=pop.col,
						font=2,cex=0.9)
				points(source.coords,
							col=pop.plot.cols,
							pch=20)
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col= pop.plot.cols,
					lwd=admix.proportions[,best],#last.params$admix.proportions,
					length=0.1)
				box(lwd=2)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/population_warbler_map_no_arrows_realpr2.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(68,101), #realpr1: c(68,101), realpr2: c(68,101), randpr: c(71,99)
					ylim=c(25,55), #realpr1: c(25,55), realpr2: c(25,55), randpr: c(20.5,56)
					xlab="long",
					ylab="lat")
				text(target.coords[c(1:k),],
						labels=pops,
						col=pop.col,
						font=2,cex=0.9)
				box(lwd=2)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/population_warbler_admix_values_nugget_realpr2.png",res=300,width=4.5*300,height=5*300,pointsize=9)
		#quartz(width=5,height=5,pointsize=9)
		order <- c(1,3:15,2,16:22)
		par(mfrow=c(2,1),mar=c(1,4,0,1),oma=c(1,1,1,1))
			plot(admix.proportions[order,best]/2,type='n',
					xlab="",
					xaxt='n',
					ylab="estimated adxmiture proportion")
				text(admix.proportions[order,best]/2,
						labels=pops[order],
						col=pop.col[order],
						font=2,cex=0.9)
				box(lwd=2)
			plot(nugget[order,best],type='n',
					xaxt='n',
					ylab="estimated population nugget")
				text(nugget[order,best],
						labels=pops[order],
						col=pop.col[order],
						font=2,cex=0.9)
				box(lwd=2)
	dev.off()

################
#	RandPrior1
################
load("~/Desktop/Dropbox/space.mix/data/warblers/warbler_spacemix/pop/warbler_pop_spaceruns/rand_prior1/warb_pop_spaceruns_randpr1_LongRun/warb_pop_spaceruns_randpr1space_MCMC_output1.Robj")
	best <- which.max(Prob)
	target.coords <- procrusteez(warbler.pop.coords,population.coordinates[[best]][1:k,],k,option=1)
	source.coords <- procrusteez(warbler.pop.coords,population.coordinates[[best]][1:k,],k,source.locs=population.coordinates[[best]][(k+1):(2*k),],option=2)
	pop.plot.cols <- fade.admixture.source.points(pop.col,admix.proportions[,best])

	png(file="~/Desktop/Dropbox/space.mix/ms/population_warbler_map_randpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(71,114), #realpr1: c(68,115), realpr2: c(53,101), randpr: c(71,114)
					ylim=c(23,55), #realpr1: c(26,53), realpr2: c(25,54), randpr: c(23,55)
					xlab="long",
					ylab="lat")
				text(target.coords[c(1:k),],
						labels=pops,
						col=pop.col,
						font=2,cex=0.9)
				points(source.coords,
							col=pop.plot.cols,
							pch=20)
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col= pop.plot.cols,
					lwd=admix.proportions[,best],#last.params$admix.proportions,
					length=0.1)
				box(lwd=2)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/population_warbler_map_no_arrows_randpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(71,99), #realpr1: c(68,101), realpr2: c(68,101), randpr: c(71,99)
					ylim=c(20.5,56), #realpr1: c(25,55), realpr2: c(25,55), randpr: c(20.5,56)
					xlab="long",
					ylab="lat")
				text(target.coords[c(1:k),],
						labels=pops,
						col=pop.col,
						font=2,cex=0.9)
				box(lwd=2)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/population_warbler_admix_values_nugget_randpr1.png",res=300,width=4.5*300,height=5*300,pointsize=9)
		#quartz(width=5,height=5,pointsize=9)
		order <- c(1,3:15,2,16:22)
		par(mfrow=c(2,1),mar=c(1,4,0,1),oma=c(1,1,1,1))
			plot(admix.proportions[order,best]/2,type='n',
					xlab="",
					xaxt='n',
					ylab="estimated adxmiture proportion")
				text(admix.proportions[order,best]/2,
						labels=pops[order],
						col=pop.col[order],
						font=2,cex=0.9)
				box(lwd=2)
			plot(nugget[order,best],type='n',
					xaxt='n',
					ylab="estimated population nugget")
				text(nugget[order,best],
						labels=pops[order],
						col=pop.col[order],
						font=2,cex=0.9)
				box(lwd=2)
	dev.off()
	
################################
#	WARBLER IND FIGS
################################

################################
#	GLOBE FIGS
################################

################################
#	SIM FIGS
################################
