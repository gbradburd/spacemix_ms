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

fade.admixture.source.points <- function(pop.cols,admix.proportions){
	faded.colors <- numeric(length(pop.cols))
	for(i in 1:length(pop.cols)){
		faded.colors[i] <- adjustcolor(pop.cols[i],admix.proportions[i])
	}
	return(faded.colors)
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

get.mean.sample.size <- function(sample.sizes){
	mean.sample.size <- mean(sample.sizes[which(sample.sizes!=0)])
	return(mean.sample.size)
}

get.transformation.matrix <- function(mean.sample.sizes){
	k <- length(mean.sample.sizes)
	transformation.matrix <- diag(k) - matrix(mean.sample.sizes/(sum(mean.sample.sizes)),nrow=k,ncol=k,byrow=TRUE)
	return(transformation.matrix)
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
#	PCA Map
################

cov_hat <- cov(t(warbler.pop.allele.counts/warbler.pop.sample.sizes))
mean.centering.matrix <- get.transformation.matrix(apply(warbler.pop.sample.sizes,1,get.mean.sample.size))
mc.cov_hat <- mean.centering.matrix %*% cov_hat %*% t(mean.centering.matrix)
pc.coords <- cbind(eigen(mc.cov_hat)$vectors[,1],eigen(mc.cov_hat)$vectors[,2])
proc.pc.coords <- fitted(vegan::procrustes(warbler.pop.coords,pc.coords))

png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/warb_pop_PC_map.png",res=200,width=5*200,height=4*200)
	#quartz(width=5,height=4)
	par(mar=c(4.5,4.5,3,1))
	plot(proc.pc.coords,type='n',
		xlab=paste("PC1 (",
					round(eigen(mc.cov_hat)$values[c(1)]/sum(eigen(mc.cov_hat)$values),3) * 100,
				"%)",sep=""),
		ylab=paste("PC2 (",
					round(eigen(mc.cov_hat)$values[c(2)]/sum(eigen(mc.cov_hat)$values),3) * 100,
				"%)",sep=""),
		main="Warbler Population PC Map",cex.axis=0.7)
	text(proc.pc.coords,col=pop.col,labels=pops,cex=0.7,font=2)
	box(lwd=2)
dev.off()
################
#	RealPrior1
################
load("~/Desktop/Dropbox/space.mix/data/warblers/warbler_spacemix/pop/warbler_pop_spaceruns/real_prior1/warb_pop_spaceruns_realpr1_LongRun/warb_pop_spaceruns_realpr1space_MCMC_output1.Robj")
	best <- which.max(Prob)
	target.coords <- procrusteez(warbler.pop.coords,population.coordinates[[best]][1:k,],k,option=1)
	source.coords <- procrusteez(warbler.pop.coords,population.coordinates[[best]][1:k,],k,source.locs=population.coordinates[[best]][(k+1):(2*k),],option=2)
	pop.plot.cols <- fade.admixture.source.points(pop.col,admix.proportions[,best])

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/population_warbler_map_realpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(68,116), #realpr1: c(68,115), realpr2: c(53,101), randpr: c(71,114)
					ylim=c(26,53), #realpr1: c(26,53), realpr2: c(25,54), randpr: c(71,114)
					xlab="Eastings",
					ylab="Northings")
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/population_warbler_map_no_arrows_realpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(68,101), #realpr1: c(68,101), realpr2: c(68,101), randpr: c(71,99)
					ylim=c(25,55), #realpr1: c(25,55), realpr2: c(25,55), randpr: c(20.5,56)
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=pop.col,
						font=2,cex=0.9)
				box(lwd=2)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/population_warbler_admix_values_nugget_realpr1.png",res=300,width=4.5*300,height=5*300,pointsize=9)
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/population_warbler_map_realpr2.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(53,101), #realpr1: c(68,115), realpr2: c(53,101), randpr: c(71,114)
					ylim=c(25,54), #realpr1: c(26,53), realpr2: c(25,54), randpr: c(71,114)
					xlab="Eastings",
					ylab="Northings")
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/population_warbler_map_no_arrows_realpr2.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(68,101), #realpr1: c(68,101), realpr2: c(68,101), randpr: c(71,99)
					ylim=c(25,55), #realpr1: c(25,55), realpr2: c(25,55), randpr: c(20.5,56)
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=pop.col,
						font=2,cex=0.9)
				box(lwd=2)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/population_warbler_admix_values_nugget_realpr2.png",res=300,width=4.5*300,height=5*300,pointsize=9)
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/population_warbler_map_randpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(71,114), #realpr1: c(68,115), realpr2: c(53,101), randpr: c(71,114)
					ylim=c(23,55), #realpr1: c(26,53), realpr2: c(25,54), randpr: c(23,55)
					xlab="Eastings",
					ylab="Northings")
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/population_warbler_map_no_arrows_randpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(71,99), #realpr1: c(68,101), realpr2: c(68,101), randpr: c(71,99)
					ylim=c(20.5,56), #realpr1: c(25,55), realpr2: c(25,55), randpr: c(20.5,56)
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=pop.col,
						font=2,cex=0.9)
				box(lwd=2)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/population_warbler_admix_values_nugget_randpr1.png",res=300,width=4.5*300,height=5*300,pointsize=9)
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
load("~/Desktop/Dropbox/space.mix/data/warblers/warbler_spacemix/ind/warb_ind_spaceruns/real_prior1/warbler_ind_dataset.Robj")
	inds <- row.names(warbler.ind.allele.counts)
	k <- length(inds)
		inds[11] <- "Vir-STvi1"
		inds[12] <- "Vir-STvi2"
		inds[13] <- "Vir-STvi3"
	ind.subspp <- unlist(strsplit(inds,"-"))[seq(1,190,2)]
	inds.col <- numeric(length(ind.subspp))
		inds.col[grepl("Vir",ind.subspp)] <- "dodgerblue2"
		inds.col[grepl("Ni",ind.subspp)] <- "slateblue4"
		inds.col[grepl("Lud",ind.subspp)] <- "mediumseagreen"
		inds.col[grepl("Tro",ind.subspp)] <- "gold"
		inds.col[grepl("Obs",ind.subspp)] <- "orange"
		inds.col[grepl("Plu",ind.subspp)] <- "red"

		plot.inds <- gsub(" ","",inds)
		plot.inds <- gsub("[[:digit:]]","",plot.inds)
		plot.inds <- gsub(c("Plu-"),"",plot.inds)
		plot.inds <- gsub(c("Vir-"),"",plot.inds)
		plot.inds <- gsub(c("Ni-"),"",plot.inds)
		plot.inds <- gsub(c("Lud-"),"",plot.inds)
		plot.inds <- gsub(c("Tro-"),"",plot.inds)
		plot.inds <- gsub(c("Obs-"),"",plot.inds)
		plot.inds <- gsub(c("vi"),"",plot.inds)

################
#	PCA Map
################

cov_hat <- cov(t(warbler.ind.allele.counts/warbler.ind.sample.sizes))
mean.centering.matrix <- get.transformation.matrix(apply(warbler.ind.sample.sizes,1,get.mean.sample.size))
mc.cov_hat <- mean.centering.matrix %*% cov_hat %*% t(mean.centering.matrix)
pc.coords <- cbind(eigen(mc.cov_hat)$vectors[,1],eigen(mc.cov_hat)$vectors[,2])
proc.pc.coords <- fitted(vegan::procrustes(warbler.ind.coords,pc.coords))

png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/warb_ind_PC_map.png",res=200,width=5*200,height=4*200)
	#quartz(width=5,height=4)
	par(mar=c(4.5,4.5,3,1))
	plot(proc.pc.coords,type='n',
		xlab=paste("PC1 (",
					round(eigen(mc.cov_hat)$values[c(1)]/sum(eigen(mc.cov_hat)$values),3) * 100,
				"%)",sep=""),
		ylab=paste("PC2 (",
					round(eigen(mc.cov_hat)$values[c(2)]/sum(eigen(mc.cov_hat)$values),3) * 100,
				"%)",sep=""),
		main="Warbler Individual PC Map",cex.axis=0.7)
	text(proc.pc.coords,col=adjustcolor(inds.col,0.8),labels=plot.inds,cex=0.5,font=2)
	box(lwd=2)
dev.off()

################
#	RealPrior1
################
load("~/Desktop/Dropbox/space.mix/data/warblers/warbler_spacemix/ind/warb_ind_spaceruns/real_prior1/warb_ind_spaceruns_realpr1_LongRun/warb_ind_spaceruns_realpr1space_MCMC_output1.Robj")
best <- which.max(Prob)
	target.coords <- procrusteez(warbler.ind.coords,population.coordinates[[best]][1:k,],k,option=1)
	source.coords <- procrusteez(warbler.ind.coords,population.coordinates[[best]][1:k,],k,source.locs=population.coordinates[[best]][(k+1):(2*k),],option=2)
	ind.plot.cols <- fade.admixture.source.points(inds.col,admix.proportions[,best])
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/individual_warbler_map_arrows_amped_admixture_realpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1],source.coords[,1]),
							max(target.coords[,1],source.coords[,1])),
					ylim=c(min(target.coords[,2],source.coords[,2]),
							max(target.coords[,2],source.coords[,2])),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=plot.inds,
						col=adjustcolor(inds.col,0.8),
						font=2,cex=0.9)
				points(source.coords,
							col=fade.admixture.source.points(inds.col,admix.proportions[,best]*10),
							pch=20)
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col= inds.col,
					lwd=admix.proportions[,best]*10,
					length=0.1)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/individual_warbler_map_arrows_realpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1],source.coords[,1]),
							max(target.coords[,1],source.coords[,1])),
					ylim=c(min(target.coords[,2],source.coords[,2]),
							max(target.coords[,2],source.coords[,2])),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=plot.inds,
						col=adjustcolor(inds.col,0.8),
						font=2,cex=0.9)
				points(source.coords,
							col=ind.plot.cols,
							pch=20)
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col= inds.col,
					lwd=admix.proportions[,best],
					length=0.1)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/individual_warbler_map_noarrows_realpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1]),max(target.coords[,1])+2),
					ylim=c(min(target.coords[,2]),max(target.coords[,2])),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=plot.inds,
						col=adjustcolor(inds.col,0.8),
						font=2,cex=0.9)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/individual_warbler_map_noarrows_closeup_realpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[-c(which(plot.inds=="TU")),1]),max(target.coords[-c(which(plot.inds=="TU")),1])),
					ylim=c(min(target.coords[-c(which(plot.inds=="TU")),2]),max(target.coords[-c(which(plot.inds=="TU")),2])),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[-c(which(plot.inds=="TU")),],
						labels=plot.inds[-c(which(plot.inds=="TU"))],
						col=adjustcolor(inds.col[-c(which(plot.inds=="TU"))],0.8),
						font=2,cex=0.9)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/individual_warbler_map_noarrows_closeup_nugget_realpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[-c(which(plot.inds=="TU")),1]),max(target.coords[-c(which(plot.inds=="TU")),1])),
					ylim=c(min(target.coords[-c(which(plot.inds=="TU")),2]),max(target.coords[-c(which(plot.inds=="TU")),2])),
					xlab="Eastings",
					ylab="Northings")
				points(target.coords[-c(which(plot.inds=="TU")),],
						col=adjustcolor(inds.col[-c(which(plot.inds=="TU"))],0.8),
						cex=last.params$nugget*5,lwd=1.5)
	dev.off()
################
#	RealPrior2
################
load("~/Desktop/Dropbox/space.mix/data/warblers/warbler_spacemix/ind/warb_ind_spaceruns/real_prior2/warb_ind_spaceruns_realpr2_LongRun/warb_ind_spaceruns_realpr2space_MCMC_output1.Robj")

best <- which.max(Prob)
	target.coords <- procrusteez(warbler.ind.coords,population.coordinates[[best]][1:k,],k,option=1)
	source.coords <- procrusteez(warbler.ind.coords,population.coordinates[[best]][1:k,],k,source.locs=population.coordinates[[best]][(k+1):(2*k),],option=2)
	ind.plot.cols <- fade.admixture.source.points(inds.col,admix.proportions[,best])
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/individual_warbler_map_arrows_amped_admixture_realpr2.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1],source.coords[,1]),
							max(target.coords[,1],source.coords[,1])),
					ylim=c(min(target.coords[,2],source.coords[,2]),
							max(target.coords[,2],source.coords[,2])),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=plot.inds,
						col=adjustcolor(inds.col,0.8),
						font=2,cex=0.9)
				points(source.coords,
							col=fade.admixture.source.points(inds.col,admix.proportions[,best]*10),
							pch=20)
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col= inds.col,
					lwd=admix.proportions[,best]*10,
					length=0.1)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/individual_warbler_map_arrows_realpr2.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1],source.coords[,1]),
							max(target.coords[,1],source.coords[,1])),
					ylim=c(min(target.coords[,2],source.coords[,2]),
							max(target.coords[,2],source.coords[,2])),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=plot.inds,
						col=adjustcolor(inds.col,0.8),
						font=2,cex=0.9)
				points(source.coords,
							col=ind.plot.cols,
							pch=20)
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col= inds.col,
					lwd=admix.proportions[,best],
					length=0.1)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/individual_warbler_map_noarrows_realpr2.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1]),max(target.coords[,1])+2),
					ylim=c(min(target.coords[,2]),max(target.coords[,2])),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=plot.inds,
						col=adjustcolor(inds.col,0.8),
						font=2,cex=0.9)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/individual_warbler_map_noarrows_closeup_realpr2.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[-c(which(plot.inds=="TU")),1]),max(target.coords[-c(which(plot.inds=="TU")),1])),
					ylim=c(min(target.coords[-c(which(plot.inds=="TU")),2]),max(target.coords[-c(which(plot.inds=="TU")),2])),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[-c(which(plot.inds=="TU")),],
						labels=plot.inds[-c(which(plot.inds=="TU"))],
						col=adjustcolor(inds.col[-c(which(plot.inds=="TU"))],0.8),
						font=2,cex=0.9)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/individual_warbler_map_noarrows_closeup_nugget_realpr2.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[-c(which(plot.inds=="TU")),1]),max(target.coords[-c(which(plot.inds=="TU")),1])),
					ylim=c(min(target.coords[-c(which(plot.inds=="TU")),2]),max(target.coords[-c(which(plot.inds=="TU")),2])),
					xlab="Eastings",
					ylab="Northings")
				points(target.coords[-c(which(plot.inds=="TU")),],
						col=adjustcolor(inds.col[-c(which(plot.inds=="TU"))],0.8),
						cex=last.params$nugget*5,lwd=1.5)
	dev.off()

################
#	RandPrior2
################
load("~/Desktop/Dropbox/space.mix/data/warblers/warbler_spacemix/ind/warb_ind_spaceruns/rand_prior2/warb_ind_spaceruns_randpr2_LongRun/warb_ind_spaceruns_randpr2space_MCMC_output1.Robj")

best <- which.max(Prob)
	target.coords <- procrusteez(warbler.ind.coords,population.coordinates[[best]][1:k,],k,option=1)
	source.coords <- procrusteez(warbler.ind.coords,population.coordinates[[best]][1:k,],k,source.locs=population.coordinates[[best]][(k+1):(2*k),],option=2)
	ind.plot.cols <- fade.admixture.source.points(inds.col,admix.proportions[,best])
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/individual_warbler_map_arrows_amped_admixture_randpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1],source.coords[,1]),
							max(target.coords[,1],source.coords[,1])),
					ylim=c(min(target.coords[,2],source.coords[,2]),
							max(target.coords[,2],source.coords[,2])),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=plot.inds,
						col=adjustcolor(inds.col,0.8),
						font=2,cex=0.9)
				points(source.coords,
							col=fade.admixture.source.points(inds.col,admix.proportions[,best]*10),
							pch=20)
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col= inds.col,
					lwd=admix.proportions[,best]*10,
					length=0.1)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/individual_warbler_map_arrows_randpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1],source.coords[,1]),
							max(target.coords[,1],source.coords[,1])),
					ylim=c(min(target.coords[,2],source.coords[,2]),
							max(target.coords[,2],source.coords[,2])),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=plot.inds,
						col=adjustcolor(inds.col,0.8),
						font=2,cex=0.9)
				points(source.coords,
							col=ind.plot.cols,
							pch=20)
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col= inds.col,
					lwd=admix.proportions[,best],
					length=0.1)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/individual_warbler_map_noarrows_randpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1]),max(target.coords[,1])+2),
					ylim=c(min(target.coords[,2]),max(target.coords[,2])),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=plot.inds,
						col=adjustcolor(inds.col,0.8),
						font=2,cex=0.9)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/individual_warbler_map_noarrows_closeup_randpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[-c(which(plot.inds=="TU")),1]),max(target.coords[-c(which(plot.inds=="TU")),1])),
					ylim=c(min(target.coords[-c(which(plot.inds=="TU")),2]),max(target.coords[-c(which(plot.inds=="TU")),2])),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[-c(which(plot.inds=="TU")),],
						labels=plot.inds[-c(which(plot.inds=="TU"))],
						col=adjustcolor(inds.col[-c(which(plot.inds=="TU"))],0.8),
						font=2,cex=0.9)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/warblers/individual_warbler_map_noarrows_closeup_nugget_randpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[-c(which(plot.inds=="TU")),1]),max(target.coords[-c(which(plot.inds=="TU")),1])),
					ylim=c(min(target.coords[-c(which(plot.inds=="TU")),2]),max(target.coords[-c(which(plot.inds=="TU")),2])),
					xlab="Eastings",
					ylab="Northings")
				points(target.coords[-c(which(plot.inds=="TU")),],
						col=adjustcolor(inds.col[-c(which(plot.inds=="TU"))],0.8),
						cex=last.params$nugget*5,lwd=1.5)
	dev.off()


################################
#	GLOBE FIGS
################################

################
#	CYOL - rand prior
################
load("~/Desktop/Dropbox/space.mix/data/globetrotter/globe_spacemix/globe_spaceruns/globe_no_admixture/rand_prior1/globe_spaceruns_NoAd_randpr1_LongRun/globe_spaceruns_NoAd_randpr1space_MCMC_output1.Robj")
load("~/Desktop/Dropbox/space.mix/data/globetrotter/globe_spacemix/globe_spaceruns/rand_prior2/globetrotter_dataset.Robj")

	globe.coords <- cbind(globetrotter.long, globetrotter.lat)
	pops <- row.names(globetrotter.counts)
	k <- length(pops)

		continent.col <- numeric(k)
			continent.col[which(globetrotter.long < -50)] <- "orange"
			continent.col[match(c("BantuKenya","BantuSouthAfrica","BiakaPygmy",
									"Egyptian","Ethiopian","EthiopianJew","Hadza","Mandenka",
									"MbutiPygmy","Moroccan","Mozabite","Sandawe","SanNamibia",
									"SanKhomani","Tunisian","Yoruba"),pops)] <- "forestgreen"
		continent.col[which(globetrotter.long > 100 &
							globetrotter.lat < 5)] <- "brown"
		continent.col[which(continent.col==0)] <- rainbow(length(continent.col[which(continent.col==0)]),
															start=4/6,end=6/6)[as.numeric(cut(globetrotter.long[which(continent.col==0)],length(which(continent.col==0))))]
		americas <- which(continent.col=="orange")
		africa <- which(continent.col=="forestgreen")
		oceania <- which(continent.col=="brown")
		east.asia <- which(globetrotter.long > 95 & 
								globetrotter.lat > 11.5)
		western.eurasia <- c(1:k)[-c(americas,africa,oceania,east.asia)]
		eurasia <- c(1:k)[-c(americas,africa,oceania)]

		# af.eff.lat <- (globetrotter.lat[africa] + abs(min(globetrotter.lat[africa])))/max(globetrotter.lat[africa] + abs(min(globetrotter.lat[africa])))
		af.loc.cols <- hsv(h = seq(0.22,0.69,length.out=length(africa))[rank(globetrotter.lat[africa])],
				s = 1,
				v = 1)
		adj.nonamaf.long <- globetrotter.long[-c(africa,americas)] + abs(min(globe.coords[western.eurasia,1]))
		eur.eff.long <- adj.nonamaf.long/max(adj.nonamaf.long)
		eur.loc.cols <- hsv(h = eur.eff.long * 0.4 + 0.6,s=1,v=1)
		am.eff.long <- (globetrotter.long[americas] + abs(min(globetrotter.long[americas])))/max(globetrotter.long[americas] + abs(min(globetrotter.long[americas])))
		am.loc.cols <- hsv(h = (am.eff.long) * 0.08 + 0.03,s=1,v=1)
		continent.col <- numeric(k)
		continent.col[americas] <- am.loc.cols
		continent.col[africa] <- af.loc.cols
		continent.col[-c(africa,americas)] <- eur.loc.cols

	best <- which.max(Prob)
	target.coords <- procrusteez(globe.coords,population.coordinates[[best]][1:k,],k,option=1)
	source.coords <- procrusteez(globe.coords,population.coordinates[[best]][1:k,],k,source.locs=population.coordinates[[best]][(k+1):(2*k),],option=2)
	eurasia.target.coords <- fitted(procrustes(globe.coords[eurasia,],population.coordinates[[best]][1:k,][eurasia,]))
	eurasia.source.coords <- fitted(procrustes(globe.coords[eurasia,],population.coordinates[[best]][(k+1):(2*k),][eurasia,]))
	require(maps)
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/globe_world_map_dots.png",res=300,width=9*300,height=5.5*300)
	#quartz(width=9,height=5.5)
	map("world",interior=FALSE,lwd=0.5,ylim=c(-60,85))
		box(lwd=2)
		points(globe.coords,pch=20,col=continent.col,cex=2)
#			legend(x = -175,y=-10,
#					legend = c("Africa","Western Eurasia","Central Eurasia","Eastern Eurasia","Oceania","Americas"),
#					text.col = c("forestgreen","blue","purple","red","brown","orange"),cex=0.8)
	dev.off()
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/globe_world_map_text.png",res=300,width=9*300,height=5.5*300)
	#quartz(width=9,height=5.5)
	map("world",interior=FALSE,lwd=0.5,ylim=c(-60,85),xlim=c(-130,180))
		box(lwd=2)
		text(globe.coords,pops,col=continent.col,cex=0.8,font=2)
#			legend(x = -175,y=-10,
#					legend = c("Africa","Western Eurasia","Central Eurasia","Eastern Eurasia","Oceania","Americas"),
#					text.col = c("forestgreen","blue","purple","red","brown","orange"),cex=0.8)
	dev.off()
	
	x.min <- min(target.coords[,1]) - 5
	x.max <- max(target.coords[,1]) + 5
	y.min <- min(target.coords[,2]) - 5
	y.max <- max(target.coords[,2]) + 5
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/globe_NoAd_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(x.min,x.max),
					ylim=c(y.min,y.max),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,0.8),
						font=2,cex=0.8)
			box(lwd=2)
#			legend(x = "bottomright",pch=NA,
#					legend = c("Africa","Western Eurasia","Central Eurasia","Eastern Eurasia","Oceania","Americas"),
#					text.col = c("forestgreen","blue","purple","red","brown","orange"))
	dev.off()
	
	x.min <- min(target.coords[africa,1]) - 5
	x.max <- max(target.coords[africa,1]) + 5
	y.min <- min(target.coords[africa,2]) - 5
	y.max <- max(target.coords[africa,2]) + 5

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/globe_Africa_NoAd_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords[africa,],type='n',
					xlim=c(x.min,x.max),
					ylim=c(y.min,y.max),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[africa,],
						labels=pops[africa],
						col=adjustcolor(continent.col[africa],0.8),
						font=2,cex=0.8)
			box(lwd=2)
	dev.off()

	x.min <- min(target.coords[-africa,1]) - 5
	x.max <- max(target.coords[-africa,1]) + 5
	y.min <- min(target.coords[-africa,2]) - 5
	y.max <- max(target.coords[-africa,2]) + 5

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/globe_Eurasia_Americas_Oceania_NoAd_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords[-africa,],type='n',
					xlim=c(x.min,x.max),
					ylim=c(y.min,y.max),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[-africa,],
						labels=pops[-africa],
						col=adjustcolor(continent.col[-africa],0.8),
						font=2,cex=0.8)
			box(lwd=2)
	dev.off()

	x.min <- min(target.coords[-c(africa,americas,oceania),1]) - 2
	x.max <- max(target.coords[-c(africa,americas,oceania),1]) + 2
	y.min <- min(target.coords[-c(africa,americas,oceania),2]) - 2
	y.max <- max(target.coords[-c(africa,americas,oceania),2]) + 2

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/globe_Eurasia_NoAd_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords[-c(africa,americas,oceania),],type='n',
					xlim=c(x.min,x.max),
					ylim=c(y.min,y.max),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[-c(africa,americas,oceania),],
						labels=pops[-c(africa,americas,oceania)],
						col=adjustcolor(continent.col[-c(africa,americas,oceania)],0.8),
						font=2,cex=0.8)
			box(lwd=2)
	dev.off()

	x.min <- min(eurasia.target.coords[,1]) - 2
	x.max <- max(eurasia.target.coords[,1]) + 2
	y.min <- min(eurasia.target.coords[,2]) - 2
	y.max <- max(eurasia.target.coords[,2]) + 2
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/globe_Eurasia_NoAd_map_indproc.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(eurasia.target.coords,type='n',
					xlim=c(x.min,x.max),
					ylim=c(y.min,y.max),
					xlab="Eastings",
					ylab="Northings")
				text(eurasia.target.coords,
						labels=pops[eurasia],
						col=adjustcolor(continent.col[eurasia],0.8),
						font=2,cex=0.8)
			box(lwd=2)
	dev.off()

	clusters <- list(	western.eurasia = western.eurasia,
						americas = americas,
						africa = africa,
						oceania = oceania,
						east.asia = east.asia)
	cluster.cols <- c("purple4","orange","forestgreen","brown","red")
	obs.D <- fields::rdist.earth(globe.coords)
	par.D <- fields::rdist.earth(target.coords)

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/globe_NoAd_dist_compare.png",res=200,height=5*200,width=12*200)
		#quartz(height=5,width=12)
		par(mfrow=c(1,2),mar=c(4,4,2,2))
		plot(obs.D,par.D,col="gray",pch=20,
#			xlim=c(0,5500),ylim=c(0,4500),
			ylab="estimated distance",
			xlab="observed distance",cex=0.7)
			for(i in 1:length(clusters)){
				use.these <- clusters[[i]]
				#points(obs.D[use.these,use.these],par.D[use.these,use.these],col=1,pch=20,cex=0.73)
				#points(obs.D[use.these,use.these],par.D[use.these,use.these],col=cluster.cols[i],pch=20,cex=0.7)
				lines(lowess(par.D[1:k,1:k][use.these,use.these] ~ obs.D[use.these,use.these]),col=cluster.cols[i],lwd=2)
			}
			rect(xleft = 0,
				 ybottom = 0,
				 xright = 5500,
				 ytop = 4500,
				 lty = 2)
			legend(x = "topleft",pch=NA,
					legend = c("Africa","Western Eurasia","East Asia","Oceania","Americas"),
					text.col = c("forestgreen","purple4","red","brown","orange"),
					cex=0.8)
		box(lwd=2)
		plot(obs.D,par.D,col="gray",pch=20,
			xlim=c(0,5500),ylim=c(0,4500),
			ylab="",
			xlab="observed distance",cex=0.7)
			for(i in 1:length(clusters)){
				use.these <- clusters[[i]]
				points(obs.D[use.these,use.these],par.D[use.these,use.these],col=1,pch=20,cex=0.73)
				points(obs.D[use.these,use.these],par.D[use.these,use.these],col=cluster.cols[i],pch=20,cex=0.7)
				lines(lowess(par.D[1:k,1:k][use.these,use.these] ~ obs.D[use.these,use.these]),col=cluster.cols[i],lwd=2)
			}
		box(lwd=2)
	dev.off()

	line.obs <- vector("list",length=length(clusters))
	line.coeffs <- numeric(length(clusters))
	centroids <- matrix(0,nrow=length(clusters),ncol=2)
	africa.centroid <- matrix(colMeans(globe.coords[clusters$africa,]),ncol=2,nrow=1)
	mean.dist.from.africa <- numeric(length(clusters))
		for(i in 1:length(line.obs)){
			use.these <- clusters[[i]]
			line.obs[[i]] <- lm(par.D[1:k,1:k][use.these,use.these][upper.tri(par.D[1:k,1:k][use.these,use.these],diag=TRUE)] ~ 
									obs.D[use.these,use.these][upper.tri(obs.D[1:k,1:k][use.these,use.these],diag=TRUE)])
			line.coeffs[i] <- line.obs[[i]]$coefficients[2]
			centroids[i,] <- colMeans(globe.coords[clusters[[i]],])
				row.names(centroids)[i] <- names(clusters[[i]])
			mean.dist.from.africa[i] <- fields::rdist.earth(africa.centroid,centroids[i,,drop=FALSE])
		}
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/globe_NoAd_dist_decay.png",res=200,height=5*200,width=12*200)
		#quartz(height=5,width=12)
		par(mfrow=c(1,2),mar=c(4,4,2,2))
		plot(obs.D,par.D,col="gray",pch=20,
#			xlim=c(0,5500),ylim=c(0,4500),
			ylab="estimated distance",
			xlab="observed distance",cex=0.7)
			for(i in 1:length(clusters)){
				use.these <- clusters[[i]]
				abline(line.obs[[i]],col=cluster.cols[i],lwd=2)
				points(obs.D[use.these,use.these],par.D[use.these,use.these],col=1,pch=20,cex=0.73)
				points(obs.D[use.these,use.these],par.D[use.these,use.these],col=cluster.cols[i],pch=20,cex=0.7)
			}
		box(lwd=2)
		plot(mean.dist.from.africa,line.coeffs,col=cluster.cols,pch=19,cex=3,
			ylab="slope of observed vs. estimated distance",
			xlab="observed distance from Africa")
			legend(x = "topright",pch=NA,
					legend = c("Africa","Western Eurasia","East Asia","Oceania","Americas"),
					text.col = c("forestgreen","purple4","red","brown","orange"),
					cex=1)
		box(lwd=2)
	dev.off()

################
#	Admixture - rand prior
################
load("~/Desktop/Dropbox/space.mix/data/globetrotter/globe_spacemix/globe_spaceruns/rand_prior2/globe_spaceruns_randpr1_LongRun/globe_spaceruns_randpr1space_MCMC_output1.Robj")
	best <- which.max(Prob)
	target.coords <- procrusteez(globe.coords,population.coordinates[[best]][1:k,],k,option=1)
	source.coords <- procrusteez(globe.coords,population.coordinates[[best]][1:k,],k,source.locs=population.coordinates[[best]][(k+1):(2*k),],option=2) 

	globe.admix.plot.cols <- fade.admixture.source.points(continent.col,admix.proportions[,best])

MAP.admix.props <- admix.proportions[,best]
MAP.nugget <- nugget[,best]

globe_ad_obj <- list(africa = africa,americas = americas,continent.col = continent.col,
						east.asia = east.asia,globe.admix.plot.cols = globe.admix.plot.cols,
						globe.coords = globe.coords,k = k,MAP.admix.props = MAP.admix.props,
						MAP.nugget = MAP.nugget,oceania = oceania,pops = pops,
						source.coords = source.coords,target.coords = target.coords,
						western.eurasia = western.eurasia)

save(globe_ad_obj,file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/globe_Ad_object.Robj")

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/globe_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(18,68),
					ylim = c(-20,50),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
				# points(source.coords[,1],
						# source.coords[,2],
							# col=globe.admix.plot.cols,
							# pch=20)
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			box(lwd=2)
#			legend(x = "bottomright",pch=NA,
#					legend = c("Africa","Western Eurasia","Central Eurasia","Eastern Eurasia","Oceania","Americas","population","source of admixture"),
#					text.col = c("forestgreen","blue","purple","red","brown","orange",1,1),
#					text.font = c(rep(1,6),2,3),family=c())
			# legend(x="topleft",
					# lwd = c(1,0.5,0.1),
					# col = c(adjustcolor(1,1),adjustcolor(1,0.5),adjustcolor(1,0.1)),
					# legend = c("w = 0.5","w = 0.25","w = 0.05"),
					# title = "Admixture proportions")
	dev.off()


	require(maps)
	x.subplot1 <- c(43,68)
	y.subplot1 <- c(-18,20.5)
	x.subplot2 <- c(15,32.5)
	y.subplot2 <- c(25.5,44)
	subplot1.x.coords <- c(19.2,30)
	subplot1.y.coords <- c(-19,-5.5)
	subplot2.x.coords <- c(-165,180)
	subplot2.y.coords <- c(-60,80)
	#rect(xleft = x.subplot2[1],ybottom = y.subplot2[1],xright = x.subplot2[2],ytop = y.subplot2[2])
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/globe_Ad_map_AfricaInset.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(15,68),
					ylim = c(-20,50),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
				# points(source.coords[,1],
						# source.coords[,2],
							# col=globe.admix.plot.cols,
							# pch=20)
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
		rect(xleft = subplot1.x.coords[1]-2,ybottom = subplot1.y.coords[1],
				xright = subplot1.x.coords[2]+1.5,ytop = subplot1.y.coords[2]+1,lty=2,border="gray",col=NA,lwd=0.6)
		lines(x = c(subplot1.x.coords[2]+1.5,x.subplot1[1]), y = c(subplot1.y.coords[2]+1,y.subplot1[2]),lwd=0.6 , lty=2, col="gray")
		lines(x = c(subplot1.x.coords[2]+1.5,x.subplot1[1]), y = c(subplot1.y.coords[1],y.subplot1[1]),lwd=0.6 , lty=2, col="gray")
		TeachingDemos::subplot(fun = {
						plot(target.coords,type='n',xlim = subplot1.x.coords,ylim = subplot1.y.coords,xlab="",ylab="",xaxt='n',yaxt='n') ; 
						text(target.coords[c(1:k),],
								labels=pops,
								col=adjustcolor(continent.col,1),
								font=2,cex=1) ; 
						text(source.coords[,1],
								source.coords[,2],
								labels=pops,
								font=3,
								col=globe.admix.plot.cols,cex=1,family="HersheySerif") ; 
						arrows(	x0 = source.coords[,1],
								y0 = source.coords[,2],
								x1 = target.coords[,1],
								y1 = target.coords[,2],
								col=globe.admix.plot.cols,
								lwd=admix.proportions[,best],
								length=0.1) ; 
							abline(v=0,lty=2,lwd=0.5) ; 
							box(lwd=1.1)
						},
					x=x.subplot1,y=y.subplot1)
		TeachingDemos::subplot(fun = {
						# par(mar=c(0.1,0.1,0.1,0.1)) ; 
						plot(0,xlim=subplot2.x.coords,ylim=subplot2.y.coords,type='n',yaxt='n',xaxt='n',xlab="",ylab="")
						map(database="world",interior=FALSE,add=TRUE,xlim=subplot2.x.coords,ylim=subplot2.y.coords,lwd=0.5); 
						points(globe.coords,pch=20,col=continent.col,cex=0.7) ; 
							box(lwd=1.1)
						},
					x=x.subplot2,y=y.subplot2)
			box(lwd=2)
	dev.off()


	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/subsaharan_africa_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(18,35),
					ylim = c(-21,0),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			box(lwd=2)
#			legend(x = "bottomright",pch=NA,
#					legend = c("Africa","Western Eurasia","Central Eurasia","Eastern Eurasia","Oceania","Americas","population","source of admixture"),
#					text.col = c("forestgreen","blue","purple","red","brown","orange",1,1),
#					text.font = c(rep(1,6),2,3))
			# legend(x="topleft",
					# lwd = c(1,0.5,0.1),
					# col = c(adjustcolor(1,1),adjustcolor(1,0.5),adjustcolor(1,0.1)),
					# legend = c("w = 0.5","w = 0.25","w = 0.05"),
					# title = "Admixture proportions")
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/eurasia_plus_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(34.5,70),
					ylim = c(25,48),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			# legend(x="topright",
					# lwd = c(1,0.5,0.1),
					# col = c(adjustcolor(1,1),adjustcolor(1,0.5),adjustcolor(1,0.1)),
					# legend = c("w = 0.5","w = 0.25","w = 0.05"),
					# title = "Admixture proportions")
			box(lwd=2)
	dev.off()
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/eurasia_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(42.5,52.5),
					ylim = c(33.5,40),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			# legend(x="topleft",
					# lwd = c(1,0.5,0.1),
					# col = c(adjustcolor(1,1),adjustcolor(1,0.5),adjustcolor(1,0.1)),
					# legend = c("w = 0.5","w = 0.25","w = 0.05"),
					# title = "Admixture proportions")
			box(lwd=2)
	dev.off()
		

	eurasia.target.coords <- procrusteez(globe.coords[eurasia,],population.coordinates[[best]][1:k,][eurasia,],k=length(eurasia),option=1)
	eurasia.source.coords <- procrusteez(globe.coords[eurasia,],population.coordinates[[best]][1:k,][eurasia,],k=length(eurasia),population.coordinates[[best]][(k+1):(2*k),][eurasia,],option=2)
	n.africa.target.coords <- procrusteez(globe.coords[eurasia,],population.coordinates[[best]][1:k,][eurasia,],k=length(africa),population.coordinates[[best]][1:k,][africa,],option=2)
	n.africa.source.coords <- procrusteez(globe.coords[eurasia,],population.coordinates[[best]][1:k,][eurasia,],k=length(africa),population.coordinates[[best]][(k+1):(2*k),][africa,],option=2)
	x.min <- 5		#min(eurasia.target.coords[,1], eurasia.source.coords[,1]) - 5
	x.max <- 125	#max(eurasia.target.coords[,1], eurasia.source.coords[,1]) + 5
	y.min <- 14.5		#min(eurasia.target.coords[,2], eurasia.source.coords[,2]) - 5
	y.max <- 54.5		#max(eurasia.target.coords[,2], eurasia.source.coords[,2]) + 5
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/eurasia_Ad_map_indproc.png",res=300,width=7*300,height=5*300)#,pointsize=9
		par(mar=c(1,1,1,1))
		#quartz(width=7,height=5,pointsize=9)
			plot(eurasia.target.coords,type='n',
					yaxt='n',
					xaxt='n',
					xlim = c(x.min,x.max),
					ylim = c(y.min,y.max),
					xlab="",
					ylab="")
				text(eurasia.target.coords,
						labels=pops[eurasia],
						col=adjustcolor(continent.col[eurasia],1),
						font=2,cex=0.8)
				text(n.africa.target.coords,
						labels=pops[africa],
						col=adjustcolor(continent.col[africa],1),
						font=2,cex=0.8)
				text(eurasia.source.coords[,1],
						eurasia.source.coords[,2],
							labels=pops[eurasia],
							font=3,
							col=globe.admix.plot.cols[eurasia],	#continent.col[eurasia],
							cex=0.8,family="HersheySerif")
				text(n.africa.source.coords[,1],
						n.africa.source.coords[,2],
							labels=pops[africa],
							font=3,
							col=globe.admix.plot.cols[africa],	#continent.col[eurasia],
							cex=0.8,family="HersheySerif")
			arrows(	x0 = eurasia.source.coords[,1],
					y0 = eurasia.source.coords[,2],
					x1 = eurasia.target.coords[,1],
					y1 = eurasia.target.coords[,2],
					col=globe.admix.plot.cols[eurasia], #continent.col[eurasia],
					lwd=admix.proportions[eurasia,best],
					length=0.1)
			arrows(	x0 = n.africa.source.coords[,1],
					y0 = n.africa.source.coords[,2],
					x1 = n.africa.target.coords[,1],
					y1 = n.africa.target.coords[,2],
					col=globe.admix.plot.cols[africa], #continent.col[eurasia],
					lwd=admix.proportions[africa,best],
					length=0.1)
			box(lwd=2)
	dev.off()
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/western_eurasia_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(42.5,48),
					ylim = c(34.5,37.5),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			# legend(x="topleft",
					# lwd = c(1,0.5,0.1),
					# col = c(adjustcolor(1,1),adjustcolor(1,0.5),adjustcolor(1,0.1)),
					# legend = c("w = 0.5","w = 0.25","w = 0.05"),
					# title = "Admixture proportions")
			box(lwd=2)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/eastern_eurasia_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(48.5,52.2),
					ylim = c(36,39.1),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			# legend(x="bottomleft",
					# lwd = c(1,0.5,0.1),
					# col = c(adjustcolor(1,1),adjustcolor(1,0.5),adjustcolor(1,0.1)),
					# legend = c("w = 0.5","w = 0.25","w = 0.05"),
					# title = "Admixture proportions")
			box(lwd=2)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/eurasiamericas_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(42.5,52.5),
					ylim = c(33.5,46.5),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			# legend(x="topleft",
					# lwd = c(1,0.5,0.1),
					# col = c(adjustcolor(1,1),adjustcolor(1,0.5),adjustcolor(1,0.1)),
					# legend = c("w = 0.5","w = 0.25","w = 0.05"),
					# title = "Admixture proportions")
			box(lwd=2)
	dev.off()
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/eurasioceania_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(42.5,55),
					ylim = c(29,40),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			# legend(x=47.5,y=32,
					# lwd = c(1,0.5,0.1),
					# col = c(adjustcolor(1,1),adjustcolor(1,0.5),adjustcolor(1,0.1)),
					# legend = c("w = 0.5","w = 0.25","w = 0.05"),
					# title = "Admixture proportions")
			box(lwd=2)
	dev.off()
	
	pop.order <- c(africa[order(globe.coords[africa,2])],
					western.eurasia[order(globe.coords[western.eurasia,1])],
					east.asia[order(globe.coords[east.asia,1])],
					oceania[rev(order(globe.coords[oceania,2]))],
					americas[rev(order(globe.coords[americas,2]))])
	admix.cred.sets <- lapply(pop.order,FUN=function(i){quantile(admix.proportions[i,]/2,c(0.025,0.975))})
		names(admix.cred.sets) <- pops[pop.order]
	nugget.cred.sets <- lapply(pop.order,FUN=function(i){quantile(nugget[i,],c(0.025,0.975))})
		names(nugget.cred.sets) <- pops[pop.order]

	make.cred.bars <- function(quantile.vector,x.coord,bar.width,color){
		lines(x = c(x.coord-bar.width/2,x.coord+bar.width/2),
				y = c(quantile.vector[1],quantile.vector[1]),col=color)
		lines(x = c(x.coord-bar.width/2,x.coord+bar.width/2),
				y = c(quantile.vector[2],quantile.vector[2]),col=color)
		lines(x = c(x.coord,x.coord),
				y = quantile.vector,col=adjustcolor(color,0.15),lwd=0.5)
	}


	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/globe_Ad_proportions.png",res=300,width=12*300,height=5*300)
		#quartz(width=12,height=5)
		plot(rowMeans(admix.proportions)[pop.order]/2,type='n',
				main = "Mean Admixture Proportions",xlab="population",
				ylab="admixture proportion (w)",ylim=c(0,max(unlist(admix.cred.sets))))
		for(i in 1:k){
			# lines(x = c(i,i),y=c(admix.cred.sets[[i]]),col=continent.col[pop.order][i])
			make.cred.bars(admix.cred.sets[[i]],i,0.5,col=continent.col[pop.order][i])
		}
			text(rowMeans(admix.proportions)[pop.order]/2,col=continent.col[pop.order],cex=0.5,labels=pops[pop.order])
	dev.off()
								
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/globe_Ad_nugget.png",res=300,width=12*300,height=5*300)
		#quartz(width=12,height=5)
		plot(rowMeans(nugget)[pop.order],type='n',
				main = "Mean Population Nuggets",
				xlab="population",ylab="nugget",ylim=c(0,max(unlist(nugget.cred.sets))))
		for(i in 1:k){
			make.cred.bars(nugget.cred.sets[[i]],i,0.5,col=continent.col[pop.order][i])
			# lines(x = c(i,i),y=c(nugget.cred.sets[[i]]),col=continent.col[pop.order][i])
		}
			text(rowMeans(nugget)[pop.order],col=continent.col[pop.order],cex=0.5,labels=pops[pop.order])
	dev.off()
	
	globe.ad.data.table <- data.frame(cbind(round(rowMeans(admix.proportions)[pop.order]/2,4),
											round(rowMeans(nugget)[pop.order],4)),row.names=pops[pop.order])
		names(globe.ad.data.table) <- c("mean_admix_prop","mean_nugget")

	write.csv(globe.ad.data.table,file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/globe_Ad_mean_pop_adprop_nugg_vals.csv")

################
#	Admixture - real prior 3
################
load("~/Desktop/Dropbox/space.mix/data/globetrotter/globe_spacemix/globe_spaceruns/real_prior3/globe_spaceruns_realpr3_LongRun/globe_spaceruns_realpr3space_MCMC_output1.Robj")
	best <- which.max(Prob)
	target.coords <- procrusteez(globe.coords,population.coordinates[[best]][1:k,],k,option=1)
	source.coords <- procrusteez(globe.coords,population.coordinates[[best]][1:k,],k,source.locs=population.coordinates[[best]][(k+1):(2*k),],option=2) 

	globe.admix.plot.cols <- fade.admixture.source.points(continent.col,admix.proportions[,best])

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/other_globe_runs/real_prior3/globe_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(-6,80),
					ylim = c(-16,60),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
				# points(source.coords[,1],
						# source.coords[,2],
							# col=globe.admix.plot.cols,
							# pch=20)
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			box(lwd=2)
#			legend(x = "bottomright",pch=NA,
#					legend = c("Africa","Western Eurasia","Central Eurasia","Eastern Eurasia","Oceania","Americas","population","source of admixture"),
#					text.col = c("forestgreen","blue","purple","red","brown","orange",1,1),
#					text.font = c(rep(1,6),2,3),family=c())
			legend(x="topleft",
					lwd = c(1,0.5,0.1),
					col = c(adjustcolor(1,1),adjustcolor(1,0.5),adjustcolor(1,0.1)),
					legend = c("w = 0.5","w = 0.25","w = 0.05"),
					title = "Admixture proportions")
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/other_globe_runs/real_prior3/subsaharan_africa_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(-2,18),
					ylim = c(-19,6),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			box(lwd=2)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/other_globe_runs/real_prior3/eurasia_plus_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(34.5,72),
					ylim = c(23,52),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			box(lwd=2)
	dev.off()
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/other_globe_runs/real_prior3/eurasia_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(38,55),
					ylim = c(26,45),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			box(lwd=2)
	dev.off()
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/other_globe_runs/real_prior3/western_eurasia_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(44.3,48),
					ylim = c(29,35),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			box(lwd=2)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/other_globe_runs/real_prior3/eastern_eurasia_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(48,54),
					ylim = c(38,42),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			box(lwd=2)
	dev.off()
	
	pop.order <- c(africa[order(globe.coords[africa,2])],
					western.eurasia[order(globe.coords[western.eurasia,1])],
					east.asia[order(globe.coords[east.asia,1])],
					oceania[rev(order(globe.coords[oceania,2]))],
					americas[rev(order(globe.coords[americas,2]))])
	admix.cred.sets <- lapply(pop.order,FUN=function(i){quantile(admix.proportions[i,]/2,c(0.025,0.975))})
		names(admix.cred.sets) <- pops[pop.order]
	nugget.cred.sets <- lapply(pop.order,FUN=function(i){quantile(nugget[i,],c(0.025,0.975))})
		names(nugget.cred.sets) <- pops[pop.order]

	make.cred.bars <- function(quantile.vector,x.coord,bar.width,color){
		lines(x = c(x.coord-bar.width/2,x.coord+bar.width/2),
				y = c(quantile.vector[1],quantile.vector[1]),col=color)
		lines(x = c(x.coord-bar.width/2,x.coord+bar.width/2),
				y = c(quantile.vector[2],quantile.vector[2]),col=color)
		lines(x = c(x.coord,x.coord),
				y = quantile.vector,col=adjustcolor(color,0.15),lwd=0.5)
	}

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/other_globe_runs/real_prior3/globe_Ad_proportions.png",res=300,width=12*300,height=5*300)
		#quartz(width=12,height=5)
		plot(rowMeans(admix.proportions)[pop.order]/2,type='n',
				main = "Mean Admixture Proportions",xlab="population",
				ylab="admixture proportion (w)",ylim=c(0,max(unlist(admix.cred.sets))))
		for(i in 1:k){
			# lines(x = c(i,i),y=c(admix.cred.sets[[i]]),col=continent.col[pop.order][i])
			make.cred.bars(admix.cred.sets[[i]],i,0.5,col=continent.col[pop.order][i])
		}
			text(rowMeans(admix.proportions)[pop.order]/2,col=continent.col[pop.order],cex=0.5,labels=pops[pop.order])
	dev.off()
								
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/other_globe_runs/real_prior3/globe_Ad_nugget.png",res=300,width=12*300,height=5*300)
		#quartz(width=12,height=5)
		plot(rowMeans(nugget)[pop.order],type='n',
				main = "Mean Population Nuggets",
				xlab="population",ylab="nugget",ylim=c(0,max(unlist(nugget.cred.sets))))
		for(i in 1:k){
			make.cred.bars(nugget.cred.sets[[i]],i,0.5,col=continent.col[pop.order][i])
			# lines(x = c(i,i),y=c(nugget.cred.sets[[i]]),col=continent.col[pop.order][i])
		}
			text(rowMeans(nugget)[pop.order],col=continent.col[pop.order],cex=0.5,labels=pops[pop.order])
	dev.off()
	
	globe.ad.data.table <- data.frame(cbind(round(rowMeans(admix.proportions)[pop.order]/2,4),
											round(rowMeans(nugget)[pop.order],4)),row.names=pops[pop.order])
		names(globe.ad.data.table) <- c("mean_admix_prop","mean_nugget")

	write.csv(globe.ad.data.table,file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/globe_Ad_mean_pop_adprop_nugg_vals.csv")

################
#	Admixture - real prior 2
################
load("~/Desktop/Dropbox/space.mix/data/globetrotter/globe_spacemix/globe_spaceruns/real_prior2/globe_spaceruns_realpr2_LongRun/globe_spaceruns_realpr2space_MCMC_output1.Robj")

	best <- which.max(Prob)
	target.coords <- procrusteez(globe.coords,population.coordinates[[best]][1:k,],k,option=1)
	source.coords <- procrusteez(globe.coords,population.coordinates[[best]][1:k,],k,source.locs=population.coordinates[[best]][(k+1):(2*k),],option=2) 

	globe.admix.plot.cols <- fade.admixture.source.points(continent.col,admix.proportions[,best])

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/other_globe_runs/real_prior2/globe_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(0,70),
					ylim = c(-16,48),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			box(lwd=2)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/other_globe_runs/real_prior2/subsaharan_africa_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(7,17),
					ylim = c(-15,0),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			box(lwd=2)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/other_globe_runs/real_prior2/eurasia_plus_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(41,60),
					ylim = c(29,46),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			box(lwd=2)
	dev.off()
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/other_globe_runs/real_prior2/eurasia_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(46,51),
					ylim = c(30,39),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			box(lwd=2)
	dev.off()
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/other_globe_runs/real_prior2/western_eurasia_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(46.2,49),
					ylim = c(32,35),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			box(lwd=2)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/other_globe_runs/real_prior2/eastern_eurasia_Ad_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(47.5,50.7),
					ylim = c(36.8,38.8),
					xlab="Eastings",
					ylab="Northings")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,1),
						font=2,cex=0.8)
				text(source.coords[,1],
						source.coords[,2],
							labels=pops,
							font=3,
							col=globe.admix.plot.cols,cex=0.8,family="HersheySerif")
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=admix.proportions[,best],
					length=0.1)
			box(lwd=2)
	dev.off()
	
	pop.order <- c(africa[order(globe.coords[africa,2])],
					western.eurasia[order(globe.coords[western.eurasia,1])],
					east.asia[order(globe.coords[east.asia,1])],
					oceania[rev(order(globe.coords[oceania,2]))],
					americas[rev(order(globe.coords[americas,2]))])
	admix.cred.sets <- lapply(pop.order,FUN=function(i){quantile(admix.proportions[i,]/2,c(0.025,0.975))})
		names(admix.cred.sets) <- pops[pop.order]
	nugget.cred.sets <- lapply(pop.order,FUN=function(i){quantile(nugget[i,],c(0.025,0.975))})
		names(nugget.cred.sets) <- pops[pop.order]

	make.cred.bars <- function(quantile.vector,x.coord,bar.width,color){
		lines(x = c(x.coord-bar.width/2,x.coord+bar.width/2),
				y = c(quantile.vector[1],quantile.vector[1]),col=color)
		lines(x = c(x.coord-bar.width/2,x.coord+bar.width/2),
				y = c(quantile.vector[2],quantile.vector[2]),col=color)
		lines(x = c(x.coord,x.coord),
				y = quantile.vector,col=adjustcolor(color,0.15),lwd=0.5)
	}

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/other_globe_runs/real_prior2/globe_Ad_proportions.png",res=300,width=12*300,height=5*300)
		#quartz(width=12,height=5)
		plot(rowMeans(admix.proportions)[pop.order]/2,type='n',
				main = "Mean Admixture Proportions",xlab="population",
				ylab="admixture proportion (w)",ylim=c(0,max(unlist(admix.cred.sets))))
		for(i in 1:k){
			# lines(x = c(i,i),y=c(admix.cred.sets[[i]]),col=continent.col[pop.order][i])
			make.cred.bars(admix.cred.sets[[i]],i,0.5,col=continent.col[pop.order][i])
		}
			text(rowMeans(admix.proportions)[pop.order]/2,col=continent.col[pop.order],cex=0.5,labels=pops[pop.order])
	dev.off()
								
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globetrotter/other_globe_runs/real_prior2/globe_Ad_nugget.png",res=300,width=12*300,height=5*300)
		#quartz(width=12,height=5)
		plot(rowMeans(nugget)[pop.order],type='n',
				main = "Mean Population Nuggets",
				xlab="population",ylab="nugget",ylim=c(0,max(unlist(nugget.cred.sets))))
		for(i in 1:k){
			make.cred.bars(nugget.cred.sets[[i]],i,0.5,col=continent.col[pop.order][i])
			# lines(x = c(i,i),y=c(nugget.cred.sets[[i]]),col=continent.col[pop.order][i])
		}
			text(rowMeans(nugget)[pop.order],col=continent.col[pop.order],cex=0.5,labels=pops[pop.order])
	dev.off()


################################
#	SIM FIGS
################################

################
#	Simulation scenarios
################
source("~/Desktop/Dropbox/space.mix/sims/spacemix_ms_sims.R")
png(file="~/Desktop/Dropbox/space.mix/ms/figs/sims/basic_lattice.png",res=200,width=6*200,height=5*200)
	migration.rate.graphic(x.pops = 5,y.pops = 6,migration.rate=1,jitter=0.25,labels=TRUE,colors=TRUE)
dev.off()

png(file="~/Desktop/Dropbox/space.mix/ms/figs/sims/barrier_lattice.png",res=200,width=6*200,height=5*200)
	migration.rate.graphic(x.pops = 5,y.pops = 6,migration.rate=1,jitter=0.25,barrier.effect=5,labels=TRUE,colors=TRUE)
dev.off()

parents <- c(78:88)
time.points <- rep(0.07,11)
expansion.list <- vector(mode="list")
	for(i in 1:length(parents)){
		expansion.list[[i]] <- list(parent=parents[i],
									daughters = parents[i]+c(11,22,33,44,55),
									time.point = time.points[i])
	}
	
png(file="~/Desktop/Dropbox/space.mix/ms/figs/sims/expansion_lattice.png",res=200,width=6*200,height=5*200)
	migration.rate.graphic(x.pop=5,y.pops=6,migration.rate=1,jitter = 0.25,expansion.list=expansion.list,labels=TRUE,colors=TRUE,ylim=c(0,10.2),curve=0.3)
dev.off()


expansion.list2 <- list(list(parent = 61,
						daughters = 105,
						time.point = 1))
	
png(file="~/Desktop/Dropbox/space.mix/ms/figs/sims/barr_indland_ad_lattice.png",res=200,width=6*200,height=5*200)
	migration.rate.graphic(x.pop=5,y.pops=6,migration.rate=1,barrier.effect=5,jitter = 0.25,expansion.list=expansion.list2,labels=TRUE,colors=TRUE,arrow.col="green")
dev.off()

expansion.list2.5 <- list(list(parent = 61,
						daughters = 83,
						time.point = 1))
	
png(file="~/Desktop/Dropbox/space.mix/ms/figs/sims/big_barr_ad_lattice.png",res=200,width=6*200,height=5*200)
	migration.rate.graphic(x.pop=5,y.pops=6,migration.rate=1,barrier.effect=5,jitter = 0.25,expansion.list=expansion.list2.5,labels=TRUE,colors=TRUE,arrow.col="green")
dev.off()

expansion.list3 <- list(list(parent = 13,
						daughters = 131,
						time.point = 1))
png(file="~/Desktop/Dropbox/space.mix/ms/figs/sims/corner_admixture_lattice.png",res=200,width=6*200,height=5*200)
	migration.rate.graphic(x.pop=5,y.pops=6,migration.rate=1,jitter = 0.25,expansion.list=expansion.list3,labels=TRUE,colors=TRUE,arrow.col="green",curve=0.2)
dev.off()

################
#	Grid
################
load("~/Desktop/Dropbox/space.mix/sims/stationary_pops/spacemix/stationary_pops_1/spacemix_ms_sim_stationary_pops_1space_MCMC_output1.Robj")
load("~/Desktop/Dropbox/space.mix/sims/stationary_pops/spacemix/stationary_pops_1/spacemix.ms.dataset_stationary_pops.Robj")
k <- last.params$k
best <- which.max(Prob)
pop.cols <- rainbow(k,start=4/6,end=6/6)[as.numeric(cut(spacemix.dataset$population.coordinates[,1],k))]
target.coords <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
							target.locs = population.coordinates[[best]][1:k,],
							k = k,
							option = 1)
pdf("~/Desktop/Dropbox/space.mix/ms/figs/sims/GeoGenMap_lattice.pdf",width=6,height=5,pointsize=9)
#quartz(width=6,height=5)
par(mar=c(4.3,4.3,3,1))
plot(target.coords,pch=1,xlim=c(0,12),ylim=c(-0.4,10),cex=3.5,
		xlab="Eastings",ylab="Northings",
		main="Inferred Population Map:\n Simple Lattice",
		col=pop.cols,lwd=2)
	box(lwd=2)
	text(target.coords,labels=paste(1:k))
dev.off()

################
#	Barrier
################
load("~/Desktop/Dropbox/space.mix/sims/barrier/spacemix/barrier_1/spacemix_ms_sim_barrier_1space_MCMC_output1.Robj")
load("~/Desktop/Dropbox/space.mix/sims/barrier/spacemix.ms.dataset.Robj")
k <- last.params$k
best <- which.max(Prob)
target.coords <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
							target.locs = population.coordinates[[best]][1:k,],
							k = k,
							option = 1)
pdf("~/Desktop/Dropbox/space.mix/ms/figs/sims/GeoGenMap_barrier.pdf",width=6,height=5,pointsize=9)
#quartz(width=6,height=5)
par(mar=c(4.3,4.3,3,1))
plot(target.coords,pch=1,xlim=c(0,12),ylim=c(-0.4,10),cex=3.5,
		xlab="Eastings",ylab="Northings",
		main="Inferred Population Map:\n Lattice with Barrier",
		col=pop.cols,lwd=2)
	box(lwd=2)
	text(target.coords,labels=paste(1:k))
dev.off()


################
#	Expansion
################

load("~/Desktop/Dropbox/space.mix/sims/expansion/spacemix/noad/rand_prior1/expansion_randpr_noad_1_LongRun/expansion_randpr_noad_1space_MCMC_output1.Robj")
load("~/Desktop/Dropbox/space.mix/sims/expansion/sim_expansion_dataset.Robj")
k <- last.params$k
best <- which.max(Prob)
target.coords <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
							target.locs = population.coordinates[[best]][1:k,],
							k = k,
							option = 1)
pdf("~/Desktop/Dropbox/space.mix/ms/figs/sims/GeoGenMap_expansion.pdf",width=6,height=5,pointsize=9)
#quartz(width=6,height=5)
par(mar=c(4.3,4.3,3,1))
plot(target.coords,pch=1,xlim=c(1.2,9.5),ylim=c(-2,12),cex=3.5,
		xlab="Eastings",ylab="Northings",
		main="Inferred Population Map:\n Lattice with Expansion",
		col=pop.cols,lwd=2)
	box(lwd=2)
	text(target.coords,labels=paste(1:k))
dev.off()


################
#	Corner Admixture
################

load("~/Desktop/Dropbox/space.mix/sims/admixture/spacemix/admixture_4/spacemix_ms_sim_admixture_4space_MCMC_output1.Robj")
load("~/Desktop/Dropbox/space.mix/sims/admixture/spacemix/admixture_4/spacemix.ms.dataset.Robj")
k <- last.params$k
best <- which.max(Prob)
target.coords <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
							target.locs = population.coordinates[[best]][1:k,],
							k = k,
							option = 1)
source.coords <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
								target.locs = population.coordinates[[best]][1:k,],
								source.locs = population.coordinates[[best]][(k+1):(2*k),],
								k = k,
								option = 2)
x.min <- min(target.coords[,1]) - 0.5
x.max <- max(target.coords[,1]) + 0.5
y.min <- min(target.coords[,2]) - 0.5
y.max <- max(target.coords[,2]) + 1
source.coord.cols <- fade.admixture.source.points(pop.cols,admix.proportions[,best])
png("~/Desktop/Dropbox/space.mix/ms/figs/sims/GeoGenMap_corner_admixture.png",res=300,width=6*300,height=5*300,pointsize=9)
#quartz(width=6,height=5)
par(mar=c(4.3,4.3,3,1))
plot(target.coords,xlim=c(x.min,x.max),ylim=c(y.min,y.max),pch=1,cex=3.5,
		xlab="Eastings",ylab="Northings",main="Inferred Population Map:\n Lattice with Admixture",
		col=pop.cols,lwd=2)
		points(source.coords,pch=20,col=source.coord.cols)
	box(lwd=2)
	text(target.coords,labels=paste(1:k))
	arrows(x0 = source.coords[,1],
			y0 = source.coords[,2],
			x1 = target.coords[,1],
			y1 = target.coords[,2],
			col= source.coord.cols,
			lwd = admix.proportions[,best],
			length=0.2)
	legend(x = "topleft",
			pch=c(20,rep(NA,5)),
			lty=c(NA,rep(1,5)),
			lwd=c(NA,0.2,0.4,0.6,0.8,1),
			col=c(1,adjustcolor(1,0.2),adjustcolor(1,0.4),adjustcolor(1,0.6),adjustcolor(1,0.8),adjustcolor(1,1)),
			legend = c("admixture source","w = 0.1","w = 0.2","w = 0.3","w = 0.4","w = 0.5"),
			cex=0.9)
dev.off()

################
#	Corner Admixture - CYOL
################

load("~/Desktop/Dropbox/space.mix/sims/admixture/spacemix/admixture_1/spacemix_ms_sim_admixture_1space_MCMC_output1.Robj")
load("~/Desktop/Dropbox/space.mix/sims/admixture/spacemix/admixture_1/spacemix.ms.dataset.Robj")
k <- last.params$k
best <- which.max(Prob)
target.coords <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
							target.locs = population.coordinates[[best]][1:k,],
							k = k,
							option = 1)
x.min <- min(target.coords[,1]) - 0.5
x.max <- max(target.coords[,1]) + 0.5
y.min <- min(target.coords[,2]) - 0.5
y.max <- max(target.coords[,2]) + 1

png("~/Desktop/Dropbox/space.mix/ms/figs/sims/GeoGenMap_corner_admixture_CYOL.png",res=300,width=6*300,height=5*300,pointsize=9)
#quartz(width=6,height=5)
par(mar=c(4.3,4.3,3,1))
plot(target.coords,xlim=c(x.min,x.max),ylim=c(y.min,y.max),pch=1,cex=3.5,
		xlab="Eastings",ylab="Northings",main="Inferred Population Map:\n Lattice with Admixture",
		col=pop.cols,lwd=2)
	box(lwd=2)
	text(target.coords,labels=paste(1:k))
dev.off()

################
#	Corner Admixture - just admixture
################

load("~/Desktop/Dropbox/space.mix/sims/admixture/spacemix/admixture_2/spacemix_ms_admixture_2space_MCMC_output1.Robj")
load("~/Desktop/Dropbox/space.mix/sims/admixture/spacemix/admixture_2/spacemix.ms.dataset.Robj")
k <- last.params$k
best <- which.max(Prob)
target.coords <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
							target.locs = population.coordinates[[best]][1:k,],
							k = k,
							option = 1)
source.coords <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
								target.locs = population.coordinates[[best]][1:k,],
								source.locs = population.coordinates[[best]][(k+1):(2*k),],
								k = k,
								option = 2)

x.min <- min(target.coords[,1]) - 0.5
x.max <- max(target.coords[,1]) + 0.5
y.min <- min(target.coords[,2]) - 0.5
y.max <- max(target.coords[,2]) + 1
source.coord.cols <- fade.admixture.source.points(pop.cols,admix.proportions[,best])
png("~/Desktop/Dropbox/space.mix/ms/figs/sims/GeoGenMap_corner_admixture_adinf.png",res=300,width=6*300,height=5*300,pointsize=9)
#quartz(width=6,height=5)
par(mar=c(4.3,4.3,3,1))
plot(target.coords,xlim=c(x.min,x.max),ylim=c(y.min,y.max),pch=1,cex=3.5,
		xlab="Eastings",ylab="Northings",main="Inferred Population Map:\n Lattice with Admixture",
		col=pop.cols,lwd=2)
		points(source.coords,pch=20,col=source.coord.cols)
	box(lwd=2)
	text(target.coords,labels=paste(1:k))
	arrows(x0 = source.coords[,1],
			y0 = source.coords[,2],
			x1 = target.coords[,1],
			y1 = target.coords[,2],
			col= source.coord.cols,
			lwd = admix.proportions[,best],
			length=0.2)
dev.off()

################
#	Barrier w/ Admixture - Inland
################

load("~/Desktop/Dropbox/space.mix/sims/bar_inland_ad/bar_inland_ad_spacemix/bar_inland_ad_spacemix1/rand_prior/bar_inland_ad_randpr__LongRun/bar_inland_ad_randpr_space_MCMC_output1.Robj")
load("~/Desktop/Dropbox/space.mix/sims/bar_inland_ad/bar_inland_ad_spacemix/bar_inland_ad_spacemix1/rand_prior/barr_inland_ad_dataset.Robj")
k <- last.params$k
best <- which.max(Prob)
target.coords <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
							target.locs = population.coordinates[[best]][1:k,],
							k = k,
							option = 1)
source.coords <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
								target.locs = population.coordinates[[best]][1:k,],
								source.locs = population.coordinates[[best]][(k+1):(2*k),],
								k = k,
								option = 2)
target.coords.list <- vector(mode="list",length = length(which(Prob!=0)))
source.coords.list <- vector(mode="list",length = length(which(Prob!=0)))
	for(i in 1:length(source.coords.list)){
		source.coords.list[[i]] <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
								target.locs = population.coordinates[[i]][1:k,],
								source.locs = population.coordinates[[i]][(k+1):(2*k),],
								k = k,
								option = 2)
		target.coords.list[[i]] <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
							target.locs = population.coordinates[[i]][1:k,],
							k = k,
							option = 1)
	}
x.min <- min(target.coords[,1]) - 0.5
x.max <- max(target.coords[,1]) + 0.5
y.min <- min(target.coords[,2]) - 0.5
y.max <- max(target.coords[,2]) + 2
scalar <- 4
source.coord.cols <- fade.admixture.source.points(pop.cols,scalar*admix.proportions[,best])
png("~/Desktop/Dropbox/space.mix/ms/figs/sims/GeoGenMap_barr_inland_admixture_1.png",res=300,width=6*300,height=5*300,pointsize=9)
#quartz(width=6,height=5)
par(mar=c(4.3,4.3,3,1))
plot(target.coords,xlim=c(x.min,x.max),ylim=c(y.min,y.max),pch=1,cex=3.5,
		xlab="Eastings",ylab="Northings",main="Inferred Population Map:\n Lattice with Barrier and Admixture",
		col=pop.cols,lwd=2)
		points(source.coords,pch=20,col=source.coord.cols)
	box(lwd=2)
	text(target.coords,labels=paste(1:k))
	arrows(x0 = source.coords[,1],
			y0 = source.coords[,2],
			x1 = target.coords[,1],
			y1 = target.coords[,2],
			col= source.coord.cols,
			lwd = scalar*admix.proportions[,best],
			length=0.2)
	legend(x = "topleft",
			pch=c(20,rep(NA,5)),
			lty=c(NA,rep(1,5)),
			lwd=c(NA,0.2,0.4,0.6,0.8,1),
			col=c(1,adjustcolor(1,0.2),adjustcolor(1,0.4),adjustcolor(1,0.6),adjustcolor(1,0.8),adjustcolor(1,1)),
			legend = c("admixture source",
								paste("w = ",round(0.1/scalar,2),sep=""),
								paste("w = ",round(0.2/scalar,2),sep=""),
								paste("w = ",round(0.3/scalar,2),sep=""),
								paste("w = ",round(0.4/scalar,2),sep=""),
								paste("w = ",round(0.5/scalar,2),sep="")),
			cex=0.9)
dev.off()


burnin <- 0.5
x <- seq(length(source.coords.list)*burnin + 1,length(source.coords.list),1)
ad.pop.source.coords.mat <- matrix(0,nrow=length(x),ncol=2)
	for(i in 1:nrow(ad.pop.source.coords.mat)){
		ad.pop.source.coords.mat[i,] <- source.coords.list[[x[i]]][23,]
	}

	#ad.pop.source.coords.cred.set <- ellipse(mu = colMeans(ad.pop.source.coords.mat), sigma = cov(ad.pop.source.coords.mat), alpha = 0.4,draw=FALSE)
require(hdrcde)
source.coord.hdr <- hdr.2d(ad.pop.source.coords.mat[,1],ad.pop.source.coords.mat[,2],prob=0.3)
source.coord.hdr2 <- source.coord.hdr
source.coord.hdr2$mode <- c(100,100)

png("~/Desktop/Dropbox/space.mix/ms/figs/sims/GeoGenMap_barr_inland_admixture_2.png",res=300,width=6*300,height=5*300,pointsize=9)
#quartz(width=6,height=5)
par(mar=c(4.3,4.3,3,1))
plot(source.coord.hdr2,shadecols=adjustcolor(pop.cols[k],0.2),show.points=FALSE,outside.points=FALSE,
		xlab="Eastings",ylab="Northings",main="Inferred Population Map:\n Lattice with Barrier and Admixture",
		xlim=c(x.min,x.max),
		ylim=c(y.min,y.max))
points(target.coords,xlim=c(x.min,x.max),ylim=c(y.min,y.max),pch=1,cex=3.5,col=pop.cols,lwd=2)
	# points(ad.pop.source.coords.mat,col=adjustcolor("blue",0.2))
	box(lwd=2)
		points(source.coords,pch=20,col=source.coord.cols)
	text(target.coords,labels=paste(1:k))
		arrows(x0 = source.coords[,1],
			y0 = source.coords[,2],
			x1 = target.coords[,1],
			y1 = target.coords[,2],
			col= source.coord.cols,
			lwd = scalar*admix.proportions[,best],
			length=0.2)
	# for(i in 1:nrow(ad.pop.source.coords.mat)){
		# arrows(x0 = ad.pop.source.coords.mat[i,1],
				# y0 = ad.pop.source.coords.mat[i,2],
				# x1 = target.coords[23,1],
				# y1 = target.coords[23,2],
				# col= adjustcolor(1,admix.proportions[23,x[i]]),
				# lwd = scalar*admix.proportions[23,x[i]],
				# length=0.2)
	# }
	legend(x = "topleft",
			pch=c(20,rep(NA,5)),
			lty=c(NA,rep(1,5)),
			lwd=c(NA,0.2,0.4,0.6,0.8,1),
			col=c(1,adjustcolor(1,0.2),adjustcolor(1,0.4),adjustcolor(1,0.6),adjustcolor(1,0.8),adjustcolor(1,1)),
			legend = c("admixture source",
								paste("w = ",round(0.1/scalar,2),sep=""),
								paste("w = ",round(0.2/scalar,2),sep=""),
								paste("w = ",round(0.3/scalar,2),sep=""),
								paste("w = ",round(0.4/scalar,2),sep=""),
								paste("w = ",round(0.5/scalar,2),sep="")),
			cex=0.9)
dev.off()


################
#	Barrier w/ Admixture - Nearest Neighbor
################

load("~/Desktop/Dropbox/space.mix/sims/big_barr_ad/big_barr_ad_spacemix/rand_prior1/big_barr_ad_randpr_1_LongRun/big_barr_ad_randpr_1space_MCMC_output1.Robj")
load("~/Desktop/Dropbox/space.mix/sims/big_barr_ad/big_barr_ad_spacemix/rand_prior1/sim_big_barr_ad_dataset.Robj")
k <- last.params$k
best <- which.max(Prob)
target.coords <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
							target.locs = population.coordinates[[best]][1:k,],
							k = k,
							option = 1)
source.coords <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
								target.locs = population.coordinates[[best]][1:k,],
								source.locs = population.coordinates[[best]][(k+1):(2*k),],
								k = k,
								option = 2)
target.coords.list <- vector(mode="list",length = length(which(Prob!=0)))
source.coords.list <- vector(mode="list",length = length(which(Prob!=0)))
	for(i in 1:length(source.coords.list)){
		source.coords.list[[i]] <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
								target.locs = population.coordinates[[i]][1:k,],
								source.locs = population.coordinates[[i]][(k+1):(2*k),],
								k = k,
								option = 2)
		target.coords.list[[i]] <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
							target.locs = population.coordinates[[i]][1:k,],
							k = k,
							option = 1)
	}
x.min <- min(target.coords[,1]) - 0.5
x.max <- max(target.coords[,1]) + 0.5
y.min <- min(target.coords[,2]) - 0.5
y.max <- max(target.coords[,2]) + 2
scalar <- 4
source.coord.cols <- fade.admixture.source.points(pop.cols,scalar*admix.proportions[,best])
png("~/Desktop/Dropbox/space.mix/ms/figs/sims/GeoGenMap_big_barr_ad_1.png",res=300,width=6*300,height=5*300,pointsize=9)
#quartz(width=6,height=5)
par(mar=c(4.3,4.3,3,1))
plot(target.coords,xlim=c(x.min,x.max),ylim=c(y.min,y.max),pch=1,cex=3.5,
		xlab="Eastings",ylab="Northings",main="Inferred Population Map:\n Lattice with Barrier and Admixture",
		col=pop.cols,lwd=2)
		points(source.coords,pch=20,col=source.coord.cols)
	box(lwd=2)
	text(target.coords,labels=paste(1:k))
	arrows(x0 = source.coords[,1],
			y0 = source.coords[,2],
			x1 = target.coords[,1],
			y1 = target.coords[,2],
			col= source.coord.cols,
			lwd = scalar*admix.proportions[,best],
			length=0.2)
	legend(x = "topleft",
			pch=c(20,rep(NA,5)),
			lty=c(NA,rep(1,5)),
			lwd=c(NA,0.2,0.4,0.6,0.8,1),
			col=c(1,adjustcolor(1,0.2),adjustcolor(1,0.4),adjustcolor(1,0.6),adjustcolor(1,0.8),adjustcolor(1,1)),
			legend = c("admixture source",
								paste("w = ",round(0.1/scalar,2),sep=""),
								paste("w = ",round(0.2/scalar,2),sep=""),
								paste("w = ",round(0.3/scalar,2),sep=""),
								paste("w = ",round(0.4/scalar,2),sep=""),
								paste("w = ",round(0.5/scalar,2),sep="")),
			cex=0.9)
dev.off()


burnin <- 0.5
x <- seq(length(source.coords.list)*burnin + 1,length(source.coords.list),1)
ad.pop.source.coords.mat <- matrix(0,nrow=length(x),ncol=2)
	for(i in 1:nrow(ad.pop.source.coords.mat)){
		ad.pop.source.coords.mat[i,] <- source.coords.list[[x[i]]][18,]
	}

	#ad.pop.source.coords.cred.set <- ellipse(mu = colMeans(ad.pop.source.coords.mat), sigma = cov(ad.pop.source.coords.mat), alpha = 0.4,draw=FALSE)
require(hdrcde)
source.coord.hdr <- hdr.2d(ad.pop.source.coords.mat[,1],ad.pop.source.coords.mat[,2],prob=0.05)
source.coord.hdr2 <- source.coord.hdr
source.coord.hdr2$mode <- c(100,100)

png("~/Desktop/Dropbox/space.mix/ms/figs/sims/GeoGenMap_barr_inland_admixture_2.png",res=300,width=6*300,height=5*300,pointsize=9)
#quartz(width=6,height=5)
par(mar=c(4.3,4.3,3,1))
plot(source.coord.hdr2,shadecols=adjustcolor(pop.cols[k],0.2),show.points=FALSE,outside.points=FALSE,
		xlab="Eastings",ylab="Northings",main="Inferred Population Map:\n Lattice with Barrier and Admixture",
		xlim=c(x.min,x.max),
		ylim=c(y.min,y.max))
points(target.coords,xlim=c(x.min,x.max),ylim=c(y.min,y.max),pch=1,cex=3.5,col=pop.cols,lwd=2)
	# points(ad.pop.source.coords.mat,col=adjustcolor("blue",0.2))
	box(lwd=2)
		points(source.coords,pch=20,col=source.coord.cols)
	text(target.coords,labels=paste(1:k))
		arrows(x0 = source.coords[,1],
			y0 = source.coords[,2],
			x1 = target.coords[,1],
			y1 = target.coords[,2],
			col= source.coord.cols,
			lwd = scalar*admix.proportions[,best],
			length=0.2)
	# for(i in 1:nrow(ad.pop.source.coords.mat)){
		# arrows(x0 = ad.pop.source.coords.mat[i,1],
				# y0 = ad.pop.source.coords.mat[i,2],
				# x1 = target.coords[23,1],
				# y1 = target.coords[23,2],
				# col= adjustcolor(1,admix.proportions[23,x[i]]),
				# lwd = scalar*admix.proportions[23,x[i]],
				# length=0.2)
	# }
	legend(x = "topleft",
			pch=c(20,rep(NA,5)),
			lty=c(NA,rep(1,5)),
			lwd=c(NA,0.2,0.4,0.6,0.8,1),
			col=c(1,adjustcolor(1,0.2),adjustcolor(1,0.4),adjustcolor(1,0.6),adjustcolor(1,0.8),adjustcolor(1,1)),
			legend = c("admixture source",
								paste("w = ",round(0.1/scalar,2),sep=""),
								paste("w = ",round(0.2/scalar,2),sep=""),
								paste("w = ",round(0.3/scalar,2),sep=""),
								paste("w = ",round(0.4/scalar,2),sep=""),
								paste("w = ",round(0.5/scalar,2),sep="")),
			cex=0.9)
dev.off()



################
#	Acceptance rate figs
################
load("~/Desktop/Dropbox/space.mix/sims/big_barr_ad/big_barr_ad_spacemix/rand_prior1/big_barr_ad_randpr_1_LongRun/big_barr_ad_randpr_1space_MCMC_output1.Robj")
load("~/Desktop/Dropbox/space.mix/sims/big_barr_ad/sim_big_barr_ad_dataset.Robj")
k <- last.params$k
pop.cols <- rainbow(k,start=4/6,end=6/6)[as.numeric(cut(spacemix.dataset$population.coordinates[,1],k))]

png(file="~/Desktop/Dropbox/space.mix/ms/figs/sims/example_acceptance_rates.png",res=300,width=12*300,height=5*300)
#quartz(width=12,height=5)
par(mfrow=c(1,2))
	plot(accept_rates$a2_accept_rate,
			type='l',ylab=expression(paste(alpha[2]," proportion accepted moves",sep="")),
			xlab="MCMC iterations",
			ylim=c(0.3,0.55))
	mtext(side=3,text="Adaptive Metropolis-within-Gibbs Proposal Mechanism",font=2,padj=-2.5,adj=8)
	matplot(t(accept_rates$nugget_accept_rate),col=pop.cols,type='l',
				ylab="nugget proportion accepted moves",
				xlab="MCMC iterations",
				ylim=c(0.3,0.55))
dev.off()





