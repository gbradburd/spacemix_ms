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

png(file="~/Desktop/Dropbox/space.mix/ms/figs/warb_pop_PC_map.png",res=200,width=5*200,height=4*200)
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/population_warbler_map_realpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/population_warbler_map_no_arrows_realpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/population_warbler_admix_values_nugget_realpr1.png",res=300,width=4.5*300,height=5*300,pointsize=9)
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/population_warbler_map_realpr2.png",res=300,width=7*300,height=5*300,pointsize=9)
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/population_warbler_map_no_arrows_realpr2.png",res=300,width=7*300,height=5*300,pointsize=9)
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/population_warbler_admix_values_nugget_realpr2.png",res=300,width=4.5*300,height=5*300,pointsize=9)
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/population_warbler_map_randpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/population_warbler_map_no_arrows_randpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/population_warbler_admix_values_nugget_randpr1.png",res=300,width=4.5*300,height=5*300,pointsize=9)
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

png(file="~/Desktop/Dropbox/space.mix/ms/figs/warb_ind_PC_map.png",res=200,width=5*200,height=4*200)
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
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/individual_warbler_map_arrows_amped_admixture_realpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1],source.coords[,1]),
							max(target.coords[,1],source.coords[,1])),
					ylim=c(min(target.coords[,2],source.coords[,2]),
							max(target.coords[,2],source.coords[,2])),
					xlab="long",
					ylab="lat")
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/individual_warbler_map_arrows_realpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1],source.coords[,1]),
							max(target.coords[,1],source.coords[,1])),
					ylim=c(min(target.coords[,2],source.coords[,2]),
							max(target.coords[,2],source.coords[,2])),
					xlab="long",
					ylab="lat")
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/individual_warbler_map_noarrows_realpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1]),max(target.coords[,1])+2),
					ylim=c(min(target.coords[,2]),max(target.coords[,2])),
					xlab="long",
					ylab="lat")
				text(target.coords[c(1:k),],
						labels=plot.inds,
						col=adjustcolor(inds.col,0.8),
						font=2,cex=0.9)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/individual_warbler_map_noarrows_closeup_realpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[-c(which(plot.inds=="TU")),1]),max(target.coords[-c(which(plot.inds=="TU")),1])),
					ylim=c(min(target.coords[-c(which(plot.inds=="TU")),2]),max(target.coords[-c(which(plot.inds=="TU")),2])),
					xlab="long",
					ylab="lat")
				text(target.coords[-c(which(plot.inds=="TU")),],
						labels=plot.inds[-c(which(plot.inds=="TU"))],
						col=adjustcolor(inds.col[-c(which(plot.inds=="TU"))],0.8),
						font=2,cex=0.9)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/individual_warbler_map_noarrows_closeup_nugget_realpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[-c(which(plot.inds=="TU")),1]),max(target.coords[-c(which(plot.inds=="TU")),1])),
					ylim=c(min(target.coords[-c(which(plot.inds=="TU")),2]),max(target.coords[-c(which(plot.inds=="TU")),2])),
					xlab="long",
					ylab="lat")
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
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/individual_warbler_map_arrows_amped_admixture_realpr2.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1],source.coords[,1]),
							max(target.coords[,1],source.coords[,1])),
					ylim=c(min(target.coords[,2],source.coords[,2]),
							max(target.coords[,2],source.coords[,2])),
					xlab="long",
					ylab="lat")
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/individual_warbler_map_arrows_realpr2.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1],source.coords[,1]),
							max(target.coords[,1],source.coords[,1])),
					ylim=c(min(target.coords[,2],source.coords[,2]),
							max(target.coords[,2],source.coords[,2])),
					xlab="long",
					ylab="lat")
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/individual_warbler_map_noarrows_realpr2.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1]),max(target.coords[,1])+2),
					ylim=c(min(target.coords[,2]),max(target.coords[,2])),
					xlab="long",
					ylab="lat")
				text(target.coords[c(1:k),],
						labels=plot.inds,
						col=adjustcolor(inds.col,0.8),
						font=2,cex=0.9)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/individual_warbler_map_noarrows_closeup_realpr2.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[-c(which(plot.inds=="TU")),1]),max(target.coords[-c(which(plot.inds=="TU")),1])),
					ylim=c(min(target.coords[-c(which(plot.inds=="TU")),2]),max(target.coords[-c(which(plot.inds=="TU")),2])),
					xlab="long",
					ylab="lat")
				text(target.coords[-c(which(plot.inds=="TU")),],
						labels=plot.inds[-c(which(plot.inds=="TU"))],
						col=adjustcolor(inds.col[-c(which(plot.inds=="TU"))],0.8),
						font=2,cex=0.9)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/individual_warbler_map_noarrows_closeup_nugget_realpr2.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[-c(which(plot.inds=="TU")),1]),max(target.coords[-c(which(plot.inds=="TU")),1])),
					ylim=c(min(target.coords[-c(which(plot.inds=="TU")),2]),max(target.coords[-c(which(plot.inds=="TU")),2])),
					xlab="long",
					ylab="lat")
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
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/individual_warbler_map_arrows_amped_admixture_randpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1],source.coords[,1]),
							max(target.coords[,1],source.coords[,1])),
					ylim=c(min(target.coords[,2],source.coords[,2]),
							max(target.coords[,2],source.coords[,2])),
					xlab="long",
					ylab="lat")
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/individual_warbler_map_arrows_randpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1],source.coords[,1]),
							max(target.coords[,1],source.coords[,1])),
					ylim=c(min(target.coords[,2],source.coords[,2]),
							max(target.coords[,2],source.coords[,2])),
					xlab="long",
					ylab="lat")
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/individual_warbler_map_noarrows_randpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[,1]),max(target.coords[,1])+2),
					ylim=c(min(target.coords[,2]),max(target.coords[,2])),
					xlab="long",
					ylab="lat")
				text(target.coords[c(1:k),],
						labels=plot.inds,
						col=adjustcolor(inds.col,0.8),
						font=2,cex=0.9)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/individual_warbler_map_noarrows_closeup_randpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[-c(which(plot.inds=="TU")),1]),max(target.coords[-c(which(plot.inds=="TU")),1])),
					ylim=c(min(target.coords[-c(which(plot.inds=="TU")),2]),max(target.coords[-c(which(plot.inds=="TU")),2])),
					xlab="long",
					ylab="lat")
				text(target.coords[-c(which(plot.inds=="TU")),],
						labels=plot.inds[-c(which(plot.inds=="TU"))],
						col=adjustcolor(inds.col[-c(which(plot.inds=="TU"))],0.8),
						font=2,cex=0.9)
	dev.off()

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/individual_warbler_map_noarrows_closeup_nugget_randpr1.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(min(target.coords[-c(which(plot.inds=="TU")),1]),max(target.coords[-c(which(plot.inds=="TU")),1])),
					ylim=c(min(target.coords[-c(which(plot.inds=="TU")),2]),max(target.coords[-c(which(plot.inds=="TU")),2])),
					xlab="long",
					ylab="lat")
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
	k <- last.params$k
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
	
	best <- which.max(Prob)
	target.coords <- procrusteez(globe.coords,population.coordinates[[best]][1:k,],k,option=1)
	source.coords <- procrusteez(globe.coords,population.coordinates[[best]][1:k,],k,source.locs=population.coordinates[[best]][(k+1):(2*k),],option=2) 

	globe.admix.plot.cols <- continent.col

	require(maps)
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globe_world_map_dots.png",res=300,width=9*300,height=5.5*300)
	#quartz(width=9,height=5.5)
	map("world")
		box(lwd=2)
		points(globe.coords,pch=20,col=continent.col,cex=2)
			legend(x = -175,y=-10,
					legend = c("Africa","Western Eurasia","Central Eurasia","Eastern Eurasia","Oceania","Americas"),
					text.col = c("forestgreen","blue","purple","red","brown","orange"),cex=0.8,
					title="Continent plotting colors",title.col=1)
	dev.off()
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globe_world_map_text.png",res=300,width=9*300,height=5.5*300)
	#quartz(width=9,height=5.5)
	map("world")
		box(lwd=2)
		text(globe.coords,pops,col=continent.col,cex=0.5,font=2)
			legend(x = -175,y=-10,
					legend = c("Africa","Western Eurasia","Central Eurasia","Eastern Eurasia","Oceania","Americas"),
					text.col = c("forestgreen","blue","purple","red","brown","orange"),cex=0.8,
					title="Continent plotting colors",title.col=1)
	dev.off()
	
	x.min <- min(target.coords[,1]) - 5
	x.max <- max(target.coords[,1]) + 5
	y.min <- min(target.coords[,2]) - 5
	y.max <- max(target.coords[,2]) + 5
	
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globe_NoAd_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(x.min,x.max),
					ylim=c(y.min,y.max),
					xlab="long",
					ylab="lat")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,0.8),
						font=2,cex=0.8)
			box(lwd=2)
			legend(x = "bottomright",pch=NA,
					legend = c("Africa","Western Eurasia","Central Eurasia","Eastern Eurasia","Oceania","Americas"),
					text.col = c("forestgreen","blue","purple","red","brown","orange"),
					title="Continent plotting colors",title.col=1)
	dev.off()
	
	x.min <- min(target.coords[africa,1]) - 5
	x.max <- max(target.coords[africa,1]) + 5
	y.min <- min(target.coords[africa,2]) - 5
	y.max <- max(target.coords[africa,2]) + 5

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globe_Africa_NoAd_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords[africa,],type='n',
					xlim=c(x.min,x.max),
					ylim=c(y.min,y.max),
					xlab="long",
					ylab="lat")
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globe_Eurasia_Americas_Oceania_NoAd_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords[-africa,],type='n',
					xlim=c(x.min,x.max),
					ylim=c(y.min,y.max),
					xlab="long",
					ylab="lat")
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globe_Eurasia_NoAd_map.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords[-c(africa,americas,oceania),],type='n',
					xlim=c(x.min,x.max),
					ylim=c(y.min,y.max),
					xlab="long",
					ylab="lat")
				text(target.coords[-c(africa,americas,oceania),],
						labels=pops[-c(africa,americas,oceania)],
						col=adjustcolor(continent.col[-c(africa,americas,oceania)],0.8),
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

	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globe_NoAd_dist_compare.png",res=200,height=5*200,width=12*200)
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
					title="Continent plotting colors",title.col=1,
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



load("~/Desktop/Dropbox/space.mix/data/globetrotter/globe_spacemix/globe_spaceruns/rand_prior2/globe_spaceruns_randpr1_LongRun/globe_spaceruns_randpr1space_MCMC_output1.Robj")
	globe.admix.plot.cols <- fade.admixture.source.points(continent.col,admix.proportions[,best])
	png(file="~/Desktop/Dropbox/space.mix/ms/figs/globe_map_arrows.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim = c(18,67),
					ylim = c(-20,50),
					xlab="long",
					ylab="lat")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,0.8),
						font=2,cex=0.8)
				points(source.coords[,1],
						source.coords[,2],
							col=globe.admix.plot.cols,
							pch=20)
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=last.params$admix.proportions,
					length=0.1)
			box(lwd=2)
			legend(x = "bottomright",pch=NA,
					legend = c("Africa","Western Eurasia","Central Eurasia","Eastern Eurasia","Oceania","Americas"),
					text.col = c("forestgreen","blue","purple","red","brown","orange"),
					title="Continent plotting colors",title.col=1)
	dev.off()

	png(file="africa_map_arrows.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(0,20),
					ylim=c(-10,10),
					xlab="long",
					ylab="lat")
				text(target.coords[c(1:k),],
						labels=pops,
						col=adjustcolor(continent.col,0.8),
						font=2,cex=0.9)
				points(source.coords[,1],
						source.coords[,2],
							col=globe.admix.plot.cols,
							pch=20)
			arrows(	x0 = source.coords[,1],
					y0 = source.coords[,2],
					x1 = target.coords[,1],
					y1 = target.coords[,2],
					col=globe.admix.plot.cols,
					lwd=last.params$admix.proportions,
					length=0.1)
	dev.off()

	png(file="Nafrica_and_europe_map_arrows.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
		x <- c(africa,europe)
			plot(target.coords[x,],type='n',
					xlim=c(40.85,48.1),
					ylim=c(29.5,31.3),
					xlab="long",
					ylab="lat")
				text(target.coords[x,],
						labels=pops[x],
						col=adjustcolor(continent.col[x],0.8),
						font=2,cex=0.9)
				points(source.coords[x,1],
						source.coords[x,2],
							col=globe.admix.plot.cols[x],
							pch=20)
			arrows(	x0 = source.coords[x,1],
					y0 = source.coords[x,2],
					x1 = target.coords[x,1],
					y1 = target.coords[x,2],
					col=globe.admix.plot.cols[x],
					lwd=last.params$admix.proportions[x],
					length=0.1)
	dev.off()

	png(file="Nafrica_and_europe_map_noarrows.png",res=300,width=7*300,height=5*300,pointsize=9)
		x <- c(africa,europe)
			plot(target.coords[x,],type='n',
					xlim=c(40.85,48.1),
					ylim=c(29.5,31.3),
					xlab="long",
					ylab="lat")
				text(target.coords[x,],
						labels=pops[x],
						col=adjustcolor(continent.col[x],0.8),
						font=2,cex=0.9)
	dev.off()

	png(file="europe_map_arrows.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
			plot(target.coords,type='n',
					xlim=c(45.7,47.75),
					ylim=c(30.35,31.35),
					xlab="long",
					ylab="lat")
				text(target.coords[europe,],
						labels=pops[europe],
						col=adjustcolor(continent.col[europe],0.8),
						font=2,cex=0.9)
				points(source.coords[europe,1],
						source.coords[europe,2],
							col=globe.admix.plot.cols[europe],
							pch=20)
			arrows(	x0 = source.coords[europe,1],
					y0 = source.coords[europe,2],
					x1 = target.coords[europe,1],
					y1 = target.coords[europe,2],
					col=globe.admix.plot.cols[europe],
					lwd=last.params$admix.proportions[europe],
					length=0.1)
	dev.off()

	png(file="eurasia_map_arrows.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
		x <- c(europe,middle.muddle,east.asia)
			plot(target.coords[x,],type='n',
					xlim=c(44,54),
					ylim=c(29,42),
					xlab="long",
					ylab="lat")
				text(target.coords[x,],
						labels=pops[x],
						col=adjustcolor(continent.col[x],0.8),
						font=2,cex=0.9)
				points(source.coords[x,1],
						source.coords[x,2],
							col=globe.admix.plot.cols[x],
							pch=20)
			arrows(	x0 = source.coords[x,1],
					y0 = source.coords[x,2],
					x1 = target.coords[x,1],
					y1 = target.coords[x,2],
					col=globe.admix.plot.cols[x],
					lwd=last.params$admix.proportions[x],
					length=0.1)
	dev.off()

	png(file="eurasia_map_arrows.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
		x <- c(europe,middle.muddle,east.asia)
			plot(target.coords[x,],type='n',
					xlim=c(44,54),
					ylim=c(29,42),
					xlab="long",
					ylab="lat")
				text(target.coords[x,],
						labels=pops[x],
						col=adjustcolor(continent.col[x],0.8),
						font=2,cex=0.9)
				points(source.coords[x,1],
						source.coords[x,2],
							col=globe.admix.plot.cols[x],
							pch=20)
			arrows(	x0 = source.coords[x,1],
					y0 = source.coords[x,2],
					x1 = target.coords[x,1],
					y1 = target.coords[x,2],
					col=globe.admix.plot.cols[x],
					lwd=last.params$admix.proportions[x],
					length=0.1)
	dev.off()

	png(file="eurasia_map_somearrows.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
		x <- c(europe,middle.muddle,east.asia)
		y <- match(c("Hazara","Uzbekistani","Uygur"),pops)
			plot(target.coords[x,],type='n',
					xlim=c(44,54),
					ylim=c(29,42),
					xlab="long",
					ylab="lat")
				text(target.coords[x,],
						labels=pops[x],
						col=adjustcolor(continent.col[x],0.8),
						font=2,cex=0.9)
				points(source.coords[y,1],
						source.coords[y,2],
							col=globe.admix.plot.cols[y],
							pch=20)
			arrows(	x0 = source.coords[y,1],
					y0 = source.coords[y,2],
					x1 = target.coords[y,1],
					y1 = target.coords[y,2],
					col=globe.admix.plot.cols[y],
					lwd=last.params$admix.proportions[y],
					length=0.1)
	dev.off()

	png(file="eastasia_map_arrows.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
		x <- c(east.asia,middle.muddle)
			plot(target.coords[x,],type='n',
					xlim=c(48.7,53.8),
					ylim=c(38.2,41),
					xlab="long",
					ylab="lat")
				text(target.coords[x,],
						labels=pops[x],
						col=adjustcolor(continent.col[x],0.8),
						font=2,cex=0.9)
				points(source.coords[x,1],
						source.coords[x,2],
							col=globe.admix.plot.cols[x],
							pch=20)
			arrows(	x0 = source.coords[x,1],
					y0 = source.coords[x,2],
					x1 = target.coords[x,1],
					y1 = target.coords[x,2],
					col=globe.admix.plot.cols[x],
					lwd=last.params$admix.proportions[x],
					length=0.1)
	dev.off()

	
	png(file="eastasia_oceania_americas_map_arrows.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
		x <- c(east.asia,oceania,americas)
			plot(target.coords[x,],type='n',
					xlim=c(42,66),
					ylim=c(35,50),
					xlab="long",
					ylab="lat")
				text(target.coords[x,],
						labels=pops[x],
						col=adjustcolor(continent.col[x],0.8),
						font=2,cex=0.9)
				points(source.coords[x,1],
						source.coords[x,2],
							col=globe.admix.plot.cols[x],
							pch=20)
			arrows(	x0 = source.coords[x,1],
					y0 = source.coords[x,2],
					x1 = target.coords[x,1],
					y1 = target.coords[x,2],
					col=globe.admix.plot.cols[x],
					lwd=last.params$admix.proportions[x],
					length=0.1)
	dev.off()
	
	png(file="eurasia_oceania_americas_map_arrows.png",res=300,width=7*300,height=5*300,pointsize=9)
		#quartz(width=7,height=5,pointsize=9)
		x <- c(europe,middle.muddle,east.asia,oceania,americas,africa)
			plot(target.coords[x,],type='n',
					xlim=c(38,74),
					ylim=c(25,50),
					xlab="long",
					ylab="lat")
				text(target.coords[x,],
						labels=pops[x],
						col=adjustcolor(continent.col[x],0.8),
						font=2,cex=0.9)
			points(source.coords[x,1],
						source.coords[x,2],
							col=globe.admix.plot.cols[x],
							pch=20)
			arrows(	x0 = source.coords[x,1],
					y0 = source.coords[x,2],
					x1 = target.coords[x,1],
					y1 = target.coords[x,2],
					col=globe.admix.plot.cols[x],
					lwd=last.params$admix.proportions[x],
					length=0.1)
	dev.off()

################################
#	SIM FIGS
################################

################
#	Simulation scenarios
################
source("~/Desktop/Dropbox/space.mix/sims/spacemix_ms_sims.R")
png(file="~/Desktop/Dropbox/space.mix/ms/figs/basic_lattice.png",res=200,width=6*200,height=5*200)
	migration.rate.graphic(x.pops = 5,y.pops = 6,migration.rate=1,jitter=0.25,labels=TRUE,colors=TRUE)
dev.off()

png(file="~/Desktop/Dropbox/space.mix/ms/figs/barrier_lattice.png",res=200,width=6*200,height=5*200)
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
	
png(file="~/Desktop/Dropbox/space.mix/ms/figs/expansion_lattice.png",res=200,width=6*200,height=5*200)
	migration.rate.graphic(x.pop=5,y.pops=6,migration.rate=1,jitter = 0.25,expansion.list=expansion.list,labels=TRUE,colors=TRUE)
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
pdf("~/Desktop/Dropbox/space.mix/ms/figs/GeoGenMap_lattice.pdf",width=6,height=5,pointsize=9)
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
pdf("~/Desktop/Dropbox/space.mix/ms/figs/GeoGenMap_barrier.pdf",width=6,height=5,pointsize=9)
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
pdf("~/Desktop/Dropbox/space.mix/ms/figs/GeoGenMap_expansion.pdf",width=6,height=5,pointsize=9)
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
png("~/Desktop/Dropbox/space.mix/ms/figs/GeoGenMap_corner_admixture.png",res=300,width=6*300,height=5*300,pointsize=9)
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

png("~/Desktop/Dropbox/space.mix/ms/figs/GeoGenMap_corner_admixture_CYOL.png",res=300,width=6*300,height=5*300,pointsize=9)
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
png("~/Desktop/Dropbox/space.mix/ms/figs/GeoGenMap_corner_admixture_adinf.png",res=300,width=6*300,height=5*300,pointsize=9)
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
#	Barrier w/ Admixture
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
png("~/Desktop/Dropbox/space.mix/ms/figs/GeoGenMap_barr_inland_admixture_1.png",res=300,width=6*300,height=5*300,pointsize=9)
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

png("~/Desktop/Dropbox/space.mix/ms/figs/GeoGenMap_barr_inland_admixture_2.png",res=300,width=6*300,height=5*300,pointsize=9)
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














