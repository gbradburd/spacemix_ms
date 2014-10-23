################################################################
################################################################
#	make figures for SpaceMix paper
################################################################
################################################################
 

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

################################
#	SIM FIGS
################################

################
#	Grid
################
load("~/Desktop/Dropbox/space.mix/sims/stationary_pops/spacemix/stationary_pops_1/spacemix_ms_sim_stationary_pops_1space_MCMC_output1.Robj")
load("~/Desktop/Dropbox/space.mix/sims/stationary_pops/spacemix/stationary_pops_1/spacemix.ms.dataset_stationary_pops.Robj")
k <- last.params$k
best <- which.max(Prob)
target.coords <- procrusteez(obs.locs = spacemix.dataset$population.coordinates,
							target.locs = population.coordinates[[best]][1:k,],
							k = k,
							option = 1)
pdf("~/Desktop/Dropbox/space.mix/ms/figs/GeoGenMap_lattice.pdf",width=6,height=5,pointsize=9)
#quartz(width=6,height=5)
par(mar=c(4.3,4.3,3,1))
plot(target.coords,pch=1,xlim=c(0,12),ylim=c(-0.4,10),cex=3.5,xlab="Eastings",ylab="Northings",main="Inferred Population Map:\n Simple Lattice")
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
plot(target.coords,pch=1,xlim=c(0,12),ylim=c(-0.4,10),cex=3.5,xlab="Eastings",ylab="Northings",main="Inferred Population Map:\n Lattice with Barrier")
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
plot(target.coords,pch=1,xlim=c(1.2,9.5),ylim=c(-2,12),cex=3.5,xlab="Eastings",ylab="Northings",main="Inferred Population Map:\n Lattice with Expansion")
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
source.coord.cols <- fade.admixture.source.points(rep(1,k),admix.proportions[,best])
png("~/Desktop/Dropbox/space.mix/ms/figs/GeoGenMap_corner_admixture.png",res=300,width=6*300,height=5*300,pointsize=9)
#quartz(width=6,height=5)
par(mar=c(4.3,4.3,3,1))
plot(target.coords,xlim=c(x.min,x.max),ylim=c(y.min,y.max),pch=1,cex=3.5,xlab="Eastings",ylab="Northings",main="Inferred Population Map:\n Lattice with Admixture")
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
plot(target.coords,xlim=c(x.min,x.max),ylim=c(y.min,y.max),pch=1,cex=3.5,xlab="Eastings",ylab="Northings",main="Inferred Population Map:\n Lattice with Admixture")
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
source.coord.cols <- fade.admixture.source.points(rep(1,k),admix.proportions[,best])
png("~/Desktop/Dropbox/space.mix/ms/figs/GeoGenMap_corner_admixture_adinf.png",res=300,width=6*300,height=5*300,pointsize=9)
#quartz(width=6,height=5)
par(mar=c(4.3,4.3,3,1))
plot(target.coords,xlim=c(x.min,x.max),ylim=c(y.min,y.max),pch=1,cex=3.5,xlab="Eastings",ylab="Northings",main="Inferred Population Map:\n Lattice with Admixture")
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
source.coord.cols <- fade.admixture.source.points(rep(1,k),scalar*admix.proportions[,best])
png("~/Desktop/Dropbox/space.mix/ms/figs/GeoGenMap_barr_inland_admixture_1.png",res=300,width=6*300,height=5*300,pointsize=9)
#quartz(width=6,height=5)
par(mar=c(4.3,4.3,3,1))
plot(target.coords,xlim=c(x.min,x.max),ylim=c(y.min,y.max),pch=1,cex=3.5,xlab="Eastings",ylab="Northings",main="Inferred Population Map:\n Lattice with Barrier and Admixture")
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
plot(source.coord.hdr2,shadecols=adjustcolor(1,0.2),show.points=FALSE,outside.points=FALSE,
		xlab="Eastings",ylab="Northings",main="Inferred Population Map:\n Lattice with Barrier and Admixture",
		xlim=c(x.min,x.max),
		ylim=c(y.min,y.max))
points(target.coords,xlim=c(x.min,x.max),ylim=c(y.min,y.max),pch=1,cex=3.5)
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














source("~/Desktop/Dropbox/space.mix/sims/spacemix_ms_sims.R")
load("~/Desktop/Dropbox/space.mix/sims/expansion/spacemix/noad/rand_prior1/sim_expansion_dataset.Robj")
load("~/Desktop/Dropbox/space.mix/sims/expansion/spacemix/noad/rand_prior1/expansion_randpr_noad_1_LongRun/expansion_randpr_noad_1space_MCMC_output1.Robj")

parents <- c(78:88)
time.points <- rep(0.07,11)
expansion.list <- vector(mode="list")
	for(i in 1:length(parents)){
		expansion.list[[i]] <- list(parent=parents[i],
									daughters = parents[i]+c(11,22,33,44,55),
									time.point = time.points[i])
	}

png(file="~/desktop/test.png",res=200,width=6*200,height=5*200,pointsize=9)
	migration.rate.graphic(x.pop=5,y.pops=6,migration.rate=1,jitter = 0.2,expansion.list=expansion.list,labels=TRUE)
dev.off()	