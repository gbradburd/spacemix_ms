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

get.procrustes.locations.posterior.list <- function(observed.coords,population.coordinates.posterior){
	target.coords.list <- vector(mode="list",length = length(population.coordinates.posterior))
	source.coords.list <- vector(mode="list",length = length(population.coordinates.posterior))
	k <- nrow(observed.coords)
	for(i in 1:length(target.coords.list)){
		target.coords.list[[i]] <- procrusteez(obs.locs = observed.coords,
												target.locs = population.coordinates.posterior[[i]][1:k,],
												k = k,
												option = 1)
		source.coords.list[[i]] <- procrusteez(obs.locs = observed.coords,
												target.locs = population.coordinates.posterior[[i]][1:k,],
												k = k,
												source.locs = population.coordinates.posterior[[i]][(k+1):(2*k),],
												option = 2)
	}
	return(list(target.coords.list= target.coords.list, source.coords.list= source.coords.list))
}

fade.admixture.source.points <- function(pop.cols,admix.proportions){
	faded.colors <- numeric(length(pop.cols))
	for(i in 1:length(pop.cols)){
		faded.colors[i] <- adjustcolor(pop.cols[i],admix.proportions[i])
	}
	return(faded.colors)
}

plot.admix.arrows <- function(source.coords,target.coords,admix.proportions=1,colors=NULL,length=0.1){
	if(is.null(colors)){
		colors <- rep("black",nrow(source.coords))
	}
	arrows(x0 = source.coords[,1],
			y0 = source.coords[,2],
			x1 = target.coords[,1],
			y1 = target.coords[,2],
			lwd = admix.proportions,
			col = colors,
			length = length)
	return(invisible("arrows!"))
}

plot.spacemix.output.map <- function(input.file,obs.coords,sample.labels,colors,pdf.file.name,pdf.width,pdf.height,admixture=FALSE,xlim=NULL,ylim=NULL,pdf=TRUE){
	load(input.file)
	best <- which.max(Prob)
	k <- last.params$k
	target.coords <- procrusteez(obs.coords,population.coordinates[[best]][1:k,],k,option=1)
	source.coords <- procrusteez(obs.coords,population.coordinates[[best]][1:k,],k,source.locs=population.coordinates[[best]][(k+1):(2*k),],option=2)
	weighted.colors <- fade.admixture.source.points(colors,admix.proportions[,best])
	if(is.null(xlim)){
		xlim=c(range(target.coords[,1],source.coords[,1]))
	}
	if(is.null(ylim)){
		ylim=c(range(target.coords[,2 ],source.coords[,2]))
	}
	if(pdf){
		pdf(file=paste(pdf.file.name,".pdf",sep=""),width=pdf.width,height=pdf.height,pointsize=9)
	}
		par(mar=c(4.5,4.5,3,1))
			plot(target.coords,type='n',
					xlim=xlim,
					ylim=ylim,
					xlab="Eastings",
					ylab="Northings",cex.lab=1.7)
				text(target.coords[c(1:k),],
						labels= sample.labels,
						col=colors,
						font=2,cex=2)
			if(admixture){
				points(source.coords,
							col=weighted.colors,
							pch=20)
				plot.admix.arrows(source.coords,target.coords,admix.proportions[,best],colors=weighted.colors,length=0.1)
			}
				box(lwd=2)
	if(pdf){			
		dev.off()
	}
}

load_MCMC_output <- function(MCMC.output.file){
    tmpenv <- environment()
	tmp <- load(MCMC.output.file,envir=tmpenv)
	mcmc.output <- lapply(tmp,get,envir=tmpenv)
	names(mcmc.output) <- tmp
	return(mcmc.output)
}

get.posterior.location.matrix.from.list <- function(posterior.list,population.index){
	post.location.matrix <- matrix(unlist(
								lapply(posterior.list,
									FUN=function(elem){elem[population.index,]})),
								nrow=length(posterior.list),ncol=2,byrow=TRUE)
	return(post.location.matrix)
}

get.credible.ellipse <- function(posterior.points,quantile){
	require(MASS)
	require(cluster)
	fit <- cov.mve(posterior.points, quantile.used = nrow(posterior.points) * quantile)
	points_in_ellipse <- posterior.points[fit$best, ]
	ellipse_boundary <- predict(ellipsoidhull(points_in_ellipse))
	return(ellipse_boundary)
}

plot.credible.ellipse <- function(ellipse_boundary,population.color,fading=0.3,lty=1){
	polygon(ellipse_boundary,col=adjustcolor(population.color,fading),border=1,lty=lty)
}

make.spacemix.map.list <- function(MCMC.output.file,observed.coords,name.vector,color.vector,quantile=0.95){
	MCMC.output <- load_MCMC_output(MCMC.output.file)
	best <- which.max(MCMC.output$Prob)
	admix.source.color.vector <- fade.admixture.source.points(color.vector,rowMeans(MCMC.output$admix.proportions))
	k <- MCMC.output$last.params$k
	target.coords <- procrusteez(observed.coords,MCMC.output$population.coordinates[[best]][1:k,],k,option=1)
	source.coords <- procrusteez(observed.coords,MCMC.output$population.coordinates[[best]][1:k,],k,
									source.locs=MCMC.output$population.coordinates[[best]][(k+1):(2*k),],option=2)
	procrustes.coord.posterior.lists <- get.procrustes.locations.posterior.list(observed.coords= observed.coords,
																				population.coordinates.posterior=MCMC.output$population.coordinates)
	posterior.target.location.matrices <- lapply(1:k,get.posterior.location.matrix.from.list,posterior.list=procrustes.coord.posterior.lists$target.coords.list)
	posterior.source.location.matrices <- lapply(1:k,get.posterior.location.matrix.from.list,posterior.list=procrustes.coord.posterior.lists$source.coords.list)
	posterior.target.ellipses <- lapply(posterior.target.location.matrices,get.credible.ellipse,quantile)
	posterior.source.ellipses <- lapply(posterior.source.location.matrices,get.credible.ellipse,quantile)
	spacemix.map.list <- c(MCMC.output,
							list(observed.coords=observed.coords),
								list(name.vector=name.vector),list(color.vector=color.vector),
								list(quantile=quantile),list(source=source),list(best = best),
								list(admix.source.color.vector = admix.source.color.vector),
								list(k = k),list(target.coords = target.coords),list(source.coords = source.coords),
								list(procrustes.coord.posterior.lists = procrustes.coord.posterior.lists),
								list(posterior.target.location.matrices = posterior.target.location.matrices),
								list(posterior.source.location.matrices = posterior.source.location.matrices),
								list(posterior.target.ellipses = posterior.target.ellipses),
								list(posterior.source.ellipses = posterior.source.ellipses))
	return(spacemix.map.list)
}

make.spacemix.map <- function(spacemix.map.list,text=FALSE,ellipses=TRUE,source.option=TRUE,xlim=NULL,ylim=NULL){
	with(spacemix.map.list,{ 
		plot(target.coords,type='n',xlim=xlim,ylim=ylim,xlab="",ylab="")
			if(ellipses){
				lapply(1:k,FUN=function(i){plot.credible.ellipse(posterior.target.ellipses[[i]],color.vector[i])})
			}
			if(text){
				text(target.coords,col=color.vector,font=2,labels=name.vector,cex=0.7)
			}
			if(source.option){
				if(ellipses){
					lapply(1:k,FUN=function(i){plot.credible.ellipse(posterior.source.ellipses[[i]],admix.source.color.vector[i],fading=1,lty=2)})
				}
				text(source.coords,col= admix.source.color.vector,font=3,labels=name.vector,cex=0.7)
				plot.admix.arrows(source.coords,target.coords,
									admix.proportions=admix.proportions[,best],
									colors=admix.source.color.vector,length=0.1)
			}
				box(lwd=2)
	})
}

query.spacemix.map <- function(focal.pops,spacemix.map.list,source.option=TRUE){
	with(spacemix.map.list,{
		# browser()
		focal.indices <- match(focal.pops,name.vector)
			for(i in 1:length(focal.indices)){
				plot.credible.ellipse(posterior.target.ellipses[[focal.indices[i]]],color.vector[focal.indices[i]],fading=1)
			}
			if(source.option){
				for(i in 1:length(focal.indices)){
					plot.credible.ellipse(posterior.source.ellipses[[focal.indices[i]]], admix.source.color.vector[focal.indices[i]],fading=1,lty=2)
				}
				text(source.coords[focal.indices,,drop=FALSE],col=1,font=3,labels=name.vector[focal.indices])
				arrows(	x0 = source.coords[focal.indices,1],
						y0 = source.coords[focal.indices,2],
						x1 = target.coords[focal.indices,1],
						y1 = target.coords[focal.indices,2],
						col= admix.source.color.vector[focal.indices[i]],
						lwd=1,
						length=0.1)
			}
			text(target.coords[focal.indices,,drop=FALSE],col=1,font=2,labels=name.vector[focal.indices],cex=1)
				box(lwd=2)
	})
}

get.credible.interval <- function(param.matrix,pop.order){
	k <- nrow(param.matrix)
	cred.intervals <- lapply(pop.order,FUN=function(i){quantile(param.matrix[i,],c(0.025,0.975))})
	return(cred.intervals)
}

make.cred.bars <- function(quantile.vector,bar.width,color.vector,vert.line.width=NULL,pop.order=NULL){
	# recover()
	if(is.null(pop.order)){
		pop.order <- 1:length(quantile.vector)
	}
	if(is.null(vert.line.width)){
		vert.line.width <- 0.5
	}
	x.coord <- 1:k
	color.vector <- color.vector[pop.order]
	for(i in 1:length(quantile.vector)){
		lines(x = c(x.coord[i]-bar.width/2,x.coord[i]+bar.width/2),
				y = c(quantile.vector[[pop.order[i]]][1],quantile.vector[[pop.order[i]]][1]),col=color.vector[i])
		lines(x = c(x.coord[i]-bar.width/2,x.coord[i]+bar.width/2),
				y = c(quantile.vector[[pop.order[i]]][2],quantile.vector[[pop.order[i]]][2]),col=color.vector[i])
		lines(x = c(x.coord[i],x.coord[i]),
				y = quantile.vector[[pop.order[i]]],col=adjustcolor(color.vector[i],0.5),lwd= vert.line.width)
	}
}

calculate.pairwise.pi <- function(ind1,ind2){
	diff.homs = sum(ind1!=ind2 & abs(ind1-ind2)!=1 )
	hets = sum(ind1==1 | ind2 ==1 )
	return((diff.homs + hets/2)/length(ind1))
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

ff.text <- function ( xy, labels, rep.fact=4, attr.fact=0.2, col="black", text.cex=1, xlim, ylim, ... ) {
    # a try at expanding labels. 
    txy <- FFieldPtRep( xy, rep.fact=4 )
    if (missing(xlim)) { xlim <- range(xy[,1],txy[,1]) }
    if (missing(ylim)) { ylim <- range(xy[,2],txy[,2]) }
    plot( xy, pch=20, col=adjustcolor(col,.5), xlim=xlim, ylim=ylim, ... )
    text( txy, labels=pops, col=col, cex=text.cex, ... )
    segments( x0=xy[,1], x1=txy[,1], y0=xy[,2], y1=txy[,2], lwd=2, col=adjustcolor(col,.25) )
}

get.sample.covariance <- function(counts,sample.sizes){
	sample.frequencies <- counts/sample.sizes
	mean.sample.sizes <- rowMeans(sample.sizes)
	mean.sample.frequencies <- matrix(apply(sample.frequencies,2,
											get.weighted.mean.frequency,
											mean.sample.sizes=mean.sample.sizes),
									nrow=length(mean.sample.sizes),ncol=ncol(sample.frequencies),byrow=TRUE)
	normalized.sample.frequencies <- sample.frequencies/sqrt(mean.sample.frequencies*(1-mean.sample.frequencies))
	sample.covariance <- cov(t(normalized.sample.frequencies),use="pairwise.complete.obs")
	loci <- ncol(sample.frequencies)
	return(list("sample.covariance" = sample.covariance,
				"mean.sample.sizes" = mean.sample.sizes,
				"loci" = loci))
}

get.weighted.mean.frequency <- function(sample.frequencies,mean.sample.sizes){
	na.pops <- which(is.na(sample.frequencies))
	if(sum(na.pops) > 0){
		sample.frequencies <- sample.frequencies[-na.pops]
		mean.sample.sizes <- mean.sample.sizes[-na.pops]
	}
	weighted.sample.frequencies <- mean.sample.sizes*sample.frequencies
	sample.frequency.weighted.mean <- sum(weighted.sample.frequencies)/sum(mean.sample.sizes)
	return(sample.frequency.weighted.mean)
}