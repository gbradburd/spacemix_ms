##################################################################
##################################################################
##	run ms on a lattice
##################################################################
##################################################################

# the dimensions of the lattice to be simulated.
#	returns list of coordinates of sampled and un-sampled demes
generate.geographic.coordinates <- function(x.pops,y.pops){
	# recover()
	total.nrow <- 2*y.pops + 1
	total.ncol <- 2*x.pops + 1
	grid.breaks.long <- matrix(seq(from=0,to=total.nrow-1,by=1),
						   nrow = total.ncol,
						   ncol = total.nrow,
						   byrow=TRUE)
	grid.breaks.lat <- matrix(seq(from=0,to=total.ncol-1,by=1),
						   nrow = total.ncol,
						   ncol = total.nrow,
						   byrow=FALSE)
	samp.and.unsamp.coords <- cbind(as.vector(grid.breaks.long), as.vector(grid.breaks.lat))
	samp.coords <- cbind(as.vector(grid.breaks.long[seq(2,total.ncol,2),seq(2,total.nrow,2)]), 
						as.vector(grid.breaks.lat[seq(2,total.ncol,2),seq(2,total.nrow,2)]))
	all.coords <- list("samp.coords" = samp.coords,
						"samp.and.unsamp.coords" = samp.and.unsamp.coords)
	return(all.coords)					
}

get.barrier.distance <- function(pop.coordinates,barrier.line){
	barrier.distance <- matrix(0,nrow=nrow(pop.coordinates),ncol=nrow(pop.coordinates))
	for(i in 1:nrow(pop.coordinates)){
		for(j in 1:nrow(pop.coordinates)){
			if(pop.coordinates[i,1] < barrier.line && pop.coordinates[j,1] > barrier.line){
					barrier.distance[i,j] <- 1
					barrier.distance[j,i] <- 1
			}
		}
	}
	return(barrier.distance)
}

#figure out which populations are on either side of a barrier down
#	the longitudinal center of the lattice
get.barrier.pops <- function(pop.coordinates,all.coords.dist,barrier.line){
	# recover()
	barrier.distance <- get.barrier.distance(pop.coordinates,barrier.line)
	barrier.index <- barrier.distance==1 & is.finite(all.coords.dist)
	return(barrier.index)
}

#this is the beating heart of the script.
#	It generates the migration rate in ms format 
#	given the specified lattice and
#	the simulation scenario
write.migration.rates <- function(x.pops,y.pops,migration.rate,barrier.effect=NULL,admixture.event=NULL,source=NULL,target=NULL,admixture.time.point=NULL,admixture.proportion=NULL,expansion.event=NULL,expansion.list=NULL){
	# recover()
	pop.coordinates <- generate.geographic.coordinates(x.pops,y.pops)
	all.coords.dist <- fields::rdist(pop.coordinates$samp.and.unsamp.coords)
	all.coords.dist[which(all.coords.dist > sqrt(2) | all.coords.dist < 1)] <- Inf
	if(!is.null(barrier.effect)){
		barrier.line <- diff(range(pop.coordinates$samp.and.unsamp.coords[,1]))/2 + 1e-10
		barrier.index <- get.barrier.pops(pop.coordinates$samp.and.unsamp.coords,all.coords.dist,barrier.line)
		all.coords.dist[barrier.index] <- all.coords.dist[barrier.index]*barrier.effect
	}
	pairwise.migration.matrix <- migration.rate/all.coords.dist
	migration.index <- which(pairwise.migration.matrix != 0,arr.ind=TRUE)
	migration.rate.vector <- c()
	for(i in 1:nrow(migration.index)){
		migration.rate.vector <- c(migration.rate.vector,
									sprintf("-m %s %s %s",
										migration.index[i,1],
										migration.index[i,2],
										pairwise.migration.matrix[migration.index[i,1],migration.index[i,2]]))
	}
	if(!is.null(admixture.event)){
		admixture.call <- write.admixture.event.call(source = source,target = target,
														npops = nrow(all.coords.dist),
														time.point = admixture.time.point,
														admixture.proportion = admixture.proportion)
		migration.rate.vector <- c(migration.rate.vector,admixture.call)
	}
	if(!is.null(expansion.event)){
		expansion.call <- write.expansion.event.call(expansion.list)
		migration.rate.vector <- c(migration.rate.vector,expansion.call)
	}
	return(migration.rate.vector)
}

#outputs the ms call for an admixture event
write.admixture.event.call <- function(source,target,npops,time.point,admixture.proportion){
# recover()
	admixture.call <- c(sprintf("-es %s %s %s",
										time.point,
										target,
										1 - admixture.proportion),
							sprintf("-ej %s %s %s",
										time.point + 0.000001,
										npops+1,
										source))
	return(admixture.call)
}

#outputs the ms call for an expansion event
write.expansion.event.call <- function(expansion.list){
	# recover()
	expansion.call <- c()
	for(i in 1:length(expansion.list)){
		for(j in 1:length(expansion.list[[i]]$daughters)){
			expansion.call <- c(expansion.call,
									sprintf("-ej %s %s %s",
									expansion.list[[i]]$time.point,
									expansion.list[[i]]$daughters[j],
									expansion.list[[i]]$parent))
		}
	}
	return(expansion.call)
}

#outputs the graphic for the specified scenario
migration.rate.graphic <- function(x.pops,y.pops,migration.rate,migration.arrows=TRUE,jitter,barrier.effect=NULL,labels=NULL,expansion.list=NULL,colors=FALSE,xlim=NULL,ylim=NULL,curve=0.5,arrow.width=1,arrow.col="gray40",pop.lab.cex=1,pop.pt.cex=1){
	# recover()
	iArrows <- igraph:::igraph.Arrows
	pop.coordinates <- generate.geographic.coordinates(x.pops,y.pops)
	all.coords.dist <- fields::rdist(pop.coordinates$samp.and.unsamp.coords)
	all.coords.dist[which(all.coords.dist > sqrt(2) | all.coords.dist < 1)] <- Inf
	if(!is.null(barrier.effect)){
		barrier.line <- diff(range(pop.coordinates$samp.and.unsamp.coords[,1]))/2 + 1e-10
		barrier.index <- get.barrier.pops(pop.coordinates$samp.and.unsamp.coords,all.coords.dist,barrier.line)
		all.coords.dist[barrier.index] <- all.coords.dist[barrier.index]*barrier.effect
	}
	if(colors){
		k <- x.pops * y.pops
		pop.cols <- rainbow(k,start=4/6,end=6/6)[as.numeric(cut(pop.coordinates[[1]][,1],k))]
	}
	pairwise.migration.matrix <- migration.rate/all.coords.dist
	par(mar=c(1,1,1,1))
	plot(pop.coordinates[[2]],main="",xlab="",ylab="",type='n',xlim=xlim,ylim=ylim,xaxt='n',yaxt='n')
	if(!colors){
		points(pop.coordinates[[2]]) ; points(pop.coordinates[[1]],col=2,pch=19,cex=2*pop.pt.cex)
	} else {
		points(pop.coordinates[[2]]) ; points(pop.coordinates[[1]],col=pop.cols,pch=19,cex=2.1*pop.pt.cex,lwd=2)
	}
	if(migration.arrows){
		for(i in 1:nrow(pop.coordinates[[2]])){
			for(j in 1:nrow(pop.coordinates[[2]])){
				if(pairwise.migration.matrix[i,j] != 0){
					arrow.coords <- arrow.graphic(pop.coordinates[[2]][i,1],pop.coordinates[[2]][i,2],pop.coordinates[[2]][j,1],pop.coordinates[[2]][j,2],jitter)
							arrows( x0 = arrow.coords$x0,
									y0 = arrow.coords$y0,
									x1 = arrow.coords$x1,
									y1 = arrow.coords$y1,
									col="purple",
									lwd=pairwise.migration.matrix[i,j],
									length=0.05,
									code=3)
				}
			}
		}
	}
	if(!is.null(labels)){
		text(pop.coordinates[[1]],labels=paste(1:(x.pops*y.pops)),col="white",cex=0.7*pop.lab.cex,font=2)
	}
			if(!is.null(expansion.list)){
			for(i in 1:length(expansion.list)){
				for(j in 1:length(expansion.list[[i]]$daughters)){
					iArrows(x1 = pop.coordinates[[2]][expansion.list[[i]]$parent,1],
							y1 = pop.coordinates[[2]][expansion.list[[i]]$parent,2], 
							x2 = pop.coordinates[[2]][expansion.list[[i]]$daughters[j],1],
							y2 = pop.coordinates[[2]][expansion.list[[i]]$daughters[j],2],
							code = 2,
							curve = curve,
							sh.col = arrow.col,
							h.col = arrow.col,
				            h.lwd = 2.5*arrow.width,
				            sh.lwd = 2.5*arrow.width,
							width = 0.8,
							size = 0.7)
				}
			}
		}
	box(lwd=2)
}

#another way to visualize pairwise Fst between the populations,
#	to check if you've simulated what you think you've simulated
fst.pop.graphic <- function(x.pops,y.pops,fst.mat,jitter,barrier.effect=NULL,labels=NULL){
	# recover()
	tmp.fst.mat <- -log(fst.mat)
	pop.coordinates <- generate.geographic.coordinates(x.pops,y.pops)[[1]]
	all.coords.dist <- fields::rdist(pop.coordinates)
	all.coords.dist[which(all.coords.dist > 2 | all.coords.dist < 1)] <- Inf
	plot(pop.coordinates,main="",xlab="",ylab="")
		for(i in 1:nrow(pop.coordinates)){
			for(j in 1:nrow(pop.coordinates)){
				if(is.finite(all.coords.dist[i,j])){
					arrow.coords <- arrow.graphic(pop.coordinates[i,1],pop.coordinates[i,2],pop.coordinates[j,1],pop.coordinates[j,2],jitter)
							arrows( x0 = arrow.coords$x0,
									y0 = arrow.coords$y0,
									x1 = arrow.coords$x1,
									y1 = arrow.coords$y1,
									col="purple",
									lwd=tmp.fst.mat[i,j],
									length=0.05,
									code=3)
				}
			}
		}
	if(!is.null(labels)){
		text(pop.coordinates,labels=paste(1:(x.pops*y.pops)),col="white",cex=0.7)
	}
}

# function for drawing pairwise migration arrows
#	on a population lattice
arrow.graphic <- function(x0,y0,x1,y1,jitter){
	coords <- list("x0" = x0, "y0" = y0, "x1" = x1, "y1" = y1)
	if(x0 == x1 | y0 == y1){
		if(x0 == x1){
			if(y0 < y1){
				coords$y0 <- y0 + jitter
				coords$y1 <- y1 - jitter
			} else {
				coords$y0 <- y0 - jitter
				coords$y1 <- y1 + jitter		
			}
		}
		if(y0 == y1){
			if(x0 < x1){
				coords$x0 <- x0 + jitter
				coords$x1 <- x1 - jitter
			} else {
				coords$x0 <- x0 - jitter
				coords$x1 <- x1 + jitter
			}
		}
	} else {
		if(y0 < y1){
			coords$y0 <- y0 + jitter
			coords$y1 <- y1 - jitter
		} else {
			coords$y0 <- y0 - jitter
			coords$y1 <- y1 + jitter		
		}
		if(x0 < x1){
			coords$x0 <- x0 + jitter
			coords$x1 <- x1 - jitter
		} else {
			coords$x0 <- x0 - jitter
			coords$x1 <- x1 + jitter
		}
	}
	return(coords)
}

# function that feeds ms the number of chromosomes to 
#	be sampled in each sampled deme
write.lattice <- function(dim.row,dim.col,nChromo){
	lattice <- matrix(0,nrow=dim.row*2+1,ncol=dim.col*2+1)
		for(i in 1:(dim.row*2+1)){
			for(j in 1:(dim.col*2+1)){
				if(i%%2 == 0 && j%%2 == 0){
					lattice[i,j] <- nChromo
				}
			}
		}
	return(lattice)
}

# code cannibalized from Dan Denison that reads ms output into R
read.ms.haplotype.matrices <- function(nsam, ndraws, ms.output.file) {
    txt <- scan(file=ms.output.file, what=character(0), sep="\n", quiet=TRUE)
    h <- list()
    ## THE OUTPUT TEXT FOR EACH DRAW SHOULD CONTAIN THE WORD "segsites"
    marker <- grep("segsites", txt)
    stopifnot(length(marker) == ndraws)
    ## GET NUMBERS OF SEGREGATING SITES IN EACH DRAW
    segsites <- sapply(strsplit(txt[marker], split=":"), function(vec) as.integer(vec[2]))
    for(draw in seq(along=marker)) {
        if(!(draw %% 100)) cat(draw, " ")
        if(segsites[draw] > 0) {
            haplotypes <- txt[(marker[draw] + 2):(marker[draw] + 2 + nsam - 1)]
            haplotypes <- strsplit(haplotypes, split="")
            h[[draw]] <- sapply(haplotypes, function(el) c(as.integer(el)))
            ## IF THERE'S 1 SEGREGATING SITE, THIS WON'T BE A MATRIX 
            if(segsites[draw] == 1) h[[draw]] <- as.matrix(h[[draw]])
            ## OTHERWISE, IT NEEDS TO BE TRANSPOSED
            else h[[draw]] <- t(h[[draw]])
        }
        else h[[draw]] <- matrix(nrow=nsam, ncol=0)
        stopifnot(all(dim(h[[draw]]) == c(nsam, segsites[draw])))  
    }
    cat("\n")
    h
}

# function that actually calls ms with the call put together by all these R functions
ms <- function(x.pops,y.pops,nChromo,theta,migration.rate,barrier.effect=NULL,admixture.event=NULL,source=NULL,target=NULL,admixture.time.point=NULL,admixture.proportion=NULL,expansion.event=NULL,expansion.list=NULL){
	# recover()
		ms.output.file <- "ms_output"
		random.seeds <- c(sample(1:100000,3,replace=TRUE))
		sampled.pops <- x.pops * y.pops
		total.pops <- (2 * x.pops + 1) * (2 * y.pops + 1)
		call <- paste(
					sprintf(
						"./ms %s 1 -t %s -s 1 -I %s %s 0.0-m %s -seeds %s %s %s",
							sampled.pops*nChromo,
							theta,
							total.pops,
							paste(write.lattice(x.pops,y.pops,nChromo),collapse = " "),
							paste(write.migration.rates(x.pops,y.pops,migration.rate,barrier.effect,admixture.event,source,target,admixture.time.point,admixture.proportion,expansion.event,expansion.list),collapse = " "),
							random.seeds[1],
							random.seeds[2],
							random.seeds[3]
					),
					"-T", ">", ms.output.file
				)
		cat(call,file="ms_call.txt")
		system(call)
		read.ms.haplotype.matrices(nsam=sampled.pops*nChromo,ndraws=1,ms.output.file=ms.output.file)
}

# sanity check on the ms command line values
	check.ms.command.line.values <- function(ms.theta,ms.r,ms.nsites,ms.m,ms.barrier.time,locus.size=ms.nsites,diploid.population.size=1000,per.bp.mu=1.5e-08,per.bp.rho=1.5e-08){
		mu						<- ms.theta/(4*diploid.population.size*locus.size)
		rho						<- ms.r/((ms.nsites-1)*4*diploid.population.size)
		migration.fraction		<- ms.m/(4*diploid.population.size)
		barrier.arrival.time 	<- ms.barrier.time/(4*diploid.population.size)
		cat(
			sprintf(
				"\t\t\t\tmu is %s\r
				rho is %s\r
				migration.fraction is %s\r
				barrier.arrival.time is %s\r",
				mu,rho,migration.fraction,barrier.arrival.time)
		)
	}

# generate times and migration rates in ms units
	generate.ms.command.line.values <- function(diploid.population.size,locus.size,per.bp.mu,migration.fraction){ #,generations.ago
			ms.command.line.values <- vector("list",length=2)
				names(ms.command.line.values) <- c("t","m") #,"time"
					ms.command.line.values$t <- 4*diploid.population.size*per.bp.mu*locus.size
					ms.command.line.values$m <- 4*diploid.population.size*migration.fraction
					# ms.command.line.values$time <- 4*diploid.population.size*generations.ago 
					
					#admixture on time scale more recent than 1/(4Nm k) won't spread
			return(ms.command.line.values)
	}

# make dataset for use by, e.g., spacemix
generate.spacemix.dataset <- function(x.pops,y.pops,nLoci,nChromo,theta,migration.rate,barrier.effect=NULL,admixture.event=NULL,source=NULL,target=NULL,admixture.time.point=NULL,admixture.proportion=NULL,expansion.event=NULL,expansion.list=NULL){
	# recover()
	#Allele Counts
		data.matrix <- do.call(cbind,replicate(nLoci,ms(x.pops,y.pops,nChromo,theta,migration.rate,barrier.effect,admixture.event,source,target,admixture.time.point,admixture.proportion,expansion.event,expansion.list)))
		sampled.pops <- x.pops*y.pops
		population.membership <- c()
			for(i in 1:sampled.pops){
				population.membership <- c(population.membership,rep(i,nChromo))
			}
		allele.counts <- matrix(0,nrow=sampled.pops,ncol=nLoci)
		for(i in 1:sampled.pops){
			allele.counts[i,] <- colSums(data.matrix[which(population.membership==i),])
		}
	#Sample sizes
		sample.sizes <- matrix(nChromo,nrow=sampled.pops,ncol=nLoci)

	#Geographic distance matrix
		population.coordinates <- generate.geographic.coordinates(x.pops,y.pops)
		D <- fields::rdist(population.coordinates[[1]])
		spacemix.dataset <- list("allele.counts" = allele.counts,
							"sample.sizes" = sample.sizes,
							"D" = D,
							"population.coordinates" = population.coordinates[[1]])
		if(!is.null(barrier.effect)){
			barrier.line <- diff(range(population.coordinates[[2]][,1]))/2 + 1e-10
			E <- get.barrier.distance(population.coordinates[[1]],barrier.line)
			spacemix.dataset <- list("allele.counts" = allele.counts,
							"sample.sizes" = sample.sizes,
							"D" = D,
							"barrier.line" = barrier.line,
							"E" = E,
							"population.coordinates" = population.coordinates[[1]])
		}
	return(spacemix.dataset)
}
