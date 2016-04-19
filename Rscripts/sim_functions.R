# Functions for simulated host-symbiont coevolution, under different genetic models
# jby 2015.08.11


#-------------------------------------------------------------------------
# functions

# to generate starting allele frequencies from an approximation of the equilibrium SFS
# 
alfreqs <- function(Nsam, prec=20){
	N = prec         # make sure N is even!
	ac = 1:(N-1)     # if sample size = 100, 

	#now implement eqn 6 from the paper (Ganapathy and Uyenoyama 2009) 
	unfolded_SFS = (1/ac) / (sum(1/ac))
	#if one wants the foldedSFS
	#folded_SFS = c(unfolded_SFS[1:(N/2-1)]+unfolded_SFS[(N-1):(N/2+1)],unfolded_SFS[N/2])
			
	return(sample(1:(N-1)/N, Nsam, replace=T, prob=unfolded_SFS))
}

# to generate Npop metapopulations of 2-locus x K haploids
make <- function(Npop=10, K=100, p1=0.5, p2=0.5){
	return(lapply(vector(mode="list", length=Npop), function(x) rbind(sample(1:0, K, replace=T, prob=c(p1,1-p1)), sample(1:0, K, replace=T, prob=c(p2,1-p2))) ))
} # validated!

# to randomly draw m*k individuals per pop, pool them, drop them into new pops
migrate <- function(mpop, m){
	k = ncol(mpop[[1]])
	np = length(mpop)
	
	tot = do.call("cbind", mpop) # make mpop into a nice looooong matrix
	kt = ncol(tot)
	Nmt = m*kt
	
	migs = sample(1:kt, Nmt) # ID migrants
	tot[,migs] = tot[,sample(migs,length(migs))] # extract, randomize, reinsert migrants
	
	return(lapply(1:np-1, function(i) tot[,(i*k+1):((i+1)*k)]))
} # validated!

# fitness outcomes from sanctions
selS <- function(ho, pa, om, Ch, Bh, Cs, Bs, tune=0.05){
	h = ho[1,] # relevant host genotypes
	p = pa[1,] # ditto symbionts
	ps = sort(sample(1:length(p), length(h), replace=T)) # sample symbionts (unsampled die)
	
	# interaction outcome matrices
	Hmat = 1 + matrix(c(Bh-Ch, Bh-Ch, -(1-om)*Ch, -Ch), 2, 2, dimnames=list(1:0, 1:0))
	Pmat = 1 + matrix(c(Bs-Cs, (1-om)*Bs, Bs-Cs, Bs), 2, 2, dimnames=list(1:0, 1:0))
	
	fits = mapply(function(x,y) c(Hmat[x,y], Pmat[y,x]), as.character(h), as.character(p[ps])) # interact pops
	
	# hard selection: is fit >= rnorm(mean=1,sd=0.1)
	# (but keep at least 5% randomly drawn individuals)
	sp1=ho[,c(which(fits[1,]>=rnorm(ncol(fits),1,tune)),sample(1:ncol(fits),0.05*ncol(fits)))]
	sp2=pa[,c(which(fits[2,]>=rnorm(ncol(fits),1,tune)),sample(1:ncol(fits),0.05*ncol(fits)))]
	
	return(list(sp1=sp1, sp2=sp2))
}

# fitness outcomes from recognition
selR <- function(ho, pa, Ch, Bh, Cs, Bs, tune=0.05){
	h = ho[1,] # relevant host genotypes
	p = apply(pa,2,function(x) paste(x,collapse="")) 
	ps = sort(sample(1:length(p), length(h), replace=T)) # sample symbionts (unsampled die)
	
	# interaction outcome matrices
	Hmat = 1 + matrix(c(Bh-Ch, 0, 0, Bh-Ch, -Ch, 0, 0, -Ch), 2, 4, dimnames=list(1:0,c(11,10,"01","00")))
	Pmat = 1 + matrix(c(Bs-Cs, 0, 0, Bs-Cs, Bs, 0, 0, Bs), 2, 4, dimnames=list(1:0,c(11,10,"01","00")))
	
	fits = mapply(function(x,y) c(Hmat[x,y], Pmat[x,y]), as.character(h), p[ps]) # interact pops
	
	# hard selection: is fit >= rnorm(mean=1,sd=0.1)
	# (but keep at least 5% randomly drawn individuals)
	sp1=ho[,c(which(fits[1,]>=rnorm(ncol(fits),1,tune)),sample(1:ncol(fits),0.05*ncol(fits)))]
	sp2=pa[,c(which(fits[2,]>=rnorm(ncol(fits),1,tune)),sample(1:ncol(fits),0.05*ncol(fits)))]
			
	return(list(sp1=sp1, sp2=sp2))
}

selSR <- function(ho, pa, om, Ch, Bh, Cs, Bs, tune=0.05){
	h = apply(ho,2,function(x) paste(x,collapse="")) # relevant host genotypes
	p = apply(pa,2,function(x) paste(x,collapse="")) 
	ps = sort(sample(1:length(p), length(h), replace=T)) # sample symbionts (unsampled die)
	
	# interaction outcome matrices
	Hmat = 1 + matrix(c(Bh-Ch, 0, Bh-Ch, 0, 0, Bh-Ch, 0, Bh-Ch, -(1-om)*Ch, 0, -Ch, 0, 0, -(1-om)*Ch, 0, -Ch), 4, 4, dimnames=list(c(11,10,"01","00"),c(11,10,"01","00")))
	Pmat = 1 + matrix(c(Bs-Cs, 0, Bs-Cs, 0, 0, Bs-Cs, 0, Bs-Cs, (1-om)*Bs, 0, Bs, 0, 0, (1-om)*Bs, 0, Bs), 4, 4, dimnames=list(c(11,10,"01","00"),c(11,10,"01","00")))
	
	fits = mapply(function(x,y) c(Hmat[x,y], Pmat[x,y]), h, p[ps]) # interact pops
	
	# hard selection: is fit >= rnorm(mean=1,sd=0.1)
	# (but keep at least 5% randomly drawn individuals)
	sp1=ho[,c(which(fits[1,]>=rnorm(ncol(fits),1,tune)),sample(1:ncol(fits),0.05*ncol(fits)))]
	sp2=pa[,c(which(fits[2,]>=rnorm(ncol(fits),1,tune)),sample(1:ncol(fits),0.05*ncol(fits)))]
		
	return(list(sp1=sp1, sp2=sp2))
}

# mate random pairs of individuals w/ probability based on fitness
mate <- function(pop,k,rec,mu){
	# draw parents, up to carrying capacity k pairings
	p1s = sample(1:dim(pop)[2],k,replace=T)
	p2s = sample(1:dim(pop)[2],k,replace=T)
	
	# function to return offspring from one pairing
	m <- function(pa1, pa2){
		return(sample(list(pa1,pa2,c(pa1[1],pa2[2]),c(pa2[1],pa1[2])), 1, prob=c((1-rec)/2,(1-rec)/2,rec/2,rec/2))[[1]])
	}
	
	# offspring
	off = mapply(function(x,y) m(pop[,x],pop[,y]), p1s, p2s)
	
	# mutation
	off[off==1] <- sample(c(1,0),length(which(off==1)),replace=T,prob=c(1-mu,mu))
	off[off==0] <- sample(c(0,1),length(which(off==0)),replace=T,prob=c(1-mu,mu))
	
	return(off)
}


# calculate allele frequencies at all loci in a pop
freqs <- function(pop){
	return(apply(pop,1,function(x) sum(x)/length(x)))
}

# calculate LD between two loci
Dprime <- function(snp.one, snp.two) {
    # @param snp.one: character vector of bases for first position
    # @param snp.two: character vector of bases for 2nd position
    # @return: r-squared value between the two positions (numeric)

    # Check that character vectors are input
    if(!is.character(snp.one)) stop("snp.one not a character vector")
    if(!is.character(snp.two)) stop("snp.two not a character vector")
    
    snp.one[snp.one=="N"] = NA
    snp.two[snp.two=="N"] = NA

    # Calculate frequencies of each genotype combination:
    a<-table(as.data.frame(cbind(snp.one, snp.two)))
    
    # Check that SNPs are biallelic
    if(sum(dim(a)) < 4) return(NA)# SNPs are not biallelic
    if(sum(dim(a)) > 4) return(NA)# SNPs are not biallelic
    
    # Calculate D prime
    # D = freq(A1)*freq(B1) - p1*q1
    # Dmax = min(p1q1, p2q2) if D < 0
    # Dmax = min(p1q2, p2q1) if D > 0
    # D' = (D/(Dmax)
    D = (a[1,1] / sum(a)) - (sum(a[,1]) / sum(a)) * (sum(a[1,]) / sum(a))
	#if(D < 0) Dmax = min(c((sum(a[,1]) / sum(a)) * (sum(a[1,]) / sum(a)), (sum(a[,2]) / sum(a)) * (sum(a[2,]) / sum(a))))
	#if(D > 0) Dmax = min(c((sum(a[,1]) / sum(a)) * (sum(a[2,]) / sum(a)), (sum(a[,2]) / sum(a)) * (sum(a[1,]) / sum(a))))
	#if(D == 0) return(0)
    
    #Dpr = D/Dmax
    
    #if(Dpr=="NaN") rsq = NA # deals with a complete linkage glitch ...

    return(D)
}

# calculate R^2 between two loci
Rsq <- function(snp.one, snp.two) {
    # @param snp.one: character vector of bases for first position
    # @param snp.two: character vector of bases for 2nd position
    # @return: r-squared value between the two positions (numeric)

    # Check that character vectors are input
    if(!is.character(snp.one)) stop("snp.one not a character vector")
    if(!is.character(snp.two)) stop("snp.two not a character vector")
    
    snp.one[snp.one=="N"] = NA
    snp.two[snp.two=="N"] = NA

    # Calculate frequencies of each genotype combination:
    a<-table(as.data.frame(cbind(snp.one, snp.two)))
    
    # Check that SNPs are biallelic
    if(sum(dim(a)) < 4) return(NA)# SNPs are not biallelic
    if(sum(dim(a)) > 4) return(NA)# SNPs are not biallelic
    
    # Calculate r-squared
    # D = freq(A1)*freq(B1) - p1*q1
    # rsq = (D/(sqrt(p1*q1*p2*q2)))^2
    rsq = ((a[1,1] / sum(a) - (sum(a[,1]) / sum(a)) * (sum(a[1,]) / sum(a))) /
            sqrt( sum((a[,1]) / sum(a))*(sum(a[1,]) / sum(a)) * (sum(a[,2]) / sum(a)) * (sum(a[2,]) /sum(a))))^2
    
    if(rsq=="NaN") rsq = NA # deals with a complete linkage glitch ...

    return(rsq)
}
