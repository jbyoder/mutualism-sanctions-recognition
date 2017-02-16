# Simulated host-symbiont coevolution, under different genetic models
# Assumes environment with multiple computing cores accessible
# jby 2017.01.05
#
# 2016.04.19 - version posted to GitHub for first bioRxiv submission
# 2016.11.27 - fixed critical bug in which requests for SR sims ran R sims
# 2016.12.02 - sourcing new functions and parameters for soft selection
# 2016.12.24 - new implementation of soft selection (mating prob based on fitness)
# 2016.12.30 - tweaks and optimization for better parallel work ...
# 2017.01.05 - and a bug fix ...

# run this as
# R --vanilla --args {cores} {sim} {rep nos} {record interval} {par file} {output file} < Rscripts/simulations_par.R
# cores = number of processes to run
# sim = N (none) | S (sanctions) | R (recognition) | SR (both)
# rep nos = range of replicate i.d. numbers, as 1-100; 101-200, &c
# record interval = time between write-out, in simulation generations
# par file = path to file with parameters for the sims
# output file = path and (base) filename for files recording simulation results
#
# e.g.: R --vanilla --args 4 S 1-100 50 Rscripts/sim_par_example.R ~/scratch/data/Symbiosis_theory/example_sim < Rscripts/simulations_par.R

setwd("~/Symbiosis_theory")
#setwd("~/Documents/Academic/Active_projects/Symbiosis_theory")

.libPaths("~/shared/Rscripts/packages")
#require("dplyr", lib.loc="~/Rpackages")
require("parallel")
require("foreach")
require("iterators")
require("doParallel", lib.loc="~/shared/Rscripts/packages")

source("Rscripts/sim_functions_soft.R")

# initialize the cluster
cl = makeCluster(as.numeric(commandArgs(trailingOnly=T)[1]))
registerDoParallel(cl)
print(paste("Made cluster with ", commandArgs(trailingOnly=T)[1] , "cores"))
print(cl)

clusterEvalQ(cl, require("parallel"))
clusterEvalQ(cl, require("foreach"))
clusterEvalQ(cl, require("iterators"))
clusterEvalQ(cl, require("doParallel", lib.loc="~/shared/Rscripts/packages"))

#-------------------------------------------------------------------------
# run setup

# "N" (none) | "S" (sanctions) | "R" (recognition) | "SR" (both)
intype = commandArgs(trailingOnly=T)[2]
cat(intype,"\n")
# intype = "S" # for interactivity

repends = as.numeric(strsplit(commandArgs(trailingOnly=T)[3],"-")[[1]])
reps = repends[1]:repends[2] # number of replicate simulations
# reps = 1:4

gens = as.numeric(commandArgs(trailingOnly=T)[4]) # time to simulate
gensout = seq(0, gens, as.numeric(commandArgs(trailingOnly=T)[5]))
# gens = 100
# gensout = seq(0, gens, 10)

# source the external parameters file
source(commandArgs(trailingOnly=T)[6])

# pull in the outfile prefix
outpre = commandArgs(trailingOnly=T)[7]

#------------------------------------------------
# a function to run one sim

Sim <- function(intype,gens,r,Npop,K1,K2,m1,m2,p11,p12,p21,p22,rec1,rec2,mu1,mu2,C1,B1,C2,B2,om)
{

# create the dataframe for results
TS = data.frame(intype=intype, 
	repl=r,
	Npop,K1,K2,m1,m2,mu1,mu2,C1,B1,C2,B2,rec1,rec2,mu1,mu2,om, 
	gen=rep(0:gens, each=2*(Npop+1)), 
	pop=rep(rep(as.character(c(1:Npop,"tot")), each=2),(gens+1)), 
	sp=rep(1:2, (Npop+1)*(gens+1)), 
	p1=NA, p2=NA, D=NA)

# starting populations
sp1 = make(Npop, K1, p11, p12)
sp2 = make(Npop, K2, p21, p22)

# keep track of what goes where ...
sp1.p = TS$pop!="tot" & TS$sp==1
sp1.t = TS$pop=="tot" & TS$sp==1
sp2.p = TS$pop!="tot" & TS$sp==2
sp2.t = TS$pop=="tot" & TS$sp==2

# starting frequencies
TS[TS$gen==0 & sp1.p, c("p1","p2")] = do.call("rbind", lapply(sp1,freqs))
TS[TS$gen==0 & sp1.t, c("p1","p2")] = freqs(do.call("cbind",sp1))
TS[TS$gen==0 & sp2.p, c("p1","p2")] = do.call("rbind", lapply(sp2,freqs))
TS[TS$gen==0 & sp2.t, c("p1","p2")] = freqs(do.call("cbind",sp2))

# starting LD
TS[TS$gen==0 & sp1.p, "D"] <- unlist(lapply(sp1, function(x) Dprime(as.character(x[1,]), as.character(x[2,]))))
TS[TS$gen==0 & sp1.t, "D"] = Dprime(as.character(do.call("cbind",sp1)[1,]), as.character(do.call("cbind",sp1)[2,]))
TS[TS$gen==0 & sp2.p, "D"] = unlist(lapply(sp2, function(x) Dprime(as.character(x[1,]), as.character(x[2,]))))
TS[TS$gen==0 & sp2.t, "D"] = Dprime(as.character(do.call("cbind",sp2)[1,]), as.character(do.call("cbind",sp2)[2,]))

# !! need to add record of each replicate's parameters !!

# generational recursion ------------------------
for(g in 1:gens){
	
	# migration among sites --- try doing this post-selection?
	sp1 = migrate(sp1, m1)
	sp2 = migrate(sp2, m2)

	# selection
	sp2 = lapply(sp2, function(x) x[,sample(1:K2,K1)]) # pick symbionts that interact
	
	if(intype=="S") Ws = mapply(function(x,y) selS(x, y, om, C1, B1, C2, B2, getfits=T), sp1, sp2, SIMPLIFY=F)
	if(intype=="R") Ws = mapply(function(x,y) selR(x, y, C1, B1, C2, B2, getfits=T), sp1, sp2, SIMPLIFY=F)
	if(intype=="SR") Ws = mapply(function(x,y) selSR(x, y, om, C1, B1, C2, B2, getfits=T), sp1, sp2, SIMPLIFY=F)
	
	# random mating and reproduction
	if(intype!="N"){
		sp1 = mapply(function(x,y) mate(x, K1, rec1, mu1, fit=y), sp1, lapply(Ws, function(w) w[1,]), SIMPLIFY=F)
		sp2 = mapply(function(x,y) mate(x, K2, rec2, mu2, fit=y), sp2, lapply(Ws, function(w) w[2,]), SIMPLIFY=F)
	}
	if(intype=="N"){
		sp1 = lapply(sp1, mate, K1, rec1, mu1, fit=rnorm(K1,1,0.1)) # non-uniform prob (test)
		sp2 = lapply(sp2, mate, K2, rec2, mu2, fit=rnorm(K1,1,0.1))
	}
	
	# record allele frequencies
	TS[TS$gen==g & sp1.p, c("p1","p2")] = do.call("rbind", lapply(sp1,freqs))
	TS[TS$gen==g & sp1.t, c("p1","p2")] = freqs(do.call("cbind",sp1))
	TS[TS$gen==g & sp2.p, c("p1","p2")] = do.call("rbind", lapply(sp2,freqs))
	TS[TS$gen==g & sp2.t, c("p1","p2")] = freqs(do.call("cbind",sp2))
	
	# record LD
	TS[TS$gen==g & sp1.p, "D"] = unlist(lapply(sp1, function(x) Dprime(as.character(x[1,]), as.character(x[2,]))))
	TS[TS$gen==g & sp1.t, "D"] = Dprime(as.character(do.call("cbind",sp1)[1,]), as.character(do.call("cbind",sp1)[2,]))
	TS[TS$gen==g & sp2.p, "D"] = unlist(lapply(sp2, function(x) Dprime(as.character(x[1,]), as.character(x[2,]))))
	TS[TS$gen==g & sp2.t, "D"] = Dprime(as.character(do.call("cbind",sp2)[1,]), as.character(do.call("cbind",sp2)[2,]))

} # end generational recursion

return(TS)
}


#-------------------------------------------------------------------------
# LOOP over replicates

# distribute necessary stuff to the cluster
clusterExport(cl, c("Sim", "intype", "gens", "reps", "Npops", "K1s", "K2s", "m1s", "m2s", "p11s", "p12s", "p21s", "p22s", "rec1s", "rec2s", "mu1s", "mu2s", "C1s", "B1s", "C2s", "B2s", "oms", "outpre"))

strt<-Sys.time()
# run the sims! ---------------------------------
foreach(i=1:length(reps)) %dopar% {

TS = Sim(intype, gens, reps[i], Npops[i], K1s[i], K2s[i], m1s[i], m2s[i], p11s[i], p12s[i], p21s[i], p22s[i], rec1s[i], rec2s[i], mu1s[i], mu2s[i], C1s[i], B1s[i], C2s[i], B2s[i], oms[i])

if(reps[i] < 10) rlab = paste("000", reps[i], sep="")
if(reps[i] >= 10 & i < 100) rlab = paste("00", reps[i], sep="")
if(reps[i] >= 100 & i < 1000) rlab = paste("0", reps[i], sep="")
if(reps[i] >= 1000) rlab = reps[i]

write.table(TS[TS$gen %in% gensout, ], file=paste(outpre, "_", intype,"_r", rlab, ".txt", sep=""), col.names=T, sep="\t")

print(paste("Finished rep", rlab, "and wrote to", paste(outpre, "_", intype,"_r", rlab, ".txt", sep="")))

} # end recursion over replicate sims

stopCluster(cl)

q()
