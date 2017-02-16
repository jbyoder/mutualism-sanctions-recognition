# parameter input for simulations_par_soft2.R
# 2016.12.27

# NOTE: the following must be accounted for before source()ing this file
# reps = vector of replicate ID numbers
# alfreqs = function from sim_functions.R, source()'d in simulations_par.R

# parameters of the metapopulation
Npops = rep(50, length(reps), replace=T) # npops
K1s = sample(20:200, length(reps), replace=T) # per-pop carrying capacity for species 1 (host)
K2s = sample(200:2000, length(reps), replace=T) # ditto, species 2 (symbiont)
m1s = runif(length(reps), 0, 0.05) # migration rates
m2s = runif(length(reps), 0, 0.05)

# genetic parameters
p11s = alfreqs(length(reps)) # starting allele freq at loc 1 for sp 1
p12s = alfreqs(length(reps))
p21s = alfreqs(length(reps)) 
p22s = alfreqs(length(reps)) 
rec1s = runif(length(reps), 0, 0.5) # recombination rate, sp 1
rec2s = runif(length(reps), 0, 0.5)
mu1s = mu2s = rep(1e-6, length(reps))  # mutation rate for sp 1, 2

# interaction terms
# new this iteration: back to the original configuration, because that issue with cost/benefits was wrong (see project notes)
oms = runif(length(reps), 0.01, 1) # sanctions effectiveness (necessary for S, SR)
C1s = runif(length(reps), 0.01, 0.1) # cost to sp 1
B1s = C1s + runif(length(reps), 0.01, 0.1) # create ratios from 2 to 10 for host
C2s = runif(length(reps), 0.001, 0.1)
B2s = C2s + runif(length(reps), 0.01, 1) # ratios from 11 to 1000 for symbiont



