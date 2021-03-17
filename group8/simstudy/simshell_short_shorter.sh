#!/bin/bash

#### Our model
# Under correctly specified model, vary data parameters
# Vary q_obs, q=10, N_samp=10 (defaults)
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=8000 --q_obs=3
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=8000 --q_obs=5
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=8000 --q_obs=7

# Vary N_samp, q=10, q_obs=7 (defaults)
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=8000 --N_samp=5
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=8000 --N_samp=15
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=8000 --N_samp=20
#Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=8000 --N_samp=50
#Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=8000 --N_samp=100

# Some combinations
# Vary N_samp, q=15, q_obs=7 (defaults)
#Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=8000 --N_samp=5 --q=15
#Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=8000 --N_samp=20 --q=15
#Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=8000 --N_samp=50 --q=15
#Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=8000 --N_samp=100 --q=15

# Vary sim-name (misspecify model)
# Default: N_subj=5, N_samp=10, q=15, q_obs=7
#Rscript model_misspec_sim.R --sim-name="misspec-gamma-gamma-poisson-mult" --corr_within=1 --iter=10000 --warmup=8000
#Rscript model_misspec_sim.R --sim-name="misspec-normal-gamma-poisson-mult" --corr_within=1 --iter=10000 --warmup=8000 
#Rscript model_misspec_sim.R --sim-name="misspec-gamma-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=8000 
