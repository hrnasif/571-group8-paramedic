#!/bin/bash

#### Our model
# Under correctly specified model, vary data parameters
# Vary q, q_obs=7, N_samp=10 (defaults)
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=25
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=50
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=200

# Vary q_obs, q=100 (fixed), N_samp=10 (default)
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --q_obs=5
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --q_obs=10
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --q_obs=15
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --q_obs=20
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --q_obs=50

# Vary N_samp, q=10, q_obs=7 (defaults)
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=5
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=5
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=20
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=50
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=100

# Some combinations
# Vary N_samp, q=100 (fixed), q_obs=7 (default)
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=5 --q=100
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=5 --q=100
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=20 --q=100
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=50 --q=100
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=100 --q=100

# Vary sim-name (misspecify model)
# Default: N_subj=5, N_samp=10, q=100, q_obs=7
Rscript model_misspec_sim.R --sim-name="misspec-gamma-gamma-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100
Rscript model_misspec_sim.R --sim-name="misspec-halft-gamma-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100
Rscript model_misspec_sim.R --sim-name="misspec-gamma-halft-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100
Rscript model_misspec_sim.R --sim-name="misspec-halft-halft-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100
Rscript model_misspec_sim.R --sim-name="misspec-normal-gamma-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100
Rscript model_misspec_sim.R --sim-name="misspec-normal-halft-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100
Rscript model_misspec_sim.R --sim-name="misspec-gamma-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100
Rscript model_misspec_sim.R --sim-name="misspec-halft-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100
Rscript model_misspec_sim.R --sim-name="misspec-normal-normal-poisson-dirimult" --corr_within=1 --iter=10000 --warmup=5000 --q=100
Rscript model_misspec_sim.R --sim-name="misspec-normal-normal-negbin-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100
Rscript model_misspec_sim.R --sim-name="misspec-normal-normal-negbin-dirimult" --corr_within=1 --iter=10000 --warmup=5000 --q=100

#### Amy's model
# Vary sim-name to change data generating process (misspecify model)
# Default: N_subj=5, N_samp=10, q=10, q_obs=7
Rscript model_misspec_sim.R --sim-name="misspec-gamma-gamma-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-halft-gamma-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-gamma-halft-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-halft-halft-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-normal-gamma-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-normal-halft-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-gamma-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-halft-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-normal-normal-poisson-dirimult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-normal-normal-negbin-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-normal-normal-negbin-dirimult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=1
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=1

# Under correctly specified model, vary data parameters
# Vary q, q_obs=7, N_samp=10 (defaults)
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=25 --N_samp=1
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=50 --N_samp=1
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --N_samp=1
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=200 --N_samp=1

# Vary q_obs, q=100 (fixed), N_samp=10 (default)
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --q_obs=5 --N_samp=1
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --q_obs=10 --N_samp=1
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --q_obs=15 --N_samp=1
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --q_obs=20 --N_samp=1
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --q_obs=50 --N_samp=1

# Vary N_samp, q=10, q_obs=7 (defaults)
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=5 --N_samp=1
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=5 --N_samp=1
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=20 --N_samp=1
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=50 --N_samp=1
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=100 --N_samp=1

# Some combinations
# Vary N_samp, q=100 (fixed), q_obs=7 (default)
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=5 --q=100
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=5 --q=100
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=20 --q=100
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=50 --q=100
Rscript model_misspec_sim.R --sim-name="spec-normal-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --N_samp=100 --q=100

# Vary sim-name (misspecify model)
# Default: N_subj=5, N_samp=10, q=100, q_obs=7
Rscript model_misspec_sim.R --sim-name="misspec-gamma-gamma-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-halft-gamma-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-gamma-halft-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-halft-halft-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-normal-gamma-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-normal-halft-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-gamma-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-halft-normal-poisson-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-normal-normal-poisson-dirimult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-normal-normal-negbin-mult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --N_samp=1
Rscript model_misspec_sim.R --sim-name="misspec-normal-normal-negbin-dirimult" --corr_within=1 --iter=10000 --warmup=5000 --q=100 --N_samp=1