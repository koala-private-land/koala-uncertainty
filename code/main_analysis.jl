# Main analysis runs - four strategies


# 1. Selected properties (sp)
# 2. Number of climate scenarios in each resolution of uncertainty (ns)
# 3. Time when uncertainty is resolved] (tt)
# 4. Koala index threshold (kt)

import Pkg
Pkg.activate(".")
Pkg.instantiate()

using CSV
using DataFrames
using Revise
using StatsBase
using Gurobi
using JuMP
using Distributions
using JLD2
using Random

Random.seed!(54815861)

# Constants
R = 30; # Total number of uncertainty realisations
rr = 30; # Sampled realisations
S = 12;


parallel_mode = false

include("sensitivity_analysis_func.jl")

dir = "results/model_runs"

println("Starting main analysis...")

sp = 1; # 1:100
tt = 6; # [1,2,3,4,5,6,7,8], 2000 - 2070
kt = 0.25; # [0.1, 0.15, 0.2, 0.25, 0.3]
ns = 12; # 1:12
sdr = 0.02;
deforestation_risk = 0.1
k = 7000;
K_pa_change = 0.0 # Change of biodiversity target per-annum (n.b. 0.024 equivalent to doubling target over 30 years)
ssb = 1.0e7 # Second-stage budget (in NPV terms)

# Ignore uncertainty
fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt, kt, ns, sdr, deforestation_risk, k, true, ssb, K_pa_change)

# Inflexible (nr)/ Restricted Flexibility (pr)/ Flexible (ar)
fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt, kt, ns, sdr, deforestation_risk, k, false, ssb, K_pa_change)

# Flexible & Learning
fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt, kt, 1, sdr, deforestation_risk, k, false, ssb, K_pa_change)