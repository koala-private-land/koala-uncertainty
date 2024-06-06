# Full sensitivity runs

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
R = 1; # Total number of uncertainty realisations
rr = 1; # Sampled realisations
S = 12;

parallel_mode = false

include("sensitivity_analysis_func.jl")

## Start sensitivity analysis

# Sensitivity parameters (base), whilst avoid rerunning base parameters again in individual iteration
sp = 1; # 1:100
sp_vec = 2:10
tt = 6; # [1,2,3,4,5,6,7,8], 2000 - 2070
tt_vec = 3:5
kt = 0.25; # [0.1, 0.15, 0.2, 0.25, 0.3]
kt_vec = [0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.275]
ns = 12; # 1:12
ns_vec = [1,3]
sdr = 0.02;
sdr_vec = [-99.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05] # -99 for green_book declining discount rates
deforestation_risk = 0.1
deforestation_risk_vec = [1.0, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.0]
k = 7000;
k_vec = [2500, 4000, 5500, 7000, 8500, 10000, 11500]
K_pa_change = 0.0 # Change of biodiversity target per-annum (n.b. 0.024 equivalent to doubling target over 30 years)
K_pa_change_vec = [0.0, 0.005, 0.01, 0.015, 0.02, 0.025]
ssb = 1.0e7 # Second-stage budget (in NPV terms)
ssb_vec = [0.0, 0.5e7, 1.0e7, 1.5e7, 2.0e7, 2.5e7, 3.0e7, 3.5e7, 4.0e7]

dir = "results/model_runs"

println("Starting sensitivity analysis...")

# Run across all discount rate choices
map((sdr_i) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt, kt, 1, sdr_i, deforestation_risk, k, false, ssb, K_pa_change), sdr_vec)
map((sdr_i) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt, kt, ns, sdr_i, deforestation_risk, k, false, ssb, K_pa_change), sdr_vec)

# Run across all values of K_pa_change, change in target per annum
map((K_pa_change) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt, kt, ns, sdr, deforestation_risk, k, false, ssb, K_pa_change), K_pa_change_vec)
map((K_pa_change) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt, kt, 1, sdr, deforestation_risk, k, false, ssb, K_pa_change), K_pa_change_vec)

# Run across all values of second-stage budget constraints
map((ssb) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt, kt, ns, sdr, deforestation_risk, k, false, ssb, K_pa_change), ssb_vec)
map((ssb) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt, kt, 1, sdr, deforestation_risk, k, false, ssb, K_pa_change), ssb_vec)

# Run across all values of k (conservation goal in area)
map((k) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt, kt, 1, sdr, deforestation_risk, k, false, ssb, K_pa_change), k_vec)
map((k) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt, kt, ns, sdr, deforestation_risk, k, false, ssb, K_pa_change), k_vec)

# Run across all deforestation risk rates
map((dr_i) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt, kt, 1, sdr, dr_i, k, false, ssb, K_pa_change), deforestation_risk_vec)
map((dr_i) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt, kt, ns, sdr, dr_i, k, false, ssb, K_pa_change), deforestation_risk_vec)

# Run across all choices of the discount rate
map((sdr_i) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt, kt, ns, sdr_i, deforestation_risk, k, false, ssb, K_pa_change), sdr_vec)

# Run across different time periods for tt
map((tt_i) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt_i, kt, 1, sdr, deforestation_risk, k, false, ssb, K_pa_change), tt_vec)
map((tt_i) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt_i, kt, ns, sdr, deforestation_risk, k, false, ssb, K_pa_change), tt_vec)

# Run across all kt
map((kt_i) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt, kt_i, 1, sdr, deforestation_risk, k, false, ssb, K_pa_change), kt_vec)
map((kt_i) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp, tt, kt_i, ns, sdr, deforestation_risk, k, false, ssb, K_pa_change), kt_vec)

# Run across different sampled properties
map((sp_i) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp_i, tt, kt, 1, sdr, deforestation_risk, k, false, ssb, K_pa_change), sp_vec)
map((sp_i) -> fcn_run_optim(kitl_index_full, stratified_samples, dir, sp_i, tt, kt, ns, sdr, deforestation_risk, k, false, ssb, K_pa_change), sp_vec)
println("Sensitivity analysis complete.")