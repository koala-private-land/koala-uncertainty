# Full sensitivity runs

# 1. Selected properties (sp)
# 2. Number of climate scenarios in each resolution of uncertainty (ns)
# 3. Time when uncertainty is resolved (tt)
# 4. Koala index threshold (kt)

using CSV
using DataFrames
using Revise
using StatsBase
using Gurobi
using JuMP
using StatsPlots
using Distributions
using JLD2
using Random

Random.seed!(54815861)

include("optim-functions.jl")

cost_df = CSV.read("data/spatial_predictions_10yr.csv", DataFrame);
kitl_index_full = CSV.read("data/kitl_prop_climate.csv", DataFrame);
stratified_samples = CSV.read("data/stratified_sample.csv", DataFrame);
climateProjList = ["CCCMA_R1", "CCCMA_R2", "CCCMA_R3", "CSIRO_R1", "CSIRO_R2", "CSIRO_R3", "ECHAM_R1", "ECHAM_R2", "ECHAM_R3", "MIROC_R1", "MIROC_R2", "MIROC_R3"];

# Constants
R = 100; # Total number of uncertainty realisations
rr = 10; # Sampled realisations
K = 7000; # Habitat protection goal

mutable struct Solution
    feasible::Bool
    model::Model
    x::AbstractArray
    y::AbstractArray
    w::AbstractArray
    cost::AbstractArray
    metric::AbstractArray
end

function cost_to_ts(cost::AbstractVector, area::AbstractVector, adopt::AbstractVector, idx_before::AbstractVector=1:30, idx_after::AbstractVector=31:60, inflation_rate::AbstractFloat=0.02)
    delta = (1 + inflation_rate) .^ (cat(idx_before, idx_after, dims=1))
    delta_before = sum(delta[idx_before])
    delta_after = sum(delta[idx_after])
    cost_ts = (cost .* 1000 .* area .* adopt ./ 10) # Cost per-year
    cost_before = cost_ts .* delta_before # Cost start at 2020
    cost_after = cost_ts .* delta_after
    return (cost_before, cost_after)
end

function subset_data(cost_full, kitl_index_full, subset_id::Vector{Int64})
    cost_subset = filter(row -> row.NewPropID ∈ subset_id, cost_full)
    kitl_subset = filter(row -> row.NewPropID ∈ subset_id, kitl_index_full)
    filter!(row -> row.NewPropID ∈ kitl_subset.NewPropID, cost_subset)
    filter!(row -> row.NewPropID ∈ cost_subset.NewPropID, kitl_subset)
    return (cost_subset, kitl_subset)
end

function metric_to_input(metric::DataFrame, area::AbstractVector, tt::Integer, kt::AbstractFloat)
    # Formats metrics to input
    t_symbols = [:t0, :t1, :t2, :t3, :t4, :t5, :t6, :t7]
    N = length(unique(metric.NewPropID))
    climateProjList = ["CCCMA_R1", "CCCMA_R2", "CCCMA_R3", "CSIRO_R1", "CSIRO_R2", "CSIRO_R3", "ECHAM_R1", "ECHAM_R2", "ECHAM_R3", "MIROC_R1", "MIROC_R2", "MIROC_R3"]
    S = length(climateProjList)
    M₁ = zeros(N, length(1:(tt-1)), S)
    M₂ = zeros(N, length(tt:length(t_symbols)), S)
    for t = eachindex(t_symbols)
        metric_unstacked = unstack(metric, :NewPropID, :climate_model, t_symbols[t])
        metric_mat = Matrix(metric_unstacked[:, climateProjList])
        metric_threshold_mat = (metric_mat .> kt)
        if (t < tt)
            M₁[:, t, :] = metric_threshold_mat .* area
        else
            M₂[:, t-tt+1, :] = metric_threshold_mat .* area
        end
    end
    return (M₁, M₂)
end

function realisation_sample(cost_subset::DataFrame, kitl_subset::DataFrame, tt::Integer=4, kt::AbstractFloat=0.25)
    (M₁, M₂) = metric_to_input(kitl_subset, vec(cost_subset.AREA), tt, kt)
    MeanAdopt = cost_subset.MeanAdopt
    SDAdopt = cost_subset.SDAdopt
    MeanProp = cost_subset.MeanProp
    SDProp = cost_subset.SDProp
    MeanWTA = cost_subset.MeanWTA
    SDWTA = cost_subset.SDWTA
    adopt_binary = (rand.(Normal.(MeanAdopt, SDAdopt))) .> rand(length(MeanAdopt))
    prop = (rand.(Normal.(MeanProp, SDProp))) .* adopt_binary
    (C₁, C₂) = cost_to_ts(MeanWTA + SDWTA .* randn(length(SDWTA)), cost_subset.AREA, cost_subset.MeanAdopt, 1:(10*(tt-1)), (1+10*(tt-1)):60, 0.02)
    sample_C₁ = C₁ .* prop
    sample_C₂ = C₂ .* prop
    sample_M₁ = M₁ .* prop
    sample_M₂ = M₂ .* prop
    area = vec(cost_subset.AREA) .* prop
    return Realisation(sample_C₁, sample_C₂, sample_M₁, sample_M₂, area, adopt_binary)
end

function fcn_tidy_two_stage_solution_sum(solution, prop_id)
    out = hcat(prop_id, solution.x, sum(solution.y, dims=2), sum(solution.w, dims=2))
    colnames = vcat(["NewPropID"; "x"; "sum_y"; "sum_w"])
    return DataFrame(out, colnames)
end

function fcn_run_optim(cost_df::DataFrame, kitl_index_full::DataFrame, stratified_samples::DataFrame, out_dir::AbstractString, sp::Integer, tt::Integer, kt::AbstractFloat, ns::Integer, baseline_conditions=false)
    Random.seed!(54815861) # Ensure all realisation samples are the same
    # Subset of properties
    run_string = "run_sp-$(sp)_tt-$(tt)_kt-$(kt)_ns-$(ns)_r-$(rr)"
    subset_id = vec(stratified_samples[:, sp])
    (cost_subset, kitl_subset) = subset_data(cost_df, kitl_index_full, subset_id)
    realisations = [realisation_sample(cost_subset, kitl_subset, tt, kt) for r in 1:R]
    realisation_areas = hcat([realisations[r].area for r in 1:R]...)
    realisation_areas_df = DataFrame(realisation_areas, :auto)
    CSV.write("$(out_dir)/area_$(run_string).csv", realisation_areas_df)
    worst_case_khab = minimum([minimum(sum(realisations[r].M₂, dims=1)) for r in 1:rr])
    metric_reshape = m -> DataFrame(reshape(permutedims(m, (1, 3, 2)), (size(m, 1) * size(m, 3), size(m, 2))), :auto)

    if (baseline_conditions)
        solution_nr = fcn_two_stage_opt_saa(realisations[1:rr]; K=7000, ns=ns, terminate_recourse=false, add_recourse=false, baseline_conditions=true)
        (cost_nr, metric_nr) = fcn_evaluate_solution(solution_nr.model, realisations)
        CSV.write("$(out_dir)/cost_baseline_$(run_string).csv", DataFrame(cost_nr, :auto))
        CSV.write("$(out_dir)/metric_baseline_$(run_string).csv", metric_reshape(metric_nr))
        out_nr = Solution(true, solution_nr.model, solution_nr.x, solution_nr.y, solution_nr.w, cost_nr, metric_nr)
        decision_nr = fcn_tidy_two_stage_solution_sum(solution_nr, cost_subset.NewPropID)
        CSV.write("$(out_dir)/decision_baseline_$(run_string).csv", decision_nr)
        @save "$(out_dir)/solution_baseline_$(run_string).jld" out_nr
        return out_nr
    end

    if (isfile("$(out_dir)/decision_fr_$(run_string).jld"))
        return
    end

    if (worst_case_khab < K)
        println("Constraint infeasible for this scenario")
        return
    else
        solution_nr = fcn_two_stage_opt_saa(realisations[1:rr]; K=7000, ns=ns, terminate_recourse=false, add_recourse=false)
        solution_ar = fcn_two_stage_opt_saa(realisations[1:rr]; K=7000, ns=ns, terminate_recourse=false, add_recourse=true)
        solution_tr = fcn_two_stage_opt_saa(realisations[1:rr]; K=7000, ns=ns, terminate_recourse=true, add_recourse=false)
        solution_fr = fcn_two_stage_opt_saa(realisations[1:rr]; K=7000, ns=ns, terminate_recourse=true, add_recourse=true)
    end

    (cost_nr, metric_nr) = fcn_evaluate_solution(solution_nr.model, realisations)
    (cost_ar, metric_ar) = fcn_evaluate_solution(solution_ar.model, realisations)
    (cost_tr, metric_tr) = fcn_evaluate_solution(solution_tr.model, realisations)
    (cost_fr, metric_fr) = fcn_evaluate_solution(solution_fr.model, realisations)

    decision_nr = fcn_tidy_two_stage_solution_sum(solution_nr, cost_subset.NewPropID)
    decision_ar = fcn_tidy_two_stage_solution_sum(solution_ar, cost_subset.NewPropID)
    decision_tr = fcn_tidy_two_stage_solution_sum(solution_tr, cost_subset.NewPropID)
    decision_fr = fcn_tidy_two_stage_solution_sum(solution_fr, cost_subset.NewPropID)

    median_ar = median((cost_ar .- cost_nr) ./ cost_nr)
    lb_ar = quantile(vec((cost_ar .- cost_nr) ./ cost_nr), 0.025)
    ub_ar = quantile(vec((cost_ar .- cost_nr) ./ cost_nr), 0.975)
    println("Value of recourse to add covenants ($(run_string)): $(round(median_ar*100))% [$(round(lb_ar*100))% - $(round(ub_ar*100))%]")

    median_fr = median((cost_fr .- cost_nr) ./ cost_nr)
    lb_fr = quantile(vec((cost_fr .- cost_nr) ./ cost_nr), 0.025)
    ub_fr = quantile(vec((cost_fr .- cost_nr) ./ cost_nr), 0.975)
    println("Value of full recourse ($(run_string)): $(round(median_fr*100))% [$(round(lb_fr*100))% - $(round(ub_fr*100))%]")

    # Save results
    CSV.write("$(out_dir)/cost_nr_$(run_string).csv", DataFrame(cost_nr, :auto))
    CSV.write("$(out_dir)/cost_ar_$(run_string).csv", DataFrame(cost_ar, :auto))
    CSV.write("$(out_dir)/cost_tr_$(run_string).csv", DataFrame(cost_tr, :auto))
    CSV.write("$(out_dir)/cost_fr_$(run_string).csv", DataFrame(cost_fr, :auto))

    CSV.write("$(out_dir)/metric_nr_$(run_string).csv", metric_reshape(metric_nr))
    CSV.write("$(out_dir)/metric_ar_$(run_string).csv", metric_reshape(metric_ar))
    CSV.write("$(out_dir)/metric_tr_$(run_string).csv", metric_reshape(metric_tr))
    CSV.write("$(out_dir)/metric_fr_$(run_string).csv", metric_reshape(metric_fr))

    CSV.write("$(out_dir)/decision_nr_$(run_string).csv", decision_nr)
    CSV.write("$(out_dir)/decision_ar_$(run_string).csv", decision_ar)
    CSV.write("$(out_dir)/decision_tr_$(run_string).csv", decision_tr)
    CSV.write("$(out_dir)/decision_fr_$(run_string).csv", decision_fr)

    out_nr = Solution(true, solution_nr.model, solution_nr.x, solution_nr.y, solution_nr.w, cost_nr, metric_nr)
    out_ar = Solution(true, solution_ar.model, solution_ar.x, solution_ar.y, solution_ar.w, cost_ar, metric_ar)
    out_tr = Solution(true, solution_tr.model, solution_tr.x, solution_tr.y, solution_tr.w, cost_tr, metric_tr)
    out_fr = Solution(true, solution_fr.model, solution_fr.x, solution_fr.y, solution_fr.w, cost_fr, metric_fr)
    #@save "$(out_dir)/solution_$(run_string).jld" out_nr out_ar out_tr out_fr
end

## Start sensitivity analysis

# Sensitivity parameters (base)
sp = 1; # 1:100
sp_vec = 1:50
tt = 6; # [0,1,2,3,4,5,6,7]
tt_vec = 3:6
kt = 0.25; # [0.1, 0.15, 0.2, 0.25, 0.3]
kt_vec = [0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3]
ns = 12; # 1:12
ns_vec = 1:12

baseline = fcn_run_optim(cost_df, kitl_index_full, stratified_samples, "results/mc_sim", sp, tt, kt, ns, true)

robust = fcn_run_optim(cost_df, kitl_index_full, stratified_samples, "results/mc_sim", sp, tt, kt, ns, false)

# Run across all ns, number of scenarios in resolved uncertainty
map((ns_i) -> fcn_run_optim(cost_df, kitl_index_full, stratified_samples, "results/mc_sim", sp, tt, kt, ns_i), ns_vec)

# Run across different time periods for tt
map((tt_i) -> fcn_run_optim(cost_df, kitl_index_full, stratified_samples, "results/mc_sim", sp, tt_i, kt, 1), tt_vec)
map((tt_i) -> fcn_run_optim(cost_df, kitl_index_full, stratified_samples, "results/mc_sim", sp, tt_i, kt, ns), tt_vec)

# Run across all kt
map((kt_i) -> fcn_run_optim(cost_df, kitl_index_full, stratified_samples, "results/mc_sim", sp, tt, kt_i, 1), kt_vec)
map((kt_i) -> fcn_run_optim(cost_df, kitl_index_full, stratified_samples, "results/mc_sim", sp, tt, kt_i, ns), kt_vec)

# Run across different sampled properties
map((sp_i) -> fcn_run_optim(cost_df, kitl_index_full, stratified_samples, "results/mc_sim", sp_i, tt, kt, 1), sp_vec)
map((sp_i) -> fcn_run_optim(cost_df, kitl_index_full, stratified_samples, "results/mc_sim", sp_i, tt, kt, ns), sp_vec)