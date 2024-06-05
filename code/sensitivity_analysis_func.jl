using CSV
using DataFrames
using Revise
using StatsBase
using Gurobi
using JuMP
using Distributions
using JLD2
using Random
using Distributed

println("Loading data files...")

cost_df = CSV.read("data/spatial_predictions_10yr.csv", DataFrame);
kitl_index_full = CSV.read("data/kitl_prop_climate.csv", DataFrame);
stratified_samples = CSV.read("data/stratified_sample.csv", DataFrame);
climateProjList = ["CCCMA_R1", "CCCMA_R2", "CCCMA_R3", "CSIRO_R1", "CSIRO_R2", "CSIRO_R3", "ECHAM_R1", "ECHAM_R2", "ECHAM_R3", "MIROC_R1", "MIROC_R2", "MIROC_R3"];
model_matrix_x = [Matrix(CSV.read("data/model_matrix/Model_MatrixX.Inf_$i.csv", DataFrame, header=false)) for i in 1:10];
model_matrix_y = [Matrix(CSV.read("data/model_matrix/Model_MatrixY.Inf_$i.csv", DataFrame, header=false)) for i in 1:10];
model_matrix_z = [Matrix(CSV.read("data/model_matrix/Model_MatrixZ.Inf_$i.csv", DataFrame, header=false)) for i in 1:10];
coefs_x = [Matrix(CSV.read("data/model_matrix/CoefsX.Inf_$i.csv", DataFrame, header=false)) for i in 1:10];
coefs_y = [Matrix(CSV.read("data/model_matrix/CoefsY.Inf_$i.csv", DataFrame, header=false)) for i in 1:10];
coefs_z = [Matrix(CSV.read("data/model_matrix/CoefsZ.Inf_$i.csv", DataFrame, header=false)) for i in 1:10];

println("Defining functions...")

addprocs(16)
@everywhere include("optim-functions.jl")
include("discounting.jl")

# Cost to time series
function cost_to_ts(cost::AbstractVector, area::AbstractVector, prop::AbstractVector, idx_before::AbstractVector=1:30, idx_after::AbstractVector=31:60, discount_rate::AbstractFloat= 0.02, perpetuity::Bool=true)
    # cost: WTA in thousands per hectare per year
    # area: Area in hectares
    # idx_before: index of the costs before covenant modification
    # idx_afer: index of the costs after covenant modification
    # discount_rate: discounting of future costs

    delta = (1 / (1+ discount_rate)) .^ (cat(idx_before, idx_after, dims=1))
    delta_before = sum(delta[idx_before])
    delta_after = sum(delta[idx_after])
    covenant_area = area .* prop;
    cost_ts = (cost .* 1000 .* covenant_area) # Cost per-year: cost per ha * area * proportion

    if (discount_rate == -99.0)
        # Use declining social discount rate per Green Book guidance
        cost_before = cost_ts .* declining_discount_annuity(1)
        cost_after = cost_ts .* declining_discount_annuity(30)
    elseif (perpetuity)
        cost_before = cost_ts ./ discount_rate # Cost start at 2020, and is a annuity payment until the end of idx_before
        cost_after = (cost_ts .* delta[idx_after[1]]) ./ discount_rate # Perpetuity costs starting at idx_after
    else
        cost_before = cost_ts .* delta # Cost start at 2020, and is a annuity payment until the end of idx_before
        cost_after = cost_ts .* delta_after # Cost start at 2020, and is a annuity payment until the end of idx_before
    end
    return (cost_before, cost_after)
end

# Subset data
function subset_propid(df, subset_id)
    return df[findall(in(subset_id),df.NewPropID),:]
end

function subset_data(cost_full, kitl_index_full, subset_id::Vector{Int64})
    cost_subset = subset_propid(cost_full, subset_id)
    kitl_subset = subset_propid(kitl_index_full, subset_id)
    cost_subset = subset_propid(cost_subset, kitl_subset.NewPropID)
    kitl_subset = subset_propid(kitl_subset, cost_subset.NewPropID)
    return (cost_subset, kitl_subset)
end

function subset_realisation(realisation, subset_id::Vector{Int64})
    C₁ = realisation.C₁[subset_id,:]
    C₂ =  realisation.C₂[subset_id,:]
    M₁ =  realisation.M₁[subset_id,:]
    M₂ =  realisation.M₂[subset_id,:]
    area =  realisation.area[subset_id]
    deforestation_risk = realisation.deforestation_risk[subset_id]
    τ =  realisation.τ[subset_id]
    p = realisation.p[subset_id]
    new_realisation = Realisation(C₁, C₂, M₁, M₂, area, deforestation_risk, τ, p)
    return new_realisation
end

function get_propid(subset_id::Vector{Int64})
    (cost_subset, kitl_subset) = subset_data(cost_df, kitl_index_full, subset_id)
    return vec(cost_subset.NewPropID)
end

# Metrics to input
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

# Realisation from MCMC sampling draws, with cost_subset constructed from fcn_get_predictions
function fcn_realisation_sample_draws(cost_subset::DataFrame, kitl_subset::DataFrame, tt::Integer=4, kt::AbstractFloat=0.25, discount_rate::Real=0.02, deforestation_risk::Real=0.03)
    (M₁, M₂) = metric_to_input(kitl_subset, cost_subset.AREA, tt, kt)
    adopt_binary = cost_subset.Adopt .> rand(length(cost_subset.Adopt))
    prop_binary = cost_subset.Prop .* adopt_binary
    (C₁, C₂) = cost_to_ts(cost_subset.WTA, cost_subset.AREA, prop_binary, 1:(10*(tt-1)), (1+10*(tt-1)):60, discount_rate)
    sample_C₁ = C₁ .* prop_binary
    sample_C₂ = C₂ .* prop_binary
    sample_M₁ = M₁ .* prop_binary
    sample_M₂ = M₂ .* prop_binary
    area = vec(cost_subset.AREA) .* prop_binary
    p = ones(size(sample_M₂, 3)) / size(sample_M₂, 3)
    #println("Generated a new realisation...")
    return Realisation(sample_C₁, sample_C₂, sample_M₁, sample_M₂, area, deforestation_risk, p)
end

function fcn_get_predictions()
    # Get cost, prop and adoption predictions from model matrix and coefficient multiplication
    imputation_idx = rand(1:10)
    draw_idx = rand(1:10000)

    # Get model matrix and coefficients
    model_matrix_x_idx = model_matrix_x[imputation_idx]
    model_matrix_y_idx = model_matrix_y[imputation_idx]
    model_matrix_z_idx = model_matrix_z[imputation_idx]
    coefs_x_idx = coefs_x[imputation_idx][draw_idx,:]
    coefs_y_idx = coefs_y[imputation_idx][draw_idx,:]
    coefs_z_idx = coefs_z[imputation_idx][draw_idx,:]

    # Get predictions
    pred_x = model_matrix_x_idx * coefs_x_idx
    pred_y = model_matrix_y_idx * coefs_y_idx
    pred_z = model_matrix_z_idx * coefs_z_idx

    # Logit transform adopt and prop 
    pred_x = exp.(pred_x) ./ (1 .+ exp.(pred_x))
    pred_z = exp.(pred_z) ./ (1 .+ exp.(pred_z))

    # Censor WTA predictions less than zero
    pred_y[pred_y .< 0] .= 0

    # Collate into DataFrame
    df = DataFrame(NewPropID = cost_df.NewPropID, WTA = pred_y, Prop = pred_z, Adopt = pred_x, AREA = cost_df.AREA)
    # Output predictions in a dataframe
    return df
end

function fcn_subset_realisation_sample_draw(kitl_index_full::DataFrame, subset_id::AbstractVector, tt::Integer=4, kt::AbstractFloat=0.25, discount_rate::Real=0.02, deforestation_risk=0.03, cost_uncertainty=false)
    # First stage costs
    cost_df = fcn_get_predictions();
    (cost_subset, kitl_subset) = subset_data(cost_df, kitl_index_full, subset_id)
    realisation = fcn_realisation_sample_draws(cost_subset, kitl_subset, tt, kt, discount_rate, deforestation_risk)

    if (cost_uncertainty)
        # Second stage costs (uncertain)
        pred_second_stage = [fcn_get_predictions() for s=1:S]
        cost_second_stage = [subset_data(pred_second_stage[s], kitl_subset, subset_id)[1] for s=1:S];
        realisations_second_stage = [fcn_realisation_sample_draws(cost_second_stage[s], kitl_subset, tt, kt, discount_rate, deforestation_risk) for s=1:S]
        C₂ = [realisations_second_stage[s].C₂ for s=1:S]
        M₂ = [realisations_second_stage[s].M₂[:,:,s] for s=1:S]
        realisation.C₂ = hcat(C₂...)
        realisation.M₂ = cat(M₂..., dims = 3)
    end

    return realisation;
end

function fcn_optim_output_model()
end

function fcn_run_optim(kitl_index_full::DataFrame, stratified_samples::DataFrame, out_dir::AbstractString, sp::Integer, tt::Integer, 
    kt::AbstractFloat, ns::Integer, discount_rate::Real=0.02, deforestation_risk::Real = 0.1, 
    K::Real = 7000, baseline_conditions=false, second_stage_budget::Float64 = 1.0e7, K_pa_change::Float64=0.0, recourses = (true,true,true,false,false))
    Random.seed!(54815861) # Ensure all realisation samples are the same
    run_string = "run_k-$(K)_sp-$(sp)_tt-$(tt)_kt-$(kt)_ns-$(ns)_r-$(rr)_sdr-$(ismissing(discount_rate) ? "ddr" : discount_rate)_dr-$(deforestation_risk)_ssb-$(second_stage_budget)_kpac-$(K_pa_change)"
    if (baseline_conditions && isfile("$(out_dir)/decision_baseline_$(run_string).csv"))
        println("Run $(run_string) skipped because output is present")
        return
    end
    if (!baseline_conditions && isfile("$(out_dir)/decision_ar_$(run_string).csv"))
        println("Run $(run_string) skipped because output is present")
        return
    end
    println("Starting run $(run_string)")
    subset_id = vec(stratified_samples[:, sp])
    realisations = [fcn_subset_realisation_sample_draw(kitl_index_full, subset_id, tt, kt, discount_rate, deforestation_risk) for r in 1:R]
    out_propid = get_propid(subset_id);
    realisation_areas = hcat([realisations[r].area for r in 1:R]...)
    realisation_areas_df = DataFrame(realisation_areas, :auto)
    CSV.write("$(out_dir)/area_$(run_string).csv", realisation_areas_df)
    realisation_cost_1 = hcat([realisations[r].C₁ for r in 1:R]...)
    realisation_cost_1_df = DataFrame(realisation_cost_1, :auto)
    CSV.write("$(out_dir)/cost1_$(run_string).csv", realisation_cost_1_df)

    # Write full costs of conservation covenant to table
    realisation_cost_full = hcat([realisations[r].C₁ + mean(realisations[r].C₂, dims = 2) for r in 1:R]...)
    realisation_cost_full_df = DataFrame(realisation_cost_full, :auto)
    CSV.write("$(out_dir)/cost_full_$(run_string).csv", realisation_cost_full_df)
    
    worst_case_khab = minimum([minimum(sum(realisations[r].M₂, dims=1)) for r in 1:rr])
    metric_reshape = m -> DataFrame(reshape(permutedims(m, (1, 3, 2)), (size(m, 1) * size(m, 3), size(m, 2))), :auto)
    (no_recourse, add_recourse, partial_recourse, terminate_recourse, full_recourse) = recourses

    if (baseline_conditions && no_recourse)
        models_nr_vec = pmap((realisation) -> fcn_two_stage_opt_deforestation(realisation; K=K, ns=ns, terminate_recourse=false, add_recourse=false, baseline_conditions=true), realisations)
        decision_nr = fcn_tidy_two_stage_solution_sum(models_nr_vec, out_propid)
        (cost_nr, metric_nr) = fcn_evaluate_solution(models_nr_vec, realisations)
        CSV.write("$(out_dir)/cost_baseline_$(run_string).csv", DataFrame(cost_nr, :auto))
        CSV.write("$(out_dir)/metric_baseline_$(run_string).csv", metric_reshape(metric_nr))
        CSV.write("$(out_dir)/decision_baseline_$(run_string).csv", decision_nr)
        #@save "$(out_dir)/solution_baseline_$(run_string).jld" out_nr
        return
    end

    if (worst_case_khab < K)
        println("Constraint infeasible for this scenario. Maximum achievable target is $(worst_case_khab)")
        return
    end

    if (no_recourse) 
        models_nr_vec = pmap((realisation) -> fcn_two_stage_opt_deforestation(realisation; K=K, ns=ns, terminate_recourse=false, add_recourse=false, K_pa_change=K_pa_change), realisations)
        (cost_nr, metric_nr) = fcn_evaluate_solution(models_nr_vec, realisations)
        decision_nr = fcn_tidy_two_stage_solution_sum(models_nr_vec, out_propid)
        CSV.write("$(out_dir)/cost_nr_$(run_string).csv", DataFrame(cost_nr, :auto))
        CSV.write("$(out_dir)/metric_nr_$(run_string).csv", metric_reshape(metric_nr))
        CSV.write("$(out_dir)/decision_nr_$(run_string).csv", decision_nr)
    end

    if (partial_recourse)
        models_pr_vec = pmap((realisation) -> fcn_two_stage_opt_deforestation(realisation; K=K, ns=ns, budget_constraint = (missing, second_stage_budget), terminate_recourse=false, add_recourse=true, K_pa_change=K_pa_change), realisations)
        (cost_pr, metric_pr) = fcn_evaluate_solution(models_pr_vec, realisations)
        decision_pr = fcn_tidy_two_stage_solution_sum(models_pr_vec, out_propid)
        CSV.write("$(out_dir)/cost_pr_$(run_string).csv", DataFrame(cost_pr, :auto))
        CSV.write("$(out_dir)/metric_pr_$(run_string).csv", metric_reshape(metric_pr))
        CSV.write("$(out_dir)/decision_pr_$(run_string).csv", decision_pr)
    end

    if (add_recourse) 
        models_ar_vec = pmap((realisation) -> fcn_two_stage_opt_deforestation(realisation; K=K, ns=ns, terminate_recourse=false, add_recourse=true, K_pa_change=K_pa_change), realisations)
        (cost_ar, metric_ar) = fcn_evaluate_solution(models_ar_vec, realisations)
        decision_ar = fcn_tidy_two_stage_solution_sum(models_ar_vec, out_propid)
        CSV.write("$(out_dir)/cost_ar_$(run_string).csv", DataFrame(cost_ar, :auto))
        CSV.write("$(out_dir)/metric_ar_$(run_string).csv", metric_reshape(metric_ar))
        CSV.write("$(out_dir)/decision_ar_$(run_string).csv", decision_ar)
        median_ar = median(skipmissing((cost_ar .- cost_pr) ./ cost_pr))
        println("Value of recourse to add covenants ($(run_string)): $(round(median_ar*100))%")
    end

    if (terminate_recourse) 
        models_tr_vec = pmap((realisation) -> fcn_two_stage_opt_deforestation(realisation; K=K, ns=ns, terminate_recourse=true, add_recourse=false, K_pa_change=K_pa_change), realisations)
        (cost_tr, metric_tr) = fcn_evaluate_solution(models_tr_vec, realisations)
        decision_tr = fcn_tidy_two_stage_solution_sum(models_tr_vec, out_propid)
        CSV.write("$(out_dir)/metric_tr_$(run_string).csv", metric_reshape(metric_tr))
        CSV.write("$(out_dir)/cost_tr_$(run_string).csv", DataFrame(cost_tr, :auto))
        CSV.write("$(out_dir)/decision_tr_$(run_string).csv", decision_tr)
    end

    if (full_recourse) 
        models_fr_vec = pmap((realisation) -> fcn_two_stage_opt_deforestation(realisation; K=K, ns=ns, terminate_recourse=true, add_recourse=false, K_pa_change=K_pa_change), realisations)
        (cost_fr, metric_fr) = fcn_evaluate_solution(solution_fr.model, realisations)
        decision_fr = fcn_tidy_two_stage_solution_sum(solution_fr, out_propid)
        median_fr = median((cost_fr .- cost_nr) ./ cost_nr)
        lb_fr = quantile(skipmissing(vec((cost_fr .- cost_nr) ./ cost_nr)), 0.025)
        ub_fr = quantile(skipmissing(vec((cost_fr .- cost_nr) ./ cost_nr)), 0.975)
        println("Value of full recourse ($(run_string)): $(round(median_fr*100))% [$(round(lb_fr*100))% - $(round(ub_fr*100))%]")
        CSV.write("$(out_dir)/cost_fr_$(run_string).csv", DataFrame(cost_fr, :auto))
        CSV.write("$(out_dir)/metric_fr_$(run_string).csv", metric_reshape(metric_fr))
        CSV.write("$(out_dir)/decision_fr_$(run_string).csv", decision_fr)
    end
end