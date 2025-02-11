# Solve a two-stage stochastic optimisation problem with JuMP where uncertainty is fully revealed in time τ
# Objective: minimise cost subject to a metric reaching a threshold across all time periods

using JuMP, Gurobi, Random, Distributed

GRB_ENV = Gurobi.Env()

Base.@kwdef mutable struct Realisation
    # C₁: first stage costs, each element of the vector is for each property unit (size: N, N: number of units)
    # C₂: second-stage costs, each element of the vector is for each property unit (length: N, N: number of units)
    # M₁: first-stage management objective indicator (size: N x 1)
    # M₂: second-stage management objective indicator, across all scenarios (N x S)
    # Area: total land area available for conservation
    # A: whether the landholder will consider a covenant or not (binary vector of length N)
    # deforestation_risk: probability that land will be deforested between now and t'
    # τ: indicator of land free of land clearing even without covenants

    C₁::AbstractArray
    C₂::AbstractArray
    M₁::AbstractArray
    M₂::AbstractArray
    area::AbstractArray
    deforestation_risk::Float64
    τ::AbstractArray
    p::AbstractArray

    function Realisation(C₁::AbstractArray, C₂::AbstractArray, M₁::AbstractArray, M₂::AbstractArray, area::AbstractArray, deforestation_risk::Float64, p::AbstractArray)
        τ = (rand(size(C₁, 1), length(p)) .< (1-deforestation_risk)) # Deforestation for each scenario
        return (new(C₁, C₂, M₁, M₂, area, deforestation_risk, τ, p))
    end
end

mutable struct Solution
    model::Model
    obj_value::Real
    x::AbstractArray
    y::AbstractArray
    w::AbstractArray
end


function get_properties(r::Realisation)
    N = size(r.C₁, 1)
    S = size(r.M₂, 3)
    return ((N, S))
end

function fcn_two_stage_opt_deforestation(realisation::Realisation; K::Real=7000, add_recourse=true, terminate_recourse=true, ns::Integer=1, baseline_conditions::Bool=false, budget_constraint=(missing, missing), min_spend=(missing,missing), K_pa_change::Float64 = 0.0)::Solution
    # Solve the deterministic equivalent of the two-stage problem with sample average approximation
    # In contrast to the SAA approach, this approach assumes that the decision maker knows the bids of landholders for certain

    # Arguments
    # realisation: a single realisation
    # K: management objective
    # add_recourse: whether to allow recourse to adding units in the second stage
    # terminate_recourse: whether to allow recourse to terminating units in the second stage
    # ns: number of scenarios in each resolution of uncertainty
    # baseline_conditions: whether or not planning is done assuming baseline conditions
    # budget_constraint: budget constraint (maximum) to the first-stage and second-stage spending, with None being unconstrained
    # min_spend: minimum amount of spending in the first-stage and second-stage
    # K_pa_change: per-annum change in the target
    println("Starting optimisation iteration...")
    N, S = get_properties(realisation) # N: number of units, S: number of climate scenarios
    s_vec = fcn_resolve_scenario_incomplete(ns)

    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_silent(model)
    set_attribute(model, "OutputFlag", 0)
    @variable(model, 0 <= x[1:N] <= 1) # First stage decision
    @variable(model, z) # Worst-case costs across all cost realisations
    @objective(model, Min, z)

    # Destructure realisation
    (; C₁, C₂, M₁, M₂, p, τ) = realisation

    
    if (add_recourse)
        @variable(model, 0 <= y[i = 1:N, s = 1:S] <= τ[i,s]) # Second stage decision: adding units, cleared land cannot be protected in the second stage
        @constraint(model, [i in 1:N, s in 1:S], x[i] + y[i, s] <= 1)
    end
    
    if (terminate_recourse)
        @variable(model, 0 <= w[1:N, 1:S] <= 1) # Second stage decision: terminating units
        @constraint(model, [i in 1:N, s in 1:S], x[i] - w[i, s] >= 0)
    end
    
    if (baseline_conditions) # Assume conditions persist according to the baseline
        # Conservation goals only need to met on average across all timesteps
        for t in axes(M₂, 2)
            @constraint(model, sum(p[s] * M₂[:, t, s]' * x for s in 1:S) .>= K)
        end
    else
        # Objective must be reached across all climate realisations before climate is known
        @constraint(model, [t in axes(M₁, 2), s in 1:S], M₁[:, t, s]' * x .>= K * (1+K_pa_change)^((t-1)*10))
        
        # After uncertainty is revealed, the objective only needs to be met at that scenario
        @constraint(model, [t in axes(M₂, 2), s in 1:S], M₂[:, t, s_vec[s]]' * (x .+ (add_recourse ? y[:, s] : 0) .- (terminate_recourse ? w[:, s] : 0)) .>= K * (1+K_pa_change)^((t+ size(M₁, 2) - 1)*10-1))
    end
    
    # Objective is to minimise costs, z
    @constraint(model, z >= sum((C₁)' * x) + sum(p[s] * sum((add_recourse ? C₂[i] * y[i, s] : 0) + (terminate_recourse ? C₂[i] * w[i, s] : 0) for i in 1:N) for s in 1:S))

    # Constrain maximum spending in Stage 1 and 2
    if (~ismissing(budget_constraint[1]))
        @constraint(model, sum((C₁)' * x) <= budget_constraint[1])
    end

    if (~ismissing(budget_constraint[2]))
        for s in 1:S
            @constraint(model, sum((add_recourse ? C₂[i] * y[i, s] : 0) + (terminate_recourse ? C₂[i] * w[i, s] : 0) for i in 1:N) <= budget_constraint[2])
        end
    end

    # Constrain minimum spending in Stage 1 and 2
    if (~ismissing(min_spend[1]))
        @constraint(model, sum((C₁)' * x) >= min_spend[1])
    end

    if (~ismissing(min_spend[2]))
        for s in 1:S
            @constraint(model, sum((add_recourse ? C₂[i] * y[i, s] : 0) + (terminate_recourse ? C₂[i] * w[i, s] : 0) for i in 1:N) <= min_spend[2])
        end
    end

    optimize!(model)
    out = Solution(
        model,
        objective_value(model),
        value.(x),
        (add_recourse ? value.(y) : zeros(N, S)),
        (terminate_recourse ? value.(w) : zeros(N, S))
    )
    return out
end

function fcn_evaluate_solution(solutions::Vector{Solution}, realisations::Vector{Realisation})
    # solutions: a vector of solved optimisation solutions
    # realisations: a vector of realisations, each realisation corresponds to each model by index

    if (length(solutions) != length(realisations))
        error("The number of models and realisations must be equal")
    end

    J = length(realisations)
    
    N, S = get_properties(realisations[1])
    T = size(realisations[1].M₁, 2) + size(realisations[1].M₂, 2)
    
    cost_mat = zeros(S, J)
    metric_mat = zeros(S, T, J)

    function fcn_evaluate_solution_model(solution::Solution, r::Realisation)
        x = solution.x
        y = solution.y
        w = solution.w

        cost = [sum(r.C₁'*x + r.C₂'*(x + y[:, s] - w[:, s])) for s = 1:S]
        metric_1 = mapreduce(permutedims, vcat, [r.M₁[:, :, s]' * x for s = 1:S])'
        metric_2 = mapreduce(permutedims, vcat, [r.M₂[:, :, s]' * (x + y[:, s] - w[:, s]) for s = 1:S])'
        metric = vcat(metric_1, metric_2)
        return (cost, metric)
    end

    for j=1:J
        (cost, metric) = fcn_evaluate_solution_model(solutions[j], realisations[j])
        cost_mat[:, j] = cost
        metric_mat[:, :, j] = metric'
    end

    return (cost_mat, metric_mat)
end


function fcn_resolve_scenario_incomplete(ns::Integer=2; S::Integer=12, N::Integer=12)
    # S: Number of possible resolutions of uncertainty
    # N: Number of scenarios
    # ns: Number of scenarios in each resolution of uncertainty
    if (ns > N)
        error("ns must be less than or equal to N")
    end
    if (ns == 1 && N <= S)
        return (1:N)
    else
        generate_scen = () -> vcat(shuffle([randperm(3) .+ (a - 1) * 3 for a in 1:4])...)
        return ([sort(generate_scen()[1:ns]) for s in 1:S])
    end
end

function fcn_tidy_two_stage_solution_sum(solution::Solution, prop_id)
    out = hcat(prop_id, solution.x, sum(solution.y, dims=2), sum(solution.w, dims=2))
    colnames = vcat(["NewPropID"; "x"; "sum_y"; "sum_w"])
    return DataFrame(out, colnames)
end

function fcn_tidy_two_stage_solution_sum(solutions::Vector{Solution}, prop_id)
    N = length(value.(solutions[1].x));
    x = mapreduce(permutedims, vcat, [soln.x for soln in solutions])
    y = mapreduce(permutedims, vcat, [sum(soln.y, dims = 2) for soln in solutions])
    w = mapreduce(permutedims, vcat, [sum(soln.w, dims = 2) for soln in solutions])

    out = hcat(prop_id, x', y', w')
    colnames = vcat(["NewPropID"; ["x$(i)" for i in range(1,length(solutions))]; ["sum_y$(i)" for i in range(1,length(solutions))]; ["sum_w$(i)" for i in range(1,length(solutions))]])
    return DataFrame(out, colnames)
end