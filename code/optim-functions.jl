# Solve a two-stage stochastic optimisation problem with JuMP where uncertainty is fully revealed in time τ
# Objective: minimise cost subject to a metric reaching a threshold across all time periods

using JuMP, Gurobi, Random

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
    feasible::Bool
    model::Model
    x::AbstractArray
    y::AbstractArray
    w::AbstractArray
    cost::AbstractArray
    metric::AbstractArray
end


function get_properties(r::Realisation)
    N = size(r.C₁, 1)
    S = size(r.M₂, 3)
    return ((N, S))
end

function fcn_two_stage_opt(C::Matrix, M::Array{Float64,3}, K::Real, p::Vector=(ones(size(M, 3)) / size(M, 3)), tt::Integer=Integer(round(size(M, 2) / 2)); β::Vector=zeros(size(M, 1)), γ::Vector=zeros(size(M, 1)), add_recourse=true, terminate_recourse=true, baseline_conditions=false)
    N, T, S = size(M) # N: number of units, T: number of timesteps, S: number of scenarios
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    @variable(model, 0 <= x[1:N] <= 1) # First stage decision
    @variable(model, 0 <= y[1:N, 1:S] <= (add_recourse ? 1 : 0)) # Second stage decision: adding units
    @variable(model, 0 <= w[1:N, 1:S] <= (terminate_recourse ? 1 : 0)) # Second stage decision: terminating units
    @objective(model, Min, sum(C[:, 1:T]' * x) + sum(p[s] * sum((β[i] + C[i, t]) * y[i, s] + (γ[i] - C[i, t]) * w[i, s] for i in 1:N, t in tt:T) for s in 1:S))

    if (baseline_conditions)
        # Require only baseline conditions be met
        for s = 1:S
            @constraint(model, sum(M[:, 1, s]' * x) >= K)
        end
    else
        for t = 1:(tt-1)
            for s = 1:S
                @constraint(model, sum(M[:, t, s]' * x) >= K)
            end
        end

        for t = tt:T
            for s = 1:S
                @constraint(model, sum(M[:, t, s]' * (x + y[:, s] - w[:, s])) >= K)
            end
        end

        for i = 1:N
            for s = 1:S
                @constraint(model, x[i] + y[i, s] <= 1)
            end
        end
    end

    for i = 1:N
        for s = 1:S
            @constraint(model, x[i] - w[i, s] >= 0)
        end
    end

    optimize!(model)

    return (
        model=model,
        obj_value=objective_value(model),
        x=value.(x),
        y=value.(y),
        w=value.(w)
    )
end

# function fcn_two_stage_opt_saa(realisations=Vector{Realisation}; K::Real=7000, β::Real=0, γ::Real=0, add_recourse=true, terminate_recourse=true)
#     # Solve the deterministic equivalent of the two-stage problem with sample average approximation

#     # Arguments
#     # realisations: a vector of Realisations
#     # K: management objective

#     N, S = get_properties(realisations[1]) # N: number of units, S: number of climate scenarios
#     J = length(realisations) # J: number of realisations of adoption behaviour
#     model = Model(Gurobi.Optimizer)
#     #set_silent(model)
#     @variable(model, 0 <= x[1:N] <= 1) # First stage decision
#     @variable(model, 0 <= y[1:N, 1:S] <= (add_recourse ? 1 : 0)) # Second stage decision: adding units
#     @variable(model, 0 <= w[1:N, 1:S] <= (terminate_recourse ? 1 : 0)) # Second stage decision: terminating units
#     @variable(model, z) # Worst-case costs across all cost realisations
#     @objective(model, Min, z)

#     for j = 1:J
#         C₁ = realisations[j].C₁ #.* realisations[j].A # First-stage cost under the J-th realisation
#         C₂ = realisations[j].C₂ #.* realisations[j].A # Second-stage cost under the J-th realisation
#         M₁ = realisations[j].M₁ #.* realisations[j].A
#         M₂ = realisations[j].M₂ #.* realisations[j].A
#         p = realisations[j].p  # Climate scenario probabilities

#         # Objective is to minimise the worst-case (highest) costs, z
#         @constraint(model, z >= sum((C₁ .+ C₂)' * x) + sum(p[s] * sum((β + C₂[i]) * y[i, s] + (γ - C₂[i]) * w[i, s] for i in 1:N) for s in 1:S))

#         # Objective must be reached across all climate realisations before climate is known

#         for t in axes(M₁, 2)
#             for s in 1:S
#                 @constraint(model, M₁[:, t, s]' * x >= K)
#             end
#         end

#         for t in axes(M₂, 2)
#             for s = 1:S
#                 # After uncertainty is revealed, the objective only needs to be met at that realised climate realisation
#                 @constraint(model, M₂[:, t, s]' * (x + y[:, s] - w[:, s]) >= K)
#             end
#         end

#     end

#     for i = 1:N
#         for s = 1:S
#             @constraint(model, x[i] + y[i, s] <= 1)
#         end
#     end

#     for i = 1:N
#         for s = 1:S
#             @constraint(model, x[i] - w[i, s] >= 0)
#         end
#     end

#     optimize!(model)

#     return (
#         model=model,
#         obj_value=objective_value(model),
#         x=value.(x),
#         y=value.(y),
#         w=value.(w)
#     )
# end

function fcn_two_stage_opt_saa(realisations=Vector{Realisation}; K::Real=7000, β::Real=0, γ::Real=0, add_recourse=true, terminate_recourse=true, ns::Integer=1, baseline_conditions=false)
    # Solve the deterministic equivalent of the two-stage problem with sample average approximation

    # Arguments
    # realisations: a vector of Realisations
    # K: management objective
    # β: cost of adding a unit in the second stage
    # γ: cost of terminating a unit in the second stage
    # add_recourse: whether to allow recourse to adding units in the second stage
    # terminate_recourse: whether to allow recourse to terminating units in the second stage
    # ns: number of scenarios in each resolution of uncertainty

    N, S = get_properties(realisations[1]) # N: number of units, S: number of climate scenarios
    J = length(realisations) # J: number of realisations of adoption behaviour
    s_vec = fcn_resolve_scenario_incomplete(ns)

    model = Model(Gurobi.Optimizer)
    set_silent(model)
    @variable(model, 0 <= x[1:N] <= 1) # First stage decision
    @variable(model, z) # Worst-case costs across all cost realisations
    @objective(model, Min, z)

    # Objective is to minimise expected costs, z
    C₁ = dropdims(mean(cat([realisations[j].C₁ for j in 1:J]..., dims=3), dims=3), dims=3) #.* realisations[j].A # First-stage cost under the J-th realisation
    C₂ = dropdims(mean(cat([realisations[j].C₂ for j in 1:J]..., dims=3), dims=3), dims=3) #.* realisations[j].A # Second-stage cost under the J-th realisation
    p = realisations[1].p  # Climate scenario probabilities

    if (add_recourse)
        @variable(model, 0 <= y[1:N, 1:S] <= 1) # Second stage decision: adding units
        @constraint(model, [i in 1:N, s in 1:S], x[i] + y[i, s] <= 1)
    end

    if (terminate_recourse)
        @variable(model, 0 <= w[1:N, 1:S] <= 1) # Second stage decision: terminating units
        @constraint(model, [i in 1:N, s in 1:S], x[i] - w[i, s] >= 0)
    end

    @constraint(model, z >= sum((C₁ .+ C₂)' * x) + sum(p[s] * sum((add_recourse ? (β + C₂[i]) * y[i, s] : 0) + (terminate_recourse ? (γ - C₂[i]) * w[i, s] : 0) for i in 1:N) for s in 1:S))

    for j = 1:J
        M₁ = realisations[j].M₁ #.* realisations[j].A
        M₂ = realisations[j].M₂ #.* realisations[j].A

        if (baseline_conditions)
            # Conservation goals only need to met under baseline conditions (i.e. third time step)
            @constraint(model, [s in 1:S], M₁[:, 3, s]' * x .>= K)
        else
            # Objective must be reached across all climate realisations before climate is known
            @constraint(model, [t in axes(M₁, 2), s in 1:S], M₁[:, t, s]' * x .>= K)

            # After uncertainty is revealed, the objective only needs to be met at that realised climate realisation
            @constraint(model, [t in axes(M₂, 2), s in 1:S], M₂[:, t, s_vec[s]]' * (x .+ (add_recourse ? y[:, s] : 0) .- (terminate_recourse ? w[:, s] : 0)) .>= K)
        end

    end

    optimize!(model)
    return (
        model=model,
        obj_value=objective_value(model),
        x=value.(x),
        y=(add_recourse ? value.(y) : zeros(N, S)),
        w=(terminate_recourse ? value.(w) : zeros(N, S))
    )
end

function fcn_two_stage_opt_deforestation(realisation::Realisation; K::Real=7000, add_recourse=true, terminate_recourse=true, ns::Integer=1, baseline_conditions=false)
    # Solve the deterministic equivalent of the two-stage problem with sample average approximation
    # In contrast to the SAA approach, this approach assumes that the decision maker knows the bids of landholders for certain

    # Arguments
    # realisation: a single realisation
    # K: management objective
    # add_recourse: whether to allow recourse to adding units in the second stage
    # terminate_recourse: whether to allow recourse to terminating units in the second stage
    # ns: number of scenarios in each resolution of uncertainty

    N, S = get_properties(realisation) # N: number of units, S: number of climate scenarios
    s_vec = fcn_resolve_scenario_incomplete(ns)

    model = Model(Gurobi.Optimizer)
    set_silent(model)
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
    
    if (baseline_conditions)
        # Conservation goals only need to met under the terminal timestep on average
        @constraint(model, sum(p[s] * M₂[:, axes(M₂, 2), s]' * x for s in 1:S) .>= K)
    else
        # Objective must be reached across all climate realisations before climate is known
        @constraint(model, [t in axes(M₁, 2), s in 1:S], M₁[:, t, s]' * x .>= K)
        
        # After uncertainty is revealed, the objective only needs to be met at that scenario
        @constraint(model, [t in axes(M₂, 2), s in 1:S], M₂[:, t, s_vec[s]]' * (x .+ (add_recourse ? y[:, s] : 0) .- (terminate_recourse ? w[:, s] : 0)) .>= K)
    end
    
    # Objective is to minimise costs, z
    @constraint(model, z >= sum((C₁)' * x) + sum(p[s] * sum((add_recourse ? C₂[i, s] * y[i, s] : 0) + (terminate_recourse ? C₂[i, s] * w[i, s] : 0) for i in 1:N) for s in 1:S))

    optimize!(model)
    return (
        model=model,
        obj_value=objective_value(model),
        x=value.(x),
        y=(add_recourse ? value.(y) : zeros(N, S)),
        w=(terminate_recourse ? value.(w) : zeros(N, S))
    )
end

function fcn_evaluate_solution(model::Model, realisations::Vector{Realisation})
    # model: Solved JuMP model
    # realisations: a vector of realisations

    J = length(realisations)
    N, S = get_properties(realisations[1])
    T = size(realisations[1].M₁, 2) + size(realisations[1].M₂, 2)

    x = value.(model[:x])
    y = haskey(model, :y) ? value.(model[:y]) : zeros(N, S)
    w = haskey(model, :w) ? value.(model[:w]) : zeros(N, S)

    cost_mat = zeros(S, J)
    metric_mat = zeros(S, T, J)

    for j = 1:J
        r = realisations[j]
        cost = [sum(r.C₁'*x + r.C₂[:,s]'*(x + realisations[s].τ .* y[:, s] - w[:, s])) for s = 1:S]
        metric_1 = mapreduce(permutedims, vcat, [r.M₁[:, :, s]' * x for s = 1:S])'
        metric_2 = mapreduce(permutedims, vcat, [r.M₂[:, :, s]' * (x + realisations[s].τ .* y[:, s] - w[:, s]) for s = 1:S])'
        metric = vcat(metric_1, metric_2)
        cost_mat[:, j] = cost
        metric_mat[:, :, j] = metric'
    end

    return (cost_mat, metric_mat)
end

function fcn_evaluate_solution(model::Vector{Model}, realisations::Vector{Realisation})
    # model: a vector of solved JuMP models
    # realisations: a vector of realisations, each realisation corresponds to each model by index

    if (length(model) != length(realisations))
        error("The number of models and realisations must be equal")
    end

    J = length(realisations)
    
    N, S = get_properties(realisations[1])
    T = size(realisations[1].M₁, 2) + size(realisations[1].M₂, 2)
    
    cost_mat = zeros(S, J)
    metric_mat = zeros(S, T, J)

    function fcn_evaluate_solution_model(model::Model, r::Realisation)
        x = value.(model[:x])
        y = haskey(model, :y) ? value.(model[:y]) : zeros(N, S)
        w = haskey(model, :w) ? value.(model[:w]) : zeros(N, S)

        cost = [sum(r.C₁'*x + r.C₂'*(x + y[:, s] - w[:, s])) for s = 1:S]
        metric_1 = mapreduce(permutedims, vcat, [r.M₁[:, :, s]' * x for s = 1:S])'
        metric_2 = mapreduce(permutedims, vcat, [r.M₂[:, :, s]' * (x + y[:, s] - w[:, s]) for s = 1:S])'
        metric = vcat(metric_1, metric_2)
        return (cost, metric)
    end

    for j=1:J
        (cost, metric) = fcn_evaluate_solution_model(model[j], realisations[j])
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

function fcn_tidy_two_stage_solution_sum(models::Vector{Model}, prop_id)
    N = length(value.(models[1][:x]));
    x = mapreduce(permutedims, vcat, [value.(model[:x]) for model in models])
    y = mapreduce(permutedims, vcat, [haskey(model, :y) ? sum(value.(model[:y]), dims=2) : zeros(N, 1) for model in models])
    w = mapreduce(permutedims, vcat, [haskey(model, :w) ? sum(value.(model[:w]), dims=2) : zeros(N, 1) for model in models])

    out = hcat(prop_id, x', y', w')
    colnames = vcat(["NewPropID"; ["x$(i)" for i in range(1,length(models))]; ["sum_y$(i)" for i in range(1,length(models))]; ["sum_w$(i)" for i in range(1,length(models))]])
    return DataFrame(out, colnames)
end