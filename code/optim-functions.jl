# Solve a two-stage stochastic optimisation problem with JuMP where uncertainty is fully revealed in time τ
# Objective: minimise cost subject to a metric reaching a threshold across all time periods
using JuMP, Gurobi

Base.@kwdef mutable struct Realisation
    # C₁: first stage costs, each element of the vector is for each property unit (size: N, N: number of units)
    # C₂: second-stage costs, each element of the vector is for each property unit (length: N, N: number of units)
    # M₁: first-stage management objective indicator (size: N x 1)
    # M₂: second-stage management objective indicator, across all scenarios (N x S)
    # A: whether the landholder will consider a covenant or not (binary vector of length N)
    C₁::AbstractArray
    C₂::AbstractArray
    M₁::AbstractArray
    M₂::AbstractArray
    A::AbstractArray = ones(size(C₁, 1))
    p::AbstractArray = ones(size(M₂, 2)) / size(M₂, 2)

    function Realisation(C₁, C₂, M₁, M₂, A=ones(size(C₁, 1)), p=ones(size(M₂, 3)) / size(M₂, 3))
        return (new(C₁, C₂, M₁, M₂, A, p))
    end
end

function get_properties(r::Realisation)
    N = size(r.C₁, 1)
    S = size(r.M₂, 3)
    return ((N, S))
end

function fcn_two_stage_opt(C::Matrix, M::Array{Float64,3}, K::Real, p::Vector=(ones(size(M, 3)) / size(M, 3)), tt::Integer=Integer(round(size(M, 2) / 2)); β::Vector=zeros(size(M, 1)), γ::Vector=zeros(size(M, 1)), add_recourse=true, terminate_recourse=true)
    N, T, S = size(M) # N: number of units, T: number of timesteps, S: number of scenarios
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    @variable(model, 0 <= x[1:N] <= 1) # First stage decision
    @variable(model, 0 <= y[1:N, 1:S] <= (add_recourse ? 1 : 0)) # Second stage decision: adding units
    @variable(model, 0 <= w[1:N, 1:S] <= (terminate_recourse ? 1 : 0)) # Second stage decision: terminating units
    @objective(model, Min, sum(C[:, 1:T]' * x) + sum(p[s] * sum((β[i] + C[i, t]) * y[i, s] + (γ[i] - C[i, t]) * w[i, s] for i in 1:N, t in tt:T) for s in 1:S))
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

function fcn_two_stage_opt_saa(realisations=Vector{Realisation}; K::Real=7000, β::Real=0, γ::Real=0, add_recourse=true, terminate_recourse=true)
    # Solve the deterministic equivalent of the two-stage problem with sample average approximation

    # Arguments
    # realisations: a vector of Realisations
    # K: management objective

    N, S = get_properties(realisations[1]) # N: number of units, S: number of climate scenarios
    J = length(realisations) # J: number of realisations of adoption behaviour
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    @variable(model, 0 <= x[1:N] <= 1) # First stage decision
    @variable(model, 0 <= y[1:N, 1:S] <= (add_recourse ? 1 : 0)) # Second stage decision: adding units
    @variable(model, 0 <= w[1:N, 1:S] <= (terminate_recourse ? 1 : 0)) # Second stage decision: terminating units
    @variable(model, z) # Worst-case costs across all cost realisations
    @objective(model, Min, z)

    for j = 1:J
        C₁ = realisations[j].C₁ #.* realisations[j].A # First-stage cost under the J-th realisation
        C₂ = realisations[j].C₂ #.* realisations[j].A # Second-stage cost under the J-th realisation
        M₁ = realisations[j].M₁ #.* realisations[j].A
        M₂ = realisations[j].M₂ #.* realisations[j].A
        p = realisations[j].p  # Climate scenario probabilities

        # Objective is to minimise the worst-case (highest) costs, z
        @constraint(model, z >= sum((C₁ .+ C₂)' * x) + sum(p[s] * sum((β + C₂[i]) * y[i, s] + (γ - C₂[i]) * w[i, s] for i in 1:N) for s in 1:S))

        # Objective must be reached across all climate realisations before climate is known

        for t in axes(M₁, 2)
            for s in 1:S
                @constraint(model, M₁[:, t, s]' * x >= K)
            end
        end

        for t in axes(M₂, 2)
            for s = 1:S
                # After uncertainty is revealed, the objective only needs to be met at that realised climate realisation
                @constraint(model, M₂[:, t, s]' * (x + y[:, s] - w[:, s]) >= K)
            end
        end

    end

    for i = 1:N
        for s = 1:S
            @constraint(model, x[i] + y[i, s] <= 1)
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

function fcn_evaluate_solution(model::Model, realisations::Vector{Realisation})
    # model: Solved JuMP model
    # realisations: a vector of realisations

    x = value.(model[:x])
    y = value.(model[:y])
    w = value.(model[:w])

    J = length(realisations)
    N, S = get_properties(realisations[1])
    T = size(realisations[1].M₁, 2) + size(realisations[1].M₂, 2)

    cost_mat = zeros(S, J)
    metric_mat = zeros(S, T, J)

    for j = 1:J
        r = realisations[j]
        cost = [r.C₁' * x + r.C₂' * (x + y[:, s] - w[:, s]) for s = 1:S]
        metric_1 = mapreduce(permutedims, vcat, [r.M₁[:, :, s]' * x for s = 1:S])'
        metric_2 = mapreduce(permutedims, vcat, [r.M₂[:, :, s]' * (x + y[:, s] - w[:, s]) for s = 1:S])'
        metric = vcat(metric_1, metric_2)

        cost_mat[:, j] = cost
        metric_mat[:, :, j] = metric'
    end

    return (cost_mat, metric_mat)
end