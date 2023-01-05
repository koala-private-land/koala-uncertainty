# Solve a two-stage stochastic optimisation problem with JuMP where uncertainty is fully revealed in time τ
# Objective: minimise cost subject to a metric reaching a threshold across all time periods
using JuMP, Gurobi

function fcn_two_stage_opt(C::Matrix, M::Array{Float64,3}, p::Vector, tt::Integer = Integer(round(size(M,2)/2)); β=0, γ=0, add_recourse = true, terminate_recourse = true)
    N, T, S = size(M); # N: number of units, T: number of timesteps, S: number of scenarios
    model = Model(Gurobi.Optimizer);
    set_silent(model);
    @variable(model, 0 <= x[1:N] <= 1); # First stage decision
    @variable(model, 0 <= y[1:N, 1:S] <= (add_recourse ? 1 : 0)); # Second stage decision: adding units
    @variable(model, 0 <= w[1:N, 1:S] <= (terminate_recourse ? 1 : 0)); # Second stage decision: terminating units
    @objective(model, Min, sum(C[:, 1:T]'*x) + sum(p[s] * sum((β+C[i,t])*y[i,s] + (γ-C[i,t])*w[i,s] for i in 1:N, t in tt:T) for s in 1:S));
    for t=1:(tt-1)
        for s=1:S
            @constraint(model, sum(M[:,t,s]' * x) >= K);
        end
    end

    for t=tt:T
        for s=1:S
            @constraint(model, sum(M[:,t,s]' * (x+y[:,s]-w[:,s])) >= K); 
        end
    end

    for i=1:N
        for s=1:S
            @constraint(model, x[i] + y[i,s] <= 1);
        end
    end

    for i=1:N
        for s=1:S
            @constraint(model, x[i] - w[i,s] >= 0);
        end
    end

    optimize!(model);

    return (
        model = model,
        obj_value = objective_value(model),
        x = value.(x),
        y = value.(y),
        w = value.(w)
    );
end