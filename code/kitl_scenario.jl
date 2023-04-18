using Distributions
using StochasticPrograms
using CSV
using DataFrames
using Gurobi

include("optim-functions.jl")


adoption_subset = CSV.read("data/adoption_subset.csv", DataFrame)

function cost_to_ts(cost, area, adopt, idx_before=1:30, idx_after=31:60, inflation_rate = 0.02)
    delta = (1+inflation_rate).^(cat(idx_before, idx_after, dims=1));
    delta_before = sum(delta[idx_before]);
    delta_after = sum(delta[idx_after]);
    cost_ts = (cost .* 1000 .* area .* adopt ./ 10); # Cost per-year
    cost_before = cost_ts .* delta_before; # Cost start at 2020
    cost_after = cost_ts .* delta_after;
    return(cost_before, cost_after)
end

@define_scenario BehaviourScenario = begin
    C₁::Array{Float64,1}
    C₂::Array{Float64,1}
    M₁::Array{Float64,3}
    M₂::Array{Float64,3}
    @zero begin
        return BehaviourScenario(0.0, 0.0, 0.0, 0.0)
    end
end

@sampler BehaviourSampler = begin
    Adopt::Array{Float64,2}
    Prop::Array{Float64,2}
    WTA::Array{Float64,2}
    Area::Array{Float64,1}
    M₁::Array{Float64,3}
    M₂::Array{Float64,3}

    BehaviourSampler(Adopt, Prop, WTA, Area, M₁, M₂) = new(Adopt, Prop, WTA, Area, M₁, M₂)

    @sample BehaviourScenario begin
        MeanAdopt = sampler.Adopt[:,1]
        SDAdopt = sampler.Adopt[:,2]
        MeanProp = sampler.Prop[:,1]
        SDProp = sampler.Prop[:,2]
        MeanWTA = sampler.WTA[:,1]
        SDWTA = sampler.WTA[:,2]
        adopt_binary = (rand.(Normal.(MeanAdopt,SDAdopt))) .> rand(length(MeanAdopt))
        (cost_before, cost_after) = cost_to_ts(MeanWTA + SDWTA .* randn(length(SDWTA)), sampler.Area, MeanProp + SDProp .* randn(size(SDProp)))
        sample_C₁ = cost_before .* adopt_binary
        sample_C₂ = cost_after .* adopt_binary
        sample_M₁ = sampler.M₁ .* adopt_binary
        sample_M₂ = sampler.M₂ .* adopt_binary
        return BehaviourScenario(sample_C₁, sample_C₂, sample_M₁, sample_M₂)
    end
end



metric_subset = CSV.read("data/metric_subset.csv", DataFrame)
tt = 4;
t_symbols = [:t0, :t1, :t2, :t3, :t4, :t5, :t6, :t7];
N = length(unique(metric_subset.NewPropID))
climateProjList = ["CCCMA_R1", "CCCMA_R2", "CCCMA_R3", "CSIRO_R1", "CSIRO_R2", "CSIRO_R3", "ECHAM_R1", "ECHAM_R2", "ECHAM_R3", "MIROC_R1", "MIROC_R2", "MIROC_R3"];
S = length(climateProjList)
M₁ = zeros(N, length(1:(tt-1)), S)
M₂ = zeros(N, length(tt:length(t_symbols)), S)
for t=1:length(t_symbols)
  for s=1:S
    metric_unstacked = unstack(metric_subset, :NewPropID, :climate_model, t_symbols[t])
    if (t < tt)
      M₁[:,t,:] = Matrix(metric_unstacked[:, climateProjList])
    else
      M₂[:,t-tt+1,:] = Matrix(metric_unstacked[:, climateProjList])
    end
  end
end

sampler = BehaviourSampler(Matrix(adoption_subset[:,[:MeanAdopt, :SDAdopt]]), Matrix(adoption_subset[:,[:MeanProp, :SDProp]]), Matrix(adoption_subset[:,[:MeanWTA, :SDWTA]]), Vector(adoption_subset.:AREA), M₁, M₂);

add_recourse = true;
terminate_recourse = false;
p = ones(S,1)/S;
β = 0;
γ = 0;
K = 7000;

sm = @stochastic_model begin
    @stage 1 begin
        @decision(model, 0 <= x[1:N] <= 1) # First stage decision
        @variable(model, y[1:N, 1:S]) # Second stage decision: adding units
        @variable(model, w[1:N, 1:S]) # Second stage decision: terminating units
        
        if (add_recourse)
            @constraint(model, 0 <= y[1:N, 1:S] <= 1) # Second stage decision: adding units    
        else
            for i = 1:N
                for s = 1:S
                    fix(y[i,s], 0);
                end
            end
        end

        if (terminate_recourse)
            @constraint(model, 0 <= w[1:N, 1:S] <= 1) # Second stage decision: terminating units
        else
            for i = 1:N
                for s = 1:S
                    fix(w[i,s], 0);
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
    end
    @stage 2 begin
        @variable(model, z)
        @uncertain C₁ C₂ M₁ M₂ from BehaviourScenario
        @objective(model, Min, z)

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
end

function ExpectationRealisation(Adopt, Prop, WTA, Area, M₁, M₂)
    MeanAdopt = Adopt[:,1]
    SDAdopt = Adopt[:,2]
    MeanProp = Prop[:,1]
    SDProp = Prop[:,2]
    MeanWTA = WTA[:,1]
    SDWTA = WTA[:,2]
    adopt_binary = MeanAdopt
    (cost_before, cost_after) = cost_to_ts(MeanWTA, Area, MeanProp)
    sample_C₁ = cost_before .* adopt_binary
    sample_C₂ = cost_after .* adopt_binary
    sample_M₁ = M₁ .* adopt_binary
    sample_M₂ = M₂ .* adopt_binary
    return Realisation(sample_C₁, sample_C₂, sample_M₁, sample_M₂)
end

r = ExpectationRealisation(Matrix(adoption_subset[:,[:MeanAdopt, :SDAdopt]]), Matrix(adoption_subset[:,[:MeanProp, :SDProp]]), Matrix(adoption_subset[:,[:MeanWTA, :SDWTA]]), Vector(adoption_subset.:AREA), M₁, M₂)

model = fcn_two_stage_opt_robust(r);

#sp = instantiate(sm, sampler, 5, optimizer = LShaped.Optimizer)
#set_optimizer(sm, LShaped.Optimizer)
#set_optimizer_attribute(sp, MasterOptimizer(), Gurobi.Optimizer)
#set_optimizer_attribute(sp, SubProblemOptimizer(), Gurobi.Optimizer)

#optimize!(sp)