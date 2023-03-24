using CSV
using DataFrames
using Revise
using StatsBase
using StochasticPrograms
using Gurobi
using JuMP
cd("/Users/frankiecho/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Documents/GitHub/koala-uncertainty/")

cost_df = CSV.read("data/meanWTA_10yr.csv", DataFrame);
habitat_df = CSV.read("data/habitat_suitw_graham.csv", DataFrame);

cost = Matrix(cost_df[:, [:cost2025, :cost2035, :cost2045, :cost2055, :cost2065, :cost2075, :cost2085]]);
N = size(cost,1);
T = size(cost,2);
S = length(unique(habitat_df.climate_model));

M = zeros(N, T, S);
for (s, cm) in enumerate(unique(habitat_df.climate_model))
    df = filter(row -> row.climate_model == cm, habitat_df);
    select!(df, Not([:climate_model, :NewPropID]))
    M[:,:,s] = Matrix(df);
end

M_end = M[:, 7, :];
ξ = [(@scenario m[i in 1:N]=Array(M_end[:, s]) probability = 1/S) for s=1:S];
two_stage_model = @stochastic_model begin
    @stage 1 begin
      @parameters begin
        K = 70.0
        tt = 4
        β = 0
        γ = 0
      end
      @decision(model, 0 <= x[i in 1:N] <= 1)
      @objective(model, Min, sum(cost[i,t]*x[i] for i in 1:N, t in 1:T))
      for t=1:(tt-1)
        for s=1:S
          @constraint(model, sum(M[i,t,s] * x[i] for i in 1:N) >= K)
        end
      end
    end
    @stage 2 begin
      @parameters begin
        K = 70.0
        tt = 4
        β = 0
        γ = 0
      end
      @uncertain m[i in 1:N]
      @recourse(model, 0 <= y[i in 1:N] <= 1)
      @recourse(model, 0 <= w[i in 1:N] <= 1)
      @objective(model, Min, sum((β + cost[i,t]) * y[i] + (γ - cost[i,t]) * w[i] for i in 1:N, t in tt:T))
      #for t=tt:T
      @constraint(model, sum(m[i] * (x[i] + y[i] - w[i]) for i in 1:N) >= K)
      #end
      for i=1:N
        @constraint(model, x[i] + y[i] <= 1)
        @constraint(model, w[i] <= x[i])
      end
    end
end

sp = instantiate(two_stage_model, ξ, optimizer = Gurobi.Optimizer)

set_time_limit_sec(sp, 10.0)
unset_silent(sp)
optimize!(sp)