using StochasticPrograms
using Gurobi

@stochastic_model simple_model begin
    @stage 1 begin
        @decision(simple_model, x₁ >= 40)
        @decision(simple_model, x₂ >= 20)
        @objective(simple_model, Min, 100*x₁ + 150*x₂)
        @constraint(simple_model, x₁ + x₂ <= 120)
    end
    @stage 2 begin
        @uncertain q₁ q₂ d₁ d₂
        @recourse(simple_model, 0 <= y₁ <= d₁)
        @recourse(simple_model, 0 <= y₂ <= d₂)
        @objective(simple_model, Max, q₁*y₁ + q₂*y₂)
        @constraint(simple_model, 6*y₁ + 10*y₂ <= 60*x₁)
        @constraint(simple_model, 8*y₁ + 5*y₂ <= 80*x₂)
    end
end

@scenario ξ[i in 1:4] = [24.0, 28.0, 500.0, 100.0] probability = 0.4
sp = instantiate(simple_model, [ξ₁, ξ₂], optimizer = Gurobi.Optimizer)
print(sp)
optimize!(sp)