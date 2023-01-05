using Distributions
using Random
using Revise;
includet("optim-functions.jl")

Random.seed!(1234);
N = 17276; # planning units
S = 36; # uncertain scenarios
years = [2025,2035,2045,2055,2065,2075,2085];
T = length(years); # timesteps
δ = (1+0.02) .^ (years .- 2022);
tt = 4; # Timestep when uncertainty is revealed

A = rand(N,1); # Area

C = repeat(A.*5 + 0.5.*rand(Normal(0,0.1), N), outer= [1,T]);
C = C .* repeat(δ', outer=[N,1]);
C[C .< 0] .= 0;

function simulate_ar(phi, n)
    dist = Normal(0,0.1)
    y = [0.0 for i = 1:n]
    noise = rand(dist, n)
    for i in 1:(n-1)
        y[i+1] = phi*y[i] + noise[i] 
    end
    return y
end

m = rand(N); # Initial metric
M = [map((i->i.+simulate_ar(0.5,T)),m) |> (i->reduce(hcat,i)') for s in 1:S ];
for s=1:S
    M[s][M[s] .< 0] .= 0.0;
    M[s][M[s] .> 1] .= 1.0;
end
Mp = zeros(N, T, S);
for s=1:S
Mp[:,:,s] = M[s];
end
M = Mp;

p = rand(S);
p = p./sum(p);

β = 0.5;
γ = 0.5;
K = 50.0;

solution = fcn_two_stage_opt(C, M, K, p, 3);