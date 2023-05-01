using JLD2
using DataFrames
using CSV

include("optim-functions.jl")

mutable struct Solution
  feasible::Bool
  model::Model
  x::AbstractArray
  y::AbstractArray
  w::AbstractArray
  cost::AbstractArray
  metric::AbstractArray
end


function extract_jld_propid(jldpath::AbstractString, sp::Integer, stratified_samples::DataFrame)
  solution = load(jldpath)
  subset_id = vec(stratified_samples[:, sp])
  function fcn_tidy_two_stage_solution_sum(solution, prop_id)
    out = hcat(prop_id, solution.x, sum(solution.y, dims=2), sum(solution.w, dims=2))
    colnames = vcat(["NewPropID"; "x"; "sum_y"; "sum_w"])
    return DataFrame(out, colnames)
  end
  tidy_solution = fcn_tidy_two_stage_solution_sum(solution, subset_id)
  return (tidy_solution)
end

path = "/Users/frankiecho/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Documents/GitHub/koala-uncertainty/results/mc_sim"

stratified_samples = CSV.read("/Users/frankiecho/Library/CloudStorage/OneDrive-TheUniversityofQueensland/Documents/GitHub/koala-uncertainty/data/stratified_sample.csv", DataFrame);
jld_files = filter(x -> endswith(x, "r-10.jld"), readdir(path))
tidy_files = map((x) -> CSV.write(joinpath(path, split(x, ".")[1] * ".csv"), extract_jld_propid(joinpath(path, x), 1, stratified_samples)), jld_files)