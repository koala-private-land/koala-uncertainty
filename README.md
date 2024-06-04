# Flexible Koala Conservation under Climate Uncertainty
Working repository containing code for analysing koala conservation decision-making under uncertainty.

## Replication instructions

### Requirements
* Julia (v1.10.3)
* Microsoft VSCode (recommended)
* data directory (available upon request)
* Gurobi optimiser

### Instructions
1. Download Julia and install it using the Windows installer ([link](https://julialang-s3.julialang.org/bin/winnt/x64/1.10/julia-1.10.3-win64.exe))
2. Download VS Code ([link](https://code.visualstudio.com/Download))
3. Ensure that your system has a valid license for Gurobi and activate the link between Gurobi Julia package and the Gurobi installation following instructions here: [link](https://github.com/jump-dev/Gurobi.jl)
4. Open a Julia terminal, cd to the repo directory and do `] activate`
5. Create a directory `results/model_runs/` directory to collect results
6. Run `sensitivity_analysis.jl` to get the model run outputs
7. Open `plot_mc_results.R` in R/ RStudio, install required packages (using `install.packages()`) and then run the whole script. This should generate the image files that replicates the paper results.