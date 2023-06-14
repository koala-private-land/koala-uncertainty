ENV["GUROBI_HOME"] = "m:\\Users\\uqfcho\\Documents\\gurobi1001\\win64"
using Pkg
Pkg.add("Gurobi")
Pkg.build("Gurobi")