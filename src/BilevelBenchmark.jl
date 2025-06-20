#module BilevelBenchmark

#using Plots
#using LaTeXStrings

using BilevelJuMP, Ipopt # External Dependencies from JSO

# Write your package code here.
#include("profiles.jl")
include("BOLIBProblems.jl")

get_bilevel_problem(2)

function BiJuMP_convertor(prob_no::Union{Int,String})
    prob = get_bilevel_problem(prob_no)
    nx, ny, nG, ng = prob.dim
    xy0 = prob.xy0
    f = prob.f_func
    F = prob.F_func
    g = prob.g_func
    G = prob.G_func

    model = BilevelModel(
    Ipopt.Optimizer,
    mode = BilevelJuMP.FortunyAmatMcCarlMode(primal_big_M = 100, dual_big_M = 100)
    )

    # Generating variables
    @variable(Upper(model), x[i = 1:nx], start = xy0[i])
    @variable(Lower(model), y[i = 1:ny], start = xy0[nx+i])

    # Generating objectives
    @objective(Upper(model), Min, F(x,y))
    @objective(Lower(model), Min, f(x,y))

    # Generating constraints
    for i in 1:nG
        @constraint(Upper(model), G(x,y)[i] ≤ 0)
    end
    for i in 1:ng
        @constraint(Upper(model), g(x,y)[i] ≤ 0)
    end
    return model
end

BiJuMP_model = BiJuMP_convertor(2)
optimize!(BiJuMP_model)

println(objective_value(BiJuMP_model))
println(value.(BiJuMP_model[:x]))
println(value.(BiJuMP_model[:y]))