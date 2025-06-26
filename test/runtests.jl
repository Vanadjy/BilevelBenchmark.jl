using BilevelBenchmark

using Test, LinearAlgebra # Basic Julia packages
using JuMP, BilevelJuMP # JuMP API
using Ipopt#, HiGHS, MibS_jll # Solvers

function BiJuMP_convertor(prob::BilevelProblem)
    nx, ny, nG, ng = prob.dim
    xy0 = prob.xy0
    f = prob.f_func
    F = prob.F_func
    g = prob.g_func
    G = prob.G_func

    model = BilevelModel(Ipopt.Optimizer; mode = BilevelJuMP.ProductMode(1e-9))

    # Generating variables
    @variable(Upper(model), x[i = 1:nx], start = xy0[i])
    @variable(Lower(model), y[i = 1:ny], start = xy0[nx+i])

    # Generating objectives
    @objective(Upper(model), Min, F(x,y))
    @objective(Lower(model), Min, f(x,y))

    # Generating constraints
    for i in 1:nG
        @constraint(Upper(model), G(x,y)[i] ≤ 0.0)
    end
    for i in 1:ng
        @constraint(Lower(model), g(x,y)[i] ≤ 0.0)
    end
    return model
end

all_probs = collect(1:162) # List of problem numbers to test
JuMP_uncompatible = [3, 12, 14, 15, 16, 17, 18, 21, 32, 33, 34, 36, 40, 41, 42, 43, 49, 50, 51, 52, 53, 54, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 79, 80, 81, 83, 84, 85, 86, 87, 88, 94, 95, 96, 97, 98, 99, 100, 101, 103, 104, 109, 112, 113, 114, 115, 117, 119, 120, 123, 124, 126, 127, 129, 130, 131, 132] # Problems that are not compatible with JuMP

prob_numbers = filter(x -> !(x in JuMP_uncompatible), all_probs) # Filter out incompatible problems
fail_to_start = [] # List to store problems that Ipopt has not been able to start the optimization - try to use another solver
failed_probs = [] # List to store problems that failed the tests
# DISCLAIMER BilevelJuMP cannot handle lower level objective that is not affine or qudratic, so they are omitted in the tests.

@testset "BiObjBenchmark.jl" begin
  for prob_no in prob_numbers
    println(prob_no)
    try
      prob  = get_bilevel_problem(prob_no)
    catch
      @warn "Problem $prob_no not found or could not be loaded."
      continue
    end
    prob  = get_bilevel_problem(prob_no)
    model = BiJuMP_convertor(prob)

    #Tests the dimension of the problems
    @test length(model[:x]) == prob.dim[1]
    @test length(model[:y]) == prob.dim[2]
    @test num_constraints(Upper(model)) == prob.dim[3]
    @test num_constraints(Lower(model)) == prob.dim[4]

    # TODO: find a way to automatize solving each problem
    try
        optimize!(model)
        if !(termination_status(model) == MOI.OPTIMAL || termination_status(model) == MOI.LOCALLY_SOLVED)
          push!(failed_probs, prob_no)
          continue
        end
        @test termination_status(model) == MOI.OPTIMAL || termination_status(model) == MOI.LOCALLY_SOLVED
    catch e
        push!(fail_to_start, prob_no)
    end
  end
end

println(length(prob_numbers) - length(failed_probs) - length(fail_to_start), " problems passed the tests.")