using BilevelBenchmark

using Test, LinearAlgebra # Basic Julia packages
using JuMP, BilevelJuMP # JuMP API
using Ipopt#, HiGHS, MibS_jll # Solvers

include("BOLIB_test.jl") # Include the test file