using BilevelBenchmark

using Test, LinearAlgebra # Basic Julia packages
using JuMP, BilevelJuMP # JuMP API
using Ipopt#, HiGHS, MibS_jll # Solvers

include("test_utils.jl")
include("BOLIB_test.jl")
include("test_solvers.jl")
include("test_struct.jl")