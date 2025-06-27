@testset "BiObjBenchmark.jl" begin
    model = get_bilevel_problem(1)
    nx = model.nx
    max_budget = 200

    D = hcat(Matrix(1.0I, nx, nx), Matrix(-1.0I, nx, nx))
    for subsolver in ["GridSearch", "NOMAD"]
        x, y, Fbest, Historics = Bilevel_DS(model, subsolver, D; max_neval_upper = max_budget, max_neval_lower = Int(max_budget/10))
        @test length(x) == model.nx
        @test length(y) == model.ny
        @test Historics[:Nhist][end] == max_budget
    end
end