@testset "BilevelSolvers.jl" begin
    for k in 1:162
        model = get_bilevel_problem(k)
        nx, ny = model.dim[1], model.dim[2]
        max_budget = 200

        D = hcat(Matrix(1.0I, nx, nx), Matrix(-1.0I, nx, nx))
        for subsolver in ["GridSearch", "NelderMead"]
            x, y, Fbest, Historics = Bilevel_DS(model, subsolver, D; max_neval_upper = max_budget, max_neval_lower = 50, verbose = false)
            @test length(x) == nx
            @test length(y) == ny
            @test Historics[:Nhist][end] == max_budget
            x_star, y_star = model.sol[1:nx], model.sol[1+nx:nx+ny]
            println(abs(model.F_func(x_star, y_star) - Historics[:Fhist][end]))
        end
    end
end