@testset "BilevelSolvers.jl" begin
    all_probs = collect(37:173)
    issued_probs = [36, 49, 50, 51, 138, 127, 131, 173]
    prob_numbers = filter(x -> !(x in issued_probs), all_probs)
    domain_errors = Int[]
    #prob_numbers = [1]
    for k in prob_numbers
        println("Problem $k")
        model = get_bilevel_problem(k)
        nx, ny = model.dim[1], model.dim[2]
        max_budget = 20

        D = hcat(Matrix(1.0I, nx, nx), Matrix(-1.0I, nx, nx))
        for subsolver in ["GridSearch", "NelderMead", "NOMAD"]
            x, y, Fbest, Historics = Bilevel_DS(model, subsolver, D; max_neval_upper = max_budget, max_neval_lower = 50, orthogonal = false, verbose = true)
            @test length(x) == nx
            @test length(y) == ny
            @test Historics[:Nhist][end] <= max_budget
            #x_star, y_star = model.sol[1:nx], model.sol[1+nx:nx+ny]
            #println(abs(model.F_func(x_star, y_star) - Historics[:Fhist][end]))
        end
    end
end