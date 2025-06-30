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