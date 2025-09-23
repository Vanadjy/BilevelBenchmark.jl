export GridSearch_subsolver, LL_subsolver, Bilevel_DS_WithoutCoupling, Biphase_feas, Bilevel_DS, CS_subsolver

using NOMAD, Printf

using Optimization, OptimizationNOMAD, Optim, OptimizationOptimJL

function GridSearch_subsolver(f,
                              g,
                              x::Vector{Float64},
                              l_bounds::Vector{Float64}, 
                              u_bounds::Vector{Float64}; 
                              num_points::Int = 500)
    
    n = length(l_bounds)
    δs = (u_bounds .- l_bounds) ./ (num_points - 1)
    Grid = zeros(n, num_points)
    fbest = Inf

    @assert length(u_bounds) == n "Lower and upper bounds must have the same dimension."

    for z = 1:num_points
        Grid[:, z] .= l_bounds .+ (z .* δs)
    end

    ybest = similar(l_bounds)
    
    for y in eachcol(Grid)
        if (any(>(0), g(x, y))) # Infeasible point
            continue
        else
            fval = f(x, y)
            if fval < fbest
                fbest = fval
                ybest .= y
            end
        end
    end

    return ybest, fbest
end

function CS_subsolver(f, g, y0, t; neval_tot = 100, δ0 = 1.0, ϵ = 1e-6, τ = 1/2, poll = "")
    xk = y0
    fbest = f(t, y0)
    neval_f = 1
    δk = δ0

    Id = Matrix(I, length(y0), length(y0))
    D = [Id -Id]
    P = [D[:,i] for i in axes(D, 2)]

    while (neval_f ≤ neval_tot) #&& (δk > ϵ)
        found_upgrade = false
        i = 1
        while !found_upgrade && (neval_f ≤ neval_tot) && (i ≤ length(P))
            Pi = xk + (δk * P[i])
            f_Pi = f(t, Pi)
            neval_f += 1
            if (f_Pi < fbest) && (all(<=(0), g(t, Pi))) # feasible upgrade of f found in the polling : success (EB here)
                fbest = f_Pi
                xk = Pi

                #reordering the poll so that the successful direction is in first
                if i > 1
                    P = vcat(circshift!(P[1:i], 1), P[i+1:end])
                end

                found_upgrade = true
            end
            i += 1
        end

        # polling step have failed : xk is a local min for the mesh
        if !found_upgrade
            δk *= τ
        end
    end
    xk, fbest
end

function LL_subsolver(model::BilevelProblem,
                      i::Int,
                      xk,
                      yk,
                      subsolver::String,
                      D::Matrix{Float64};
                      Δk::Float64 = 1.0,
                      δk::Float64 = 1/2,
                      nomad_options::NOMADOptions = NOMADOptions(),
                      max_neval_lower::Int = 100
    )
    subsolver_avail = ["GridSearch", "NOMAD", "Ipopt", "CS"]
    nx = model.dim[1]
    ny = model.dim[2]
    yk_new = similar(yk)
    fk_new = zero(eltype(xk))
    # Compute optimal answer of lower level on each poll point
    if subsolver == "GridSearch"
        t = xk + δk * D[:, i]
        # Apply Grid Search on the local frame of size Δk
        yk_new, fk_new = GridSearch_subsolver(model.f_func, model.g_func, t, yk .- Δk*ones(eltype(xk), ny), yk .+ Δk*ones(eltype(xk), ny); num_points = max_neval_lower)
    elseif subsolver == "NelderMead"
        y0 = model.xy0[nx+1:nx+ny]
        t = xk + δk * D[:, i]
        f_bb(y, t) = (any(>(0), g(t, y))) ? 1e16 : f(t, y)
        bb_func = OptimizationFunction(f_bb)
        prob = OptimizationProblem(bb_func, y0, t)
        sol = Optimization.solve(prob, Optim.NelderMead())
        yk_new .= sol.u
        fk_new = sol.objective
    elseif subsolver == "CS"
        t = xk + δk * D[:, i]
        yk_new, fk_new = CS_subsolver(model.f_func, model.g_func, yk, t; neval_tot = max_neval_lower, δ0 = Δk)
    elseif subsolver == "NOMAD"
        # Apply NOMAD solver
        t = xk + δk * D[:, i]
        f = model.f_func
        function bb(y)
            fy = f(t, y)
            if model.dim[4] > 0
                g = model.g_func
                gy = g(t, y)
                bb_outputs = [fy; gy]
            else
                bb_outputs = [fy]
            end
            success = true
            count_eval = true
            return (success, count_eval, bb_outputs)
        end
        if model.dim[4] > 0
            pb = NomadProblem(ny, 2, ["OBJ", nomad_options.cons_handle], bb)
        else
            pb = NomadProblem(ny, 1, ["OBJ"], bb)
        end

        pb.options.max_bb_eval = max_neval_lower
        pb.options.quad_model_search = nomad_options.quad_model_search
        pb.options.sgtelib_model_search = false
        pb.options.speculative_search = false
        pb.options.nm_search = false
        pb.options.vns_mads_search = false
        pb.options.direction_type = nomad_options.direction_type
        pb.options.eval_queue_sort = nomad_options.eval_queue_sort # deactivate use of quadratic ordering
        pb.options.max_time = nomad_options.max_time # fix maximum execution time
        #pb.options.display_stats = nomad_options.display_stats # some display options
        pb.options.display_degree = nomad_options.display_degree # removing intermediate logs of NOMAD

        # Always solve the subproblem with NOMAD by starting at the same y0
        if nomad_options.start_points == "y0"
            result = NOMAD.solve(pb, model.xy0[nx+1:nx+ny])
            if result.x_best_feas !== nothing
                yk_new .= result.x_best_feas
                fk_new = result.bbo_best_feas[1]
            else
                yk_new .= yk
                fk_new = model.f_func(xk, yk)
            end
        elseif nomad_options.start_points == "yk-1"
            result = NOMAD.solve(pb, yk)
            if result.x_best_feas !== nothing
                yk_new .= result.x_best_feas
                fk_new = result.bbo_best_feas[1]
            else
                yk_new .= yk
                fk_new = model.f_func(xk, yk)
            end
        else
            @error "Start points must be either 'y0' or 'yk-1'. Other start points are not supported yet."
        end
    else
        @error "Subsolver $subsolver is not known or implemented. Try one of the subsolvers among $subsolver_avail"
    end
    return yk_new, fk_new
end

function Bilevel_DS_WithoutCoupling(model::BilevelProblem,
                                    F,
                                    f,
                                    g,
                                    subsolver::String,
                                    D::Matrix{Float64};
                                    bilevel_options::BilevelOptions = BilevelOptions(),
                                    nomad_options::NOMADOptions = NOMADOptions()
    )
    start = time()
    elapsed_time = time() - start

    #Initialization from the bilevel options
    max_neval_upper = bilevel_options.max_neval_upper
    max_neval_upper_cons = bilevel_options.max_neval_upper_cons
    max_neval_lower = bilevel_options.max_neval_lower
    Δ0 = bilevel_options.Δ0
    γ = bilevel_options.γ
    oppportunistic = bilevel_options.oppportunistic
    ordered = bilevel_options.ordered
    search = bilevel_options.search
    orthogonal = bilevel_options.orthogonal
    max_time = bilevel_options.max_time
    verbose = bilevel_options.verbose

    @assert γ > 0 "The parameter γ must be positive."
    @assert γ < 1 "The parameter γ must be less than 1."
    @assert Δ0 > 0 "The initial step size Δ0 must be positive."
    @assert ordered ≤ oppportunistic "The poll cannot be ordered if we don't apply oppportunistic scheme."
    @assert size(D, 2) == 2*model.dim[1] "Set of poll directions need to be a maximal positive basis"

    #Initialization from the model
    nx = model.dim[1]
    ny = model.dim[2]
    x0y0 = model.xy0
    xk = x0y0[1:nx]
    yk = x0y0[nx+1:nx+ny]

    Fk = F(xk, yk)
    fk = f(xk, yk)

    #Initialize optimization parameters
    Δk = Δ0
    δk = min(Δk, Δk^2)
    Fbest = Inf
    poll_improvement = false
    neval_upper = 1
    neval_upper_cons = 1

    if verbose > 0
        #! format: off
        @info @sprintf "%6s %8s %8s %7s %1s %7s" "N" "F(x,y)" "f(x,y)" "Δk" "poll" "time [s]"
        #! format: on
    end

    k = 1

    while !(neval_upper ≥ max_neval_upper || neval_upper_cons ≥ max_neval_upper_cons || elapsed_time ≥ max_time)
        ## ------------------ Updating Mesh parameter ------------------ ##
        δk = min(Δk, Δk^2)

        ## ------------------ Poll step ------------------ ##
        i = 0
        poll_improvement = false
        stop_poll = false

        ## ------------------ Generate Poll directions ------------------ ##

        H = zeros(eltype(xk), nx, nx)
        yk_new = similar(yk)
        fk_new = zero(eltype(xk))

        if orthogonal
            # Generate a random vector
            v = rand(nx)
            v /= norm(v)

            # Generate Househodler matrix
            Householder!(H, v)
            for j = 1:nx
                D[:, j] .= round.((Δk/δk*norm(H[:, j], Inf)) * H[:, j])
                D[:, j + nx] .= (-1.0) * round.((Δk/(δk*norm(H[:, j], Inf))) * H[:, j])
            end
        else
            #@warn "No other way to build dense directions has been implemented yet."
            for j = 1:nx
                D[j,j] = 1.0
                D[j,j+nx] = -1.0
            end
        end
        while (i < size(D, 2)) && !(stop_poll) && (neval_upper < max_neval_upper) && (neval_upper_cons < max_neval_upper_cons)
            i += 1

            t = xk + δk * D[:, i]
            yk_new, fk_new = LL_subsolver(model, i, xk, yk, subsolver, D; Δk = Δk, δk = δk, nomad_options = nomad_options, max_neval_lower = ny * max_neval_lower)

            neval_upper += 1
            Fk_new = F(t, yk_new)
            if (Fk_new < Fk) # Successful iteration
                poll_improvement = true
                xk .= t
                yk .= yk_new
                Fk = Fk_new
                fk = fk_new
                Δk /= γ #Increase Mesh size parameter

                # Apply poll strategies for successful iterations
                if oppportunistic
                    if ordered
                        d_temp = similar(xk)
                        d_temp .= D[:, i]
                        for j in 2:i
                            D[:, j] .= D[:, j-1]
                        end
                        D[:, 1] .= d_temp
                    end
                    stop_poll = true
                    continue
                end
            end
        end # end of Poll
        
        if !poll_improvement # Unsuccessful iteration
            Δk *= γ
        end

        poll_status = poll_improvement ? "succ" : "fail"
        elapsed_time = time() - start

        if verbose > 0
            #! format: off
            @info @sprintf "%6d %8.2e %8.2e %7.1e %1s %7.1e" neval_upper Fk fk Δk poll_status elapsed_time
            #! format: on
        end

        # Update iterate
        k += 1
    end
    return xk, yk, Fbest
end

function Biphase_LL(model, nomad_options, tol)
    g = model.g_func
    t = model.xy0[1:model.dim[1]]
    function bb(y)
        hy = sum(max.(g(t, y), 0).^2)
        bb_outputs = [hy]
        success = true
        count_eval = true
        return (success, count_eval, bb_outputs)
    end
    pb = NomadProblem(model.dim[2], 1, ["OBJ"], bb)

    pb.options.max_bb_eval = nomad_options.max_bb_eval
    pb.options.quad_model_search = nomad_options.quad_model_search
    pb.options.sgtelib_model_search = false
    pb.options.speculative_search = false
    pb.options.nm_search = false
    pb.options.vns_mads_search = false
    pb.options.direction_type = nomad_options.direction_type
    pb.options.eval_queue_sort = nomad_options.eval_queue_sort # deactivate use of quadratic ordering
    pb.options.max_time = nomad_options.max_time # fix maximum execution time
    #pb.options.display_stats = nomad_options.display_stats # some display options
    pb.options.display_degree = nomad_options.display_degree # removing intermediate logs of NOMAD

    result = NOMAD.solve(pb, model.xy0[model.dim[1]+1:model.dim[1]+model.dim[2]])
    if result.bbo_best_feas[1] > tol
        @warn "No feasible point found for the first phase for the lower level biphase."
    end
    return result.x_best_feas, result.bbo_best_feas[1]
end

function Biphase_feas(model, bilevel_options, nomad_options, D, tol)
    @info "Starting bilevel first phase to find a feasible point."
    G = model.G_func
    g = model.g_func
    h(x,y) = sum(max.(G(x,y), 0).^2)
    xk, yk, Fbest = Bilevel_DS_WithoutCoupling(model, h, model.f_func, g, "NOMAD", D;
                                                bilevel_options = bilevel_options,
                                                nomad_options = nomad_options
    )
    if Fbest > tol
        @warn "No feasible point found for the first phase."
    end
    return xk, yk
end


function Bilevel_DS(model::BilevelProblem,
                    subsolver::String,
                    D::Matrix{Float64};
                    bilevel_options::BilevelOptions = BilevelOptions(),
                    nomad_options::NOMADOptions = NOMADOptions()
    )
    start = time()
    elapsed_time = time() - start

    #Initialization from the bilevel options
    max_neval_upper = bilevel_options.max_neval_upper
    max_neval_upper_cons = bilevel_options.max_neval_upper_cons
    max_neval_lower = bilevel_options.max_neval_lower
    Δ0 = bilevel_options.Δ0
    γ = bilevel_options.γ
    oppportunistic = bilevel_options.oppportunistic
    ordered = bilevel_options.ordered
    search = bilevel_options.search
    orthogonal = bilevel_options.orthogonal
    max_time = bilevel_options.max_time
    biphase = bilevel_options.biphase
    verbose = bilevel_options.verbose

    @assert γ > 0 "The parameter γ must be positive."
    @assert γ < 1 "The parameter γ must be less than 1."
    @assert Δ0 > 0 "The initial step size Δ0 must be positive."
    @assert ordered ≤ oppportunistic "The poll cannot be ordered if we don't apply oppportunistic scheme."
    @assert size(D, 2) == 2*model.dim[1] "Set of poll directions need to be a maximal positive basis"

    #Initialization
    nx = model.dim[1]
    ny = model.dim[2]
    x0y0 = model.xy0
    xk = x0y0[1:nx]
    yk = x0y0[nx+1:nx+ny]

    #Define functions relative to the model
    F = model.F_func
    f = model.f_func
    G = model.G_func
    g = model.g_func

    Fk = F(xk, yk)
    fk = f(xk, yk)

    #Initialize optimization parameters
    Δk = Δ0
    δk = min(Δk, Δk^2)
    Fbest = Inf
    poll_improvement = false
    neval_upper = 1
    neval_upper_cons = 1

    #Declare historics
    Neval_upper_hist = zeros(max(max_neval_upper, max_neval_upper_cons))
    F_hist = zeros(max(max_neval_upper, max_neval_upper_cons))
    f_hist = zeros(max(max_neval_upper, max_neval_upper_cons))
    x_hist = zeros(nx, max(max_neval_upper, max_neval_upper_cons))
    y_hist = zeros(ny, max(max_neval_upper, max_neval_upper_cons))

    #Initialize historics
    Neval_upper_hist[1] = neval_upper
    F_hist[1] = Fk
    f_hist[1] = fk
    x_hist[:, 1] .= xk
    y_hist[:, 1] .= yk


    if (model.dim[4] > 0) && any(>(0), g(xk, yk)) # Need a single level first phase to find a feasible LL point
        if biphase
            @info "Starting first phase to find a feasible LL point."
            yk, fk = Biphase_LL(model, nomad_options, 1e-6)
        end
    end

    if (model.dim[3] > 0) && (any(>(0), G(xk, yk))) # If the model has upper level couppling constraints and are violated
        if biphase
            xk, yk = Biphase_feas(model, bilevel_options, nomad_options, D, 1e-6)
        end
    end

    if verbose > 0
        #! format: off
        @info @sprintf "%6s %8s %8s %7s %1s %7s" "N" "F(x,y)" "f(x,y)" "Δk" "poll" "time [s]"
        #! format: on
    end

    k = 1

    while !(neval_upper ≥ max_neval_upper || neval_upper_cons ≥ max_neval_upper_cons || elapsed_time ≥ max_time)
        ## ------------------ Updating Mesh parameter ------------------ ##
        δk = min(Δk, Δk^2)

        ## ------------------ Poll step ------------------ ##
        i = 0
        poll_improvement = false
        stop_poll = false

        ## ------------------ Generate Poll directions ------------------ ##

        H = zeros(eltype(xk), nx, nx)
        yk_new = similar(yk)
        fk_new = zero(eltype(xk))

        if orthogonal
            # Generate a random vector
            v = rand(nx)
            v /= norm(v)

            # Generate Househodler matrix
            Householder!(H, v)
            for j = 1:nx
                D[:, j] .= round.((Δk/δk*norm(H[:, j], Inf)) * H[:, j])
                D[:, j + nx] .= (-1.0) * round.((Δk/(δk*norm(H[:, j], Inf))) * H[:, j])
            end
        else
            #@warn "No other way to build dense directions has been implemented yet."
            for j = 1:nx
                D[j,j] = 1.0
                D[j,j+nx] = -1.0
            end
        end
        while (i < size(D, 2)) && !(stop_poll) && (neval_upper < max_neval_upper) && (neval_upper_cons < max_neval_upper_cons)
            i += 1

            yk_new, fk_new = LL_subsolver(model, i, xk, yk, subsolver, D; Δk = Δk, δk = δk, nomad_options = nomad_options, max_neval_lower = ny * max_neval_lower)

            neval_upper += 1
            t = xk + δk * D[:, i]
            Fk_new = F(t, yk_new)

            # Check feasibility
            Gk = G(t, yk_new)
            neval_upper_cons += 1
            if (any(>(0), Gk))
                #@info "Infeasible point : don't call the blackbox"
                continue # Infeasible point
            end
            if (Fk_new < Fk) # Successful iteration
                poll_improvement = true
                xk .= t
                yk .= yk_new
                Fk = Fk_new
                fk = fk_new
                Δk /= γ #Increase Mesh size parameter

                # Apply poll strategies for successful iterations
                if oppportunistic
                    if ordered
                        d_temp = similar(xk)
                        d_temp .= D[:, i]
                        for j in 2:i
                            D[:, j] .= D[:, j-1]
                        end
                        D[:, 1] .= d_temp
                    end
                    stop_poll = true
                    continue
                end
            end
        end # end of Poll
        
        if !poll_improvement # Unsuccessful iteration
            Δk *= γ
        end

        poll_status = poll_improvement ? "succ" : "fail"
        elapsed_time = time() - start

        if verbose > 0
            #! format: off
            @info @sprintf "%6d %8.2e %8.2e %7.1e %1s %7.1e" neval_upper Fk fk Δk poll_status elapsed_time
            #! format: on
        end

        # Update historics
        k += 1
        Neval_upper_hist[k] = neval_upper
        F_hist[k] = Fk
        f_hist[k] = fk
        x_hist[:, k] .= xk
        y_hist[:, k] .= yk

    end
    Historics = Dict(:Nhist => Neval_upper_hist[1:k],
                     :Fhist => F_hist[1:k],
                     :fhist => f_hist[1:k],
                     :xhist => x_hist[:, 1:k],
                     :yhist => x_hist[:, 1:k]
        )
    return xk, yk, Fbest, Historics
end