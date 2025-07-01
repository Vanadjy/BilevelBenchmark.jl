export GridSearch_subsolver, Bilevel_DS

using NOMAD, Printf

#using Optimization, OptimizationNOMAD, Optim, OptimizationOptimJL

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


function Bilevel_DS(model::BilevelProblem,
                    subsolver::String,
                    D::Matrix{Float64};
                    Δ0::Float64 = 1.0,
                    γ::Float64 = 1/2,
                    oppportunistic::Bool = true,
                    ordered::Bool = false,
                    search::Bool = false,
                    orthogonal::Bool = true,
                    max_neval_upper::Int = 1000,
                    max_neval_upper_cons::Int = 1000,
                    max_neval_lower::Int = 100,
                    tol_upper::Float64 = 1e-6,
                    tol_lower::Float64 = 1e-6,
                    max_time::Float64 = 3600.0,
                    verbose::Bool = true
    )
    subsolver_avail = ["GridSearch", "NOMAD", "Ipopt", "NelderMead"]
    start = time()
    elapsed_time = time() - start

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
                D[:, j + nx] .= (-1.0) * round.((Δk/δk*norm(H[:, j], Inf)) * H[:, j])
            end
            display(D) #TODO: check back OrthoMADS: directions just become bigger and bigger
        else
            #@warn "No other way to build dense directions has been implemented yet."
            for j = 1:nx
                D[j,j] = 1.0
                D[j,j+nx] = -1.0
            end
        end
        while (i < size(D, 2)) && !(stop_poll) && (neval_upper < max_neval_upper) && (neval_upper_cons < max_neval_upper_cons)
            i += 1

            # Compute optimal answer of lower level on each poll point
            if subsolver == "GridSearch"
                t = xk + δk * D[:, i]
                # Apply Grid Search on the local frame of size Δk
                yk_new, fk_new = GridSearch_subsolver(f, g, t, yk .- Δk*ones(eltype(xk), ny), yk .+ Δk*ones(eltype(xk), ny); num_points = max_neval_lower)
            elseif subsolver == "NelderMead"
                y0 = model.xy0[nx+1:nx+ny]
                t = xk + δk * D[:, i]
                f_bb(y, t) = (any(>(0), g(t, y))) ? 1e16 : f(t, y)
                bb_func = OptimizationFunction(f_bb)
                prob = OptimizationProblem(bb_func, y0, t)
                sol = Optimization.solve(prob, Optim.NelderMead())
                yk_new .= sol.u
                fk_new = sol.objective
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
                    pb = NomadProblem(ny, 2, ["OBJ", "PB"], bb)
                else
                    pb = NomadProblem(ny, 1, ["OBJ"], bb)
                end


                # Always solve the subproblem with NOMAD by starting at the same y0
                pb.options.max_bb_eval = max_neval_lower
                pb.options.quad_model_search = false # deactivate quadratic model subproblem resolution
                pb.options.direction_type = "ORTHO N+1 NEG"
                pb.options.eval_queue_sort = "DIR_LAST_SUCCESS" # deactivate use of quadratic ordering
                pb.options.max_time = max_time # fix maximum execution time
                #pb.options.display_stats = ["EVAL", "SOL", "OBJ"] # some display options
                pb.options.display_degree = 0  # removing intermediate logs of NOMAD

                result = NOMAD.solve(pb, x0y0[nx+1:nx+ny])
                if result.x_best_feas !== nothing
                    yk_new .= result.x_best_feas
                    fk_new = result.bbo_best_feas[1]
                end
            #elseif solver == "Ipopt"
            else
                @error "Subsolver $subsolver is not known or implemented. Rather try one of the subsolvers among $subsolver_avail"
            end

            neval_upper += 1
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
                     :xhist => x_hist[1:k],
                     :yhist => x_hist[1:k]
        )
    return xk, yk, Fbest, Historics
end