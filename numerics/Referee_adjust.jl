## REFEREE ONLY ON THE LAST ITERATE - JUST ONE ALGORITHM AS SELECTED REFEREE ##

function referee_final!(model::BilevelProblem, xHists, yHists, fHists, prob::Int, algo::String, referee_name::String, nomad_options::NOMADOptions; tol_ref::R = 1e-3, starter::String = "y0") where {R <: Float64}
    referee_flag = false
    x_star = xHists[prob][algo][:, end]
    f_star = fHists[prob][algo][end]

    # Reoptimization using NOMAD - our referee
    f = model.f_func
    x0y0 = model.xy0
    nx = model.dim[1]
    ny = model.dim[2]

    if referee_name == "NOMAD"
        function bb(y)
            fy = f(x_star, y)
            if model.dim[4] > 0
                g = model.g_func
                gy = g(x_star, y)
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

        pb.options.max_bb_eval = nomad_options.max_bb_eval
        pb.options.quad_model_search = nomad_options.quad_model_search
        pb.options.direction_type = nomad_options.direction_type
        pb.options.eval_queue_sort = nomad_options.eval_queue_sort # deactivate use of quadratic ordering
        pb.options.max_time = nomad_options.max_time # fix maximum execution time
        #pb.options.display_stats = nomad_options.display_stats # some display options
        pb.options.display_degree = nomad_options.display_degree # removing intermediate logs of NOMAD
        # Always solve the subproblem with NOMAD by starting at the same y0
        result = NOMAD.solve(pb, x0y0[nx+1:nx+ny])


        if result.bbo_best_feas !== nothing
            new_f = result.bbo_best_feas[1]
        else # The referee returned nothing : algo wins automatically
            new_f = Inf
        end
    elseif referee_name == "Ipopt"
    end
    
    if new_f < f_star - tol_ref # The referee found a strictly better solution than the algo
        @info "The referee found a better final solution than $algo on problem $prob"
        referee_flag = true
    end
    return referee_flag
end

function Referee_adjust(algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, referee_name::String; max_budget::Int = 300, cons_handle::String = "PB")
    algos_blames = Dict(algo_names .=> [Int[] for i in 1:length(algo_names)])
    for prob in eachindex(prob_numbers)
        model = get_bilevel_problem(prob_numbers[prob])
        options = NOMADOptions(max_bb_eval = max_budget, quad_model_search = true, cons_handle = cons_handle)
        for algo in algo_names
            if x_all_hists[prob][algo][end] !== x_all_hists[prob][algo][1] # If the algorithm did not moved from the starting point, ignore it
                flag = referee_final!(model, x_all_hists, y_all_hists, f_all_hists, prob, algo, referee_name, options)
                if flag #referee found a better LL solution than algo
                    push!(algos_blames[algo], prob_numbers[prob])
                end
            end
        end
    end
    return algos_blames
end

## REFEREE ONLY ON THE LAST ITERATE - ALL COMPARED ALGORITHMS SELECTED AS REFEREES ##

function intern_all_referee_final(model::BilevelProblem, xHists, yHists, fHists, prob::Int, algo::String, referees::Vector{String}, algos_nomad_options; tol_ref::R = 1e-3, ref_kknown::Vector{String} = ["Algo1", "Algo2", "Algo3"]) where {R <: Float64}
    referee_flag = false
    x_star = xHists[prob][algo][:, end]
    f_star = fHists[prob][algo][end]

    # Reoptimization using NOMAD - our referee
    new_f = Inf
    f = model.f_func
    x0y0 = model.xy0
    nx = model.dim[1]
    ny = model.dim[2]

    for ref_index in eachindex(referees)
        referee_name = referees[ref_index]
        (referee_name != algo) || continue
        if referee_name == "Algo$(ref_index)"
            nomad_options = algos_nomad_options[ref_index]
            function bb(y)
                fy = f(x_star, y)
                if model.dim[4] > 0
                    g = model.g_func
                    gy = g(x_star, y)
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

            pb.options.max_bb_eval = nomad_options.max_bb_eval
            pb.options.quad_model_search = nomad_options.quad_model_search
            pb.options.direction_type = nomad_options.direction_type
            pb.options.eval_queue_sort = nomad_options.eval_queue_sort # deactivate use of quadratic ordering
            pb.options.max_time = nomad_options.max_time # fix maximum execution time
            #pb.options.display_stats = nomad_options.display_stats # some display options
            pb.options.display_degree = nomad_options.display_degree # removing intermediate logs of NOMAD
            # Always solve the subproblem with NOMAD by starting at the same y0
            result = NOMAD.solve(pb, x0y0[nx+1:nx+ny])

            if result.bbo_best_feas !== nothing && new_f > (result.bbo_best_feas[1]) # Check if the referee returned a better feasible point
                new_f = result.bbo_best_feas[1]
            end
        else
            if referee_name ∉ ref_kknown
                error("Referee Error: Unknown referee name. Select a name among $(ref_kknown)")
            else
                error("Referee Error: Referee name given not supported in this version of the code.")
            end
        end
    end
    
    if new_f < f_star - tol_ref # The BEST referee found a strictly better solution than the algo
        @info "A referee found a better final solution than $algo on problem $prob"
        referee_flag = true
    end
    return referee_flag
end

function Referee_all_adjust(algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, referees::Vector{String}; max_budget::Int = 300, cons_handle::String = "PB")
    algos_blames = Dict(algo_names .=> [Int[] for i in 1:length(algo_names)])
    for prob in eachindex(prob_numbers)
        model = get_bilevel_problem(prob_numbers[prob])
        all_options = [NOMADOptions(max_bb_eval = max_budget, cons_handle = cons_handle),
                       NOMADOptions(max_bb_eval = max_budget, direction_type = "ORTHO N+1 NEG", cons_handle = cons_handle),
                       NOMADOptions(max_bb_eval = max_budget, quad_model_search = true, cons_handle = cons_handle)]
        for algo in algo_names
            if x_all_hists[prob][algo][end] !== x_all_hists[prob][algo][1] # If the algorithm did not moved from the starting point, ignore it
                flag = intern_all_referee_final(model, x_all_hists, y_all_hists, f_all_hists, prob, algo, referees, all_options)
                if flag #referee found a better LL solution than algo
                    push!(algos_blames[algo], prob_numbers[prob])
                end
            end
        end
    end
    return algos_blames
end

## REFEREE ON ALL THE ITERATES HISTORIC - ALL COMPARED ALGORITHMS SELECTED AS REFEREES ##

function intern_all_referee(k::Int, model::BilevelProblem, xHists, yHists, fHists, prob::Int, algo::String, referees::Vector{String}, algos_nomad_options; tol_ref::R = 1e-3, ref_kknown::Vector{String} = ["Algo1", "Algo2", "Algo3"]) where {R <: Float64}
    referee_flag = false
    x_star = xHists[prob][algo][:, k]
    f_star = fHists[prob][algo][k]

    # Reoptimization using NOMAD - our referee
    new_f = Inf
    f = model.f_func
    x0y0 = model.xy0
    nx = model.dim[1]
    ny = model.dim[2]

    for ref_index in eachindex(referees)
        referee_name = referees[ref_index]
        (referee_name != algo) || continue
        if referee_name == "Algo$(ref_index)"
            nomad_options = algos_nomad_options[ref_index]
            function bb(y)
                fy = f(x_star, y)
                if model.dim[4] > 0
                    g = model.g_func
                    gy = g(x_star, y)
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

            pb.options.max_bb_eval = nomad_options.max_bb_eval
            pb.options.quad_model_search = nomad_options.quad_model_search
            pb.options.direction_type = nomad_options.direction_type
            pb.options.eval_queue_sort = nomad_options.eval_queue_sort # deactivate use of quadratic ordering
            pb.options.max_time = nomad_options.max_time # fix maximum execution time
            #pb.options.display_stats = nomad_options.display_stats # some display options
            pb.options.display_degree = nomad_options.display_degree # removing intermediate logs of NOMAD
            # Always solve the subproblem with NOMAD by starting at the same y0
            result = NOMAD.solve(pb, x0y0[nx+1:nx+ny])

            if result.bbo_best_feas !== nothing && new_f > (result.bbo_best_feas[1]) # Check if the referee returned a better feasible point
                new_f = result.bbo_best_feas[1]
            end
        else
            if referee_name ∉ ref_kknown
                error("Referee Error: Unknown referee name. Select a name among $(ref_kknown)")
            else
                error("Referee Error: Referee name given not supported in this version of the code.")
            end
        end
    end
    
    if new_f < f_star - tol_ref # The BEST referee found a strictly better solution than the algo
        @info "A referee found a better final solution than $algo on problem $prob at iterate $k"
        referee_flag = true
    end
    return referee_flag
end

function Referee_all_historic_adjust!(F_all_hists_adjusted, algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, referees::Vector{String}; max_budget::Int = 300, cons_handle::String = "PB")
    for prob in eachindex(prob_numbers)
        model = get_bilevel_problem(prob_numbers[prob])
        all_options = [NOMADOptions(max_bb_eval = max_budget, cons_handle = cons_handle),
                       NOMADOptions(max_bb_eval = max_budget, direction_type = "ORTHO N+1 NEG", cons_handle = cons_handle),
                       NOMADOptions(max_bb_eval = max_budget, quad_model_search = true, cons_handle = cons_handle)]
        for algo in algo_names
            if x_all_hists[prob][algo][end] !== x_all_hists[prob][algo][1] # If the algorithm did not moved from the starting point, ignore it
                for k in eachindex(f_all_hists[prob][algo])
                    flag = intern_all_referee(k, model, x_all_hists, y_all_hists, f_all_hists, prob, algo, referees, all_options)
                    if flag #referee found at least once a better LL solution than algo
                        F_all_hists_adjusted[prob][algo][k] = Inf # Set at Inf the corresponding value in the upper objective historic 
                    end
                end
            end
        end
    end
    return F_all_hists_adjusted
end

## REFEREE FROM THE LAST ITERATE HISTORIC AND GOING BACKWARD ##

function Referee_backward_historic_adjust!(F_all_hists_adjusted, algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, referees::Vector{String}; max_budget::Int = 300, cons_handle::String = "PB")
    for prob in eachindex(prob_numbers)
        model = get_bilevel_problem(prob_numbers[prob])
        all_options = [NOMADOptions(max_bb_eval = max_budget, cons_handle = cons_handle),
                       NOMADOptions(max_bb_eval = max_budget, direction_type = "ORTHO N+1 NEG", cons_handle = cons_handle),
                       NOMADOptions(max_bb_eval = max_budget, quad_model_search = true, cons_handle = cons_handle)]
        for algo in algo_names
            if x_all_hists[prob][algo][end] !== x_all_hists[prob][algo][1] # If the algorithm did not moved from the starting point, ignore it
                flag = true
                k = length(f_all_hists[prob][algo])
                while flag && k >= 1 # Once we found an admissible point in the historic, we stop
                    flag = intern_all_referee(k, model, x_all_hists, y_all_hists, f_all_hists, prob, algo, referees, all_options)
                    if flag #referee found at least once a better LL solution than algo
                        F_all_hists_adjusted[prob][algo][k] = Inf # Set at Inf the corresponding value in the upper objective historic 
                    end
                    k -= 1
                end
            end
        end
    end
    return F_all_hists_adjusted
end