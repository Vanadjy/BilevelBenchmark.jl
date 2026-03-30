## REFEREE ON ALL THE ITERATES HISTORIC - ALL COMPARED ALGORITHMS SELECTED AS REFEREES ##

function referee_challenge(k::Int, model::BilevelProblem, xHists, yHists, fHists, prob::Int, algo::String, referees::Vector{Union{String, Int}}, algos_nomad_options; tol_ref::R = 1e-3, ref_kknown::Vector{String} = ["Algo1", "Algo2", "Algo3"]) where {R <: Float64}
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
        LL_solver_name = ""
        if referee_name == "Algo1"
            LL_solver_name = "CS"
        elseif referee_name == "Algo4"
            LL_solver_name = "COBYLA"
        elseif referee_name == "Algo5"
            LL_solver_name = "CMAES"
         elseif referee_name == "Algo2" || referee_name == "Algo3"
            LL_solver_name = "NOMAD"
         else
            error("Referee Error: Referee name given not supported in this version of the code.")
        end

        if LL_solver_name == "NOMAD"
            y_new, new_f, neval_lower = LL_subsolver(model, x_star, x0y0[nx+1:nx+ny], LL_solver_name; nomad_options = algos_nomad_options[referee_name], max_neval_lower = ny * 100)
        else
            y_new, new_f, neval_lower = LL_subsolver(model, x_star, x0y0[nx+1:nx+ny], LL_solver_name; max_neval_lower = ny * 100)
        end
        if new_f < f_star - tol_ref # The BEST referee found a strictly better solution than the algo
            @info "A referee found a better final solution than $algo on problem $prob at iterate $k"
            referee_flag = true
            break # If one referee found a better solution, we stop and consider that the algo is invalidated at this iterate
        end
    end
    return referee_flag
end

function EndPoint_Referee(algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, referees::Vector{Union{String, Int}}; max_budget::Int = 300, cons_handle::String = "PB")
    algos_blames = Dict(algo_names .=> [Int[] for i in 1:length(algo_names)])
    res_matrix = [Vector{Int64}[] for i in 1:length(prob_numbers), j in 1:length(algo_names)]
    ignored_pbs = Int64[]
    for prob in eachindex(prob_numbers)
        model = get_bilevel_problem(prob_numbers[prob])
        nomad_options = [NOMADOptions(max_bb_eval = max_budget, cons_handle = cons_handle),
                       NOMADOptions(max_bb_eval = max_budget, direction_type = "ORTHO N+1 NEG", cons_handle = cons_handle),
                       NOMADOptions(max_bb_eval = max_budget, quad_model_search = true, cons_handle = cons_handle),
                       NOMADOptions(max_bb_eval = max_budget, cons_handle = cons_handle),
                       NOMADOptions(max_bb_eval = max_budget, cons_handle = cons_handle)]
        for a in eachindex(algo_names)
            algo = algo_names[a]
            if x_all_hists[prob][algo][end] !== x_all_hists[prob][algo][1] # If the algorithm did not moved from the starting point, ignore it
                Random.seed!(seed)
                k = length(f_all_hists[prob][algo])
                flag = referee_challenge(k, model, x_all_hists, y_all_hists, f_all_hists, prob, algo, referees, nomad_options)
                if flag #referee found a better LL solution than algo
                    push!(algos_blames[algo], prob_numbers[prob])
                    #push!(res_matrix[prob, a], successful_refs)
                end
            else
                push!(ignored_pbs, prob)
            end
        end
    end
    #kept_pbs = setdiff(collect(1:length(prob_numbers)), ignored_pbs)
    #JLD2.save_object(joinpath("/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark/numerics/logs", "res_matrix-FinalRef.jld2"), res_matrix[kept_pbs, :])
    #JLD2.save_object(joinpath("/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark/numerics/logs", "ignored_pbs.jld2"), ignored_pbs)
    return algos_blames
end

function Complete_Referee!(F_all_hists_adjusted, N_all_hists_adjusted, algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, referees::Vector{Union{String, Int}}; max_budget::Int = 300, cons_handle::String = "PB")
    for prob in eachindex(prob_numbers)
        model = get_bilevel_problem(prob_numbers[prob])
        nomad_options = Dict("Algo2" => NOMADOptions(max_bb_eval = max_budget, direction_type = "ORTHO N+1 NEG", cons_handle = cons_handle),
                              "Algo3" => NOMADOptions(max_bb_eval = max_budget, quad_model_search = true, cons_handle = cons_handle))
        for algo in algo_names
            if x_all_hists[prob][algo][end] !== x_all_hists[prob][algo][1] # If the algorithm did not moved from the starting point, ignore it
                k = length(f_all_hists[prob][algo])
                while k >= 1
                    Random.seed!(seed)
                    flag = referee_challenge(k, model, x_all_hists, y_all_hists, f_all_hists, prob, algo, referees, nomad_options)
                    if flag #referee found at least once a better LL solution than algo
                        F_all_hists_adjusted[prob][algo][k] = NaN # Set at Inf the corresponding value in the upper objective historic
                        N_all_hists_adjusted[prob][algo][k] = NaN # Set at Inf the corresponding value in the lower objective historic
                    end
                    k -= 1
                end
            end
        end
    end
    orphans = Dict{String, Vector{Int}}(key => Int[] for key in algo_names)
    # After the process, remove all Inf entries from the historics
    for prob in eachindex(prob_numbers)
        for algo in algo_names
            filter!(!isnan, F_all_hists_adjusted[prob][algo])
            filter!(!isnan, N_all_hists_adjusted[prob][algo])
            if length(F_all_hists_adjusted[prob][algo]) == 0 # If the referee invalidated all the historic, count it as an orphaned run
                push!(orphans[algo], prob_numbers[prob])
            end
        end
    end
    #JLD2.save_object(joinpath("/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark/numerics/logs", "orphans-AllRef.jld2"), orphans)
    return F_all_hists_adjusted, N_all_hists_adjusted
end

## REFEREE FROM THE LAST ITERATE HISTORIC AND GOING BACKWARD ##

function Reverse_Referee!(F_all_hists_adjusted, N_all_hists_adjusted, algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, referees::Vector{Union{String, Int}}; max_budget::Int = 300, cons_handle::String = "PB")
    for prob in eachindex(prob_numbers)
        model = get_bilevel_problem(prob_numbers[prob])
        nomad_options = [NOMADOptions(max_bb_eval = max_budget, cons_handle = cons_handle),
                       NOMADOptions(max_bb_eval = max_budget, direction_type = "ORTHO N+1 NEG", cons_handle = cons_handle),
                       NOMADOptions(max_bb_eval = max_budget, quad_model_search = true, cons_handle = cons_handle),
                       NOMADOptions(max_bb_eval = max_budget, cons_handle = cons_handle),
                       NOMADOptions(max_bb_eval = max_budget, cons_handle = cons_handle)]
        for algo in algo_names
            if x_all_hists[prob][algo][end] !== x_all_hists[prob][algo][1] # If the algorithm did not moved from the starting point, ignore it
                flag = true
                k = length(f_all_hists[prob][algo])
                while flag && k >= 1 # Once we found an admissible point in the historic, we stop
                    Random.seed!(seed)
                    flag = referee_challenge(k, model, x_all_hists, y_all_hists, f_all_hists, prob, algo, referees, nomad_options)
                    if flag #referee found at least once a better LL solution than algo
                        F_all_hists_adjusted[prob][algo][k] = NaN # Set at Inf the corresponding value in the upper objective historic 
                        N_all_hists_adjusted[prob][algo][k] = NaN # Set at Inf the corresponding value in the lower objective historic
                    end
                    k -= 1
                end
            end
        end
    end
    orphans = Dict{String, Vector{Int}}(key => Int[] for key in algo_names)
    # After the process, remove all Inf entries from the historics
    for prob in eachindex(prob_numbers)
        for algo in algo_names
            filter!(!isnan, F_all_hists_adjusted[prob][algo])
            filter!(!isnan, N_all_hists_adjusted[prob][algo])
            if length(F_all_hists_adjusted[prob][algo]) == 0 # If the referee invalidated all the historic, count it as an orphaned run
                push!(orphans[algo], prob_numbers[prob])
            end
        end
    end
    #JLD2.save_object(joinpath("/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark/numerics/logs", "orphans-BackwardRef.jld2"), orphans)
    return F_all_hists_adjusted, N_all_hists_adjusted
end