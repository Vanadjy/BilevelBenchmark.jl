export Profile_Historics

using JLD2

function Profile_Historics(prob_numbers::Vector{I}, algo_names::Vector{Union{S, I}}, bilevel_options::BilevelOptions, subsolver_options; cons_handle::String = "PB", start_point::String = "y0") where {I <: Int, S <:String}
    ## Listing the problems ##
    n_probs = length(prob_numbers)
    n_algos = length(algo_names)

    # Instantiate historic storages ##
    N_UL_all_hists = [Dict(algo_names .=> [zeros(bilevel_options.max_neval_upper) for i in 1:n_algos]) for k in 1:n_probs] # Storing the historic of BB evaluations of each algo for each problem
    N_LL_all_hists = [Dict(algo_names .=> [zeros(bilevel_options.max_neval_upper) for i in 1:n_algos]) for k in 1:n_probs] # Storing the historic of BB evaluations of each algo for each problem
    F_all_hists = [Dict(algo_names .=> [zeros(bilevel_options.max_neval_upper) for i in 1:n_algos]) for k in 1:n_probs] # Storing the historic of upper objective of each algo for each problem
    f_all_hists = [Dict(algo_names .=> [zeros(bilevel_options.max_neval_upper) for i in 1:n_algos]) for k in 1:n_probs] # Storing the historic of lower objective of each algo for each problem

    x_all_hists = [Dict(algo_names .=> [zeros(2, bilevel_options.max_neval_upper) for i in 1:n_algos]) for k in 1:n_probs] # Storing the historic of x (upper variables) of each algo for each problem
    y_all_hists = [Dict(algo_names .=> [zeros(2, bilevel_options.max_neval_upper) for i in 1:n_algos]) for k in 1:n_probs] # Storing the historic of y (lower variables) of each algo for each problem
    t_all_hists = [Dict(algo_names .=> [zeros(bilevel_options.max_neval_upper) for i in 1:n_algos]) for k in 1:n_probs] # Storing the historic of elapsed time of each algo for each problem

    lambda_list = load_object("numerics/logs/lambda_list.jld2")
    for prob_iter in eachindex(prob_numbers)
        k = prob_numbers[prob_iter]
        model = get_bilevel_problem(k)
        ωk = lambda_list[prob_iter]
        nx, ny = model.dim[1], model.dim[2]
        D = zeros(nx, 2*nx)

        for algo in algo_names
            @info "Running $algo on problem $k"
            options = subsolver_options[algo]
            if options == "CS"
                if algo == "Algo11"
                    bilevel_options.max_neval_lower = 1000*100
                end
                x, y, Fbest, Historics = Bilevel_DS(model,
                                    options,
                                    ωk,
                                    D;
                                    bilevel_options = bilevel_options
                )

                N_UL_all_hists[prob_iter][algo] = Historics[:N_UL_hist]
                N_LL_all_hists[prob_iter][algo] = Historics[:N_LL_hist]
                F_all_hists[prob_iter][algo] = Historics[:Fhist]
                f_all_hists[prob_iter][algo] = Historics[:fhist]
                x_all_hists[prob_iter][algo] = Historics[:xhist]
                y_all_hists[prob_iter][algo] = Historics[:yhist]
                t_all_hists[prob_iter][algo] = Historics[:thist]
            elseif options == "COBYLA"
                if algo == "Algo10"
                    bilevel_options.max_neval_lower = 100*lower_budget
                end
                x, y, Fbest, Historics = Bilevel_DS(model,
                                    options,
                                    ωk,
                                    D;
                                    bilevel_options = bilevel_options
                )

                N_UL_all_hists[prob_iter][algo] = Historics[:N_UL_hist]
                N_LL_all_hists[prob_iter][algo] = Historics[:N_LL_hist]
                F_all_hists[prob_iter][algo] = Historics[:Fhist]
                f_all_hists[prob_iter][algo] = Historics[:fhist]
                x_all_hists[prob_iter][algo] = Historics[:xhist]
                y_all_hists[prob_iter][algo] = Historics[:yhist]
                t_all_hists[prob_iter][algo] = Historics[:thist]
            elseif options == "CMAES"
                if algo == "Algo9"
                    bilevel_options.max_neval_lower = 100*lower_budget
                end
                x, y, Fbest, Historics = Bilevel_DS(model,
                                    options,
                                    ωk,
                                    D;
                                    bilevel_options = bilevel_options
                )

                N_UL_all_hists[prob_iter][algo] = Historics[:N_UL_hist]
                N_LL_all_hists[prob_iter][algo] = Historics[:N_LL_hist]
                F_all_hists[prob_iter][algo] = Historics[:Fhist]
                f_all_hists[prob_iter][algo] = Historics[:fhist]
                x_all_hists[prob_iter][algo] = Historics[:xhist]
                y_all_hists[prob_iter][algo] = Historics[:yhist]
                t_all_hists[prob_iter][algo] = Historics[:thist]
            elseif options == "LHS"
                x, y, Fbest, Historics = Bilevel_DS(model,
                                    options,
                                    ωk,
                                    D;
                                    bilevel_options = bilevel_options
                )

                N_UL_all_hists[prob_iter][algo] = Historics[:N_UL_hist]
                N_LL_all_hists[prob_iter][algo] = Historics[:N_LL_hist]
                F_all_hists[prob_iter][algo] = Historics[:Fhist]
                f_all_hists[prob_iter][algo] = Historics[:fhist]
                x_all_hists[prob_iter][algo] = Historics[:xhist]
                y_all_hists[prob_iter][algo] = Historics[:yhist]
                t_all_hists[prob_iter][algo] = Historics[:thist]
            else
                x, y, Fbest, Historics = Bilevel_DS(model,
                                                    "NOMAD",
                                                    ωk,
                                                    D;
                                                    nomad_options = options,
                                                    bilevel_options = bilevel_options
                )

                N_UL_all_hists[prob_iter][algo] = Historics[:N_UL_hist]
                N_LL_all_hists[prob_iter][algo] = Historics[:N_LL_hist]                
                F_all_hists[prob_iter][algo] = Historics[:Fhist]
                f_all_hists[prob_iter][algo] = Historics[:fhist]
                x_all_hists[prob_iter][algo] = Historics[:xhist]
                y_all_hists[prob_iter][algo] = Historics[:yhist]
                t_all_hists[prob_iter][algo] = Historics[:thist]
            end
        end
    end
    return N_UL_all_hists, N_LL_all_hists, F_all_hists, f_all_hists, x_all_hists, y_all_hists, t_all_hists
end