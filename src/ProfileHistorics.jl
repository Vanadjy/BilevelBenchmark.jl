export Profile_Historics

function Profile_Historics(prob_numbers::Vector{I}, algo_names::Vector{S}, bilevel_options::BilevelOptions; cons_handle::String = "PB", start_point::String = "y0") where {I <: Int, S <:String}
    ## Listing the problems ##
    n_probs = length(prob_numbers)

    ## Initialize solver options ##
    n_algos = length(algo_names)
    options1 = "CS"

    options2 = NOMADOptions(max_bb_eval = bilevel_options.max_neval_lower, quad_model_search = bilevel_options.search, direction_type = "ORTHO N+1 NEG", cons_handle = cons_handle, start_points = start_point) # ORTHO N+1 NEG and no search

    options3 = NOMADOptions(max_bb_eval = bilevel_options.max_neval_lower) # Default NOMAD (ORTHO 2N and quadratic search)
    All_options = Dict(algo_names .=> [options1, options2, options3])

    # Instantiate historic storages ##
    N_all_hists = [Dict(algo_names .=> [zeros(bilevel_options.max_neval_upper) for i in 1:n_algos]) for k in 1:n_probs] # Storing the historic of BB evaluations of each algo for each problem
    F_all_hists = [Dict(algo_names .=> [zeros(bilevel_options.max_neval_upper) for i in 1:n_algos]) for k in 1:n_probs] # Storing the historic of upper objective of each algo for each problem
    f_all_hists = [Dict(algo_names .=> [zeros(bilevel_options.max_neval_upper) for i in 1:n_algos]) for k in 1:n_probs] # Storing the historic of lower objective of each algo for each problem

    x_all_hists = [Dict(algo_names .=> [zeros(2, bilevel_options.max_neval_upper) for i in 1:n_algos]) for k in 1:n_probs] # Storing the historic of x (upper variables) of each algo for each problem
    y_all_hists = [Dict(algo_names .=> [zeros(2, bilevel_options.max_neval_upper) for i in 1:n_algos]) for k in 1:n_probs] # Storing the historic of y (lower variables) of each algo for each problem

    for prob_iter in eachindex(prob_numbers)
        k = prob_numbers[prob_iter]
        model = get_bilevel_problem(k)
        nx, ny = model.dim[1], model.dim[2]
        D = zeros(nx, 2*nx)

        for (algo, options) in All_options
            if options == "CS"
                x, y, Fbest, Historics = Bilevel_DS(model,
                                    "CS",
                                    D;
                                    bilevel_options = bilevel_options
                )

                N_all_hists[prob_iter][algo] = Historics[:Nhist]
                F_all_hists[prob_iter][algo] = Historics[:Fhist]
                f_all_hists[prob_iter][algo] = Historics[:fhist]
                x_all_hists[prob_iter][algo] = Historics[:xhist]
                y_all_hists[prob_iter][algo] = Historics[:yhist]
            else
                x, y, Fbest, Historics = Bilevel_DS(model,
                                                    "NOMAD",
                                                    D;
                                                    nomad_options = options,
                                                    bilevel_options = bilevel_options
                )

                N_all_hists[prob_iter][algo] = Historics[:Nhist]
                F_all_hists[prob_iter][algo] = Historics[:Fhist]
                f_all_hists[prob_iter][algo] = Historics[:fhist]
                x_all_hists[prob_iter][algo] = Historics[:xhist]
                y_all_hists[prob_iter][algo] = Historics[:yhist]
            end
        end
    end
    return N_all_hists, F_all_hists, f_all_hists, x_all_hists, y_all_hists
end