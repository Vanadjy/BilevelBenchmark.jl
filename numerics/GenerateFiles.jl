function generate_files!(prob_numbers::Vector{Int}, algo_names::Vector{String}; cons_handle::String = "PB", save::Bool = true, path::String = "/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/JLD2saves", path_package::String = "/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark")
    N_all_hists, F_all_hists, f_all_hists, x_all_hists, y_all_hists = Profile_Historics(prob_numbers, algo_names; cons_handle = cons_handle)

    if save
        cd(path)
        JLD2.save_object("N_all_hists-cons=$cons_handle.jld2", N_all_hists)
        JLD2.save_object("F_all_hists-cons=$cons_handle.jld2", F_all_hists)
        JLD2.save_object("f_all_hists-cons=$cons_handle.jld2", f_all_hists)
        JLD2.save_object("x_all_hists-cons=$cons_handle.jld2", x_all_hists)
        JLD2.save_object("y_all_hists-cons=$cons_handle.jld2", y_all_hists)
        cd(path_package)
    end
end