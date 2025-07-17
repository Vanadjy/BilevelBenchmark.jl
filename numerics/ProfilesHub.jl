using BilevelBenchmark
using PGFPlots, JLD2, NOMAD

# To check the profiles
using BenchmarkProfiles
using Plots

include("plot_settings.jl")
include("plot-utils.jl")
include("GenerateFiles.jl")
include("DrawProfiles.jl")
include("Referee_adjust.jl")

all_probs = collect(1:173)
issued_probs = [36, 49, 50, 51, 138, 127, 131, 173]
prob_numbers = filter(x -> !(x in issued_probs), all_probs)

algo_names = ["Algo1", "Algo2", "Algo3"]
cons_handle = "PB"
log_scaling = true

generate_files = false
draw_profiles = false
referee_please = false
typeof_referee = "All_Final"
draw_profiles_adjusted = false
confirm_profiles = true

@assert draw_profiles_adjusted <= referee_please "Cannot adjust the profiles if no referee."
path_jld2 = "/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/JLD2saves"

if generate_files
    generate_files!(prob_numbers, algo_names; cons_handle = cons_handle, path = path_jld2)
end

cd(path_jld2)
N_all_hists = load_object("N_all_hists-cons=$cons_handle.jld2")
F_all_hists = load_object("F_all_hists-cons=$cons_handle.jld2")
f_all_hists = load_object("f_all_hists-cons=$cons_handle.jld2")
x_all_hists = load_object("x_all_hists-cons=$cons_handle.jld2")
y_all_hists = load_object("y_all_hists-cons=$cons_handle.jld2")
cd("/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark")

αs = collect(1:100)
ks = collect(1:100)
y_perf = zeros(Float64, length(αs), length(algo_names))
y_data = zeros(Float64, length(ks), length(algo_names))

if confirm_profiles
    for τ in [1e-1, 1e-2]
        NapMatrix = Nap_Matrix(F_all_hists, N_all_hists, algo_names, prob_numbers, τ)
        perf_prof = performance_profile(PlotsBackend(), NapMatrix, algo_names, title="Performance Profile τ = $(τ*100)%")
        NpVector = ScaleDataDim(prob_numbers)
        data_prof = data_profile(PlotsBackend(), NapMatrix, NpVector, algo_names, title="Data Profile τ = $(τ*100)%")

        display(perf_prof)
        display(data_prof)
    end
end

if draw_profiles
    draw_profiles!([1e-1, 1e-2], αs, ks, algo_names, prob_numbers, F_all_hists, N_all_hists; cons_handle = cons_handle, log_scaling = log_scaling)
end

if referee_please
    if typeof_referee == "Single_Final"
        referee_name = "NOMAD"
        algo_blames = Referee_adjust(algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, referee_name)
        F_all_hists_adjusted = copy(F_all_hists)

        if draw_profiles_adjusted
            # Change F historics for the accuracy computation
            for algo in algo_names
                for issued_prob in algo_blames[algo]
                    prob_index = findall(x->x==issued_prob, prob_numbers)

                    # If the referee found a better final solution, invalidate ALL the historic of this algo for this problem
                    F_all_hists_adjusted[prob_index][1][algo] .= fill(Inf, length(F_all_hists[prob_index][1][algo]))
                end
            end
            draw_profiles!([1e-1, 1e-2], αs, ks, algo_names, prob_numbers, F_all_hists_adjusted, N_all_hists; adjusted = draw_profiles_adjusted, cons_handle = cons_handle, log_scaling = log_scaling)
        end
    elseif typeof_referee == "All_Final"
        algo_blames = Referee_all_adjust(algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, algo_names)
        F_all_hists_adjusted = copy(F_all_hists)

        if draw_profiles_adjusted
            # Change F historics for the accuracy computation
            for algo in algo_names
                for issued_prob in algo_blames[algo]
                    prob_index = findall(x->x==issued_prob, prob_numbers)

                    # If the referee found a better final solution, invalidate ALL the historic of this algo for this problem
                    F_all_hists_adjusted[prob_index][1][algo] .= fill(Inf, length(F_all_hists[prob_index][1][algo]))
                end
            end
            draw_profiles!([1e-1, 1e-2], αs, ks, algo_names, prob_numbers, F_all_hists_adjusted, N_all_hists; adjusted = draw_profiles_adjusted, cons_handle = cons_handle, log_scaling = log_scaling, type_of_ref = typeof_referee)
        end
    end
end