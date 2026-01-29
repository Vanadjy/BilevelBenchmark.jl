using BilevelBenchmark
using PGFPlots, JLD2, NOMAD
using Logging

# To check the profiles
using BenchmarkProfiles
using Plots

io = open("Output_intern_total_referee.txt", "w+")

include("plot_settings.jl")
include("plot-utils.jl")
include("GenerateFiles.jl")
include("DrawProfiles.jl")
include("Referee_adjust.jl")

#test_prob = 2
all_probs = collect(1:173)
issued_probs = [36, 49, 50, 51, 138, 127, 131, 173]
prob_numbers = filter(x -> !(x in issued_probs), all_probs)
conv_problems = [59, 82, 124]
#prob_numbers= [test_prob]

## displays all the problems with an infeasible starting point

#=infeasible_starting = []
for prob in prob_numbers
    model = get_bilevel_problem(prob)
    xk = model.xy0[1:model.dim[1]]
    yk = model.xy0[model.dim[1] + 1:model.dim[1] + model.dim[2]]
    if (any(>(0), model.g_func(xk, yk)))
        push!(infeasible_starting, prob)
    end
end=#

algo_names = ["Algo1", "Algo2", "Algo3"]

cons_handle = "EB"
log_scaling = true
starter = "y0" # "y0" or "yk-1"

upper_budget = 300
lower_budget = 100

bilevel_options = BilevelOptions(;
        subsolver_name = "NOMAD",
        γ = 1/2,
        oppportunistic = true,
        ordered = false,
        search = false,
        orthogonal = true,
        max_neval_upper = upper_budget,
        max_neval_upper_cons = upper_budget,
        max_neval_lower = lower_budget,
        Δ0 = 1.0,
        tol_upper = 1e-6,
        tol_lower = 1e-6,
        max_time = 3600.0,
        biphase = false,
        verbose = false
    )


hub_options = HubOPtions(
    generate_files = false,
    draw_conv = false,
    draw_profiles = true,
    referee_please = true,
    typeof_referee = "All_All", # Possibles : "Single_Final", "All_Final", "All_All", "All_Backward"
    generate_F_adjusted = true,
    draw_conv_adjusted = false,
    draw_profiles_adjusted = true,
    confirm_profiles = false,
    save_logs = false,
    ω_toggle = true
)

if hub_options.save_logs
    logger = SimpleLogger(io)
    global_logger(logger)
end

path_jld2 = "/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/JLD2saves"

#= # Open the file for writing (or appending with "a")
io = open("logfile_all_probs.txt", "w")

# Create a logger that writes to the file
logger = ConsoleLogger(io)

# Temporarily set the logger for the current scope
with_logger(logger) do

end =#

if hub_options.generate_files
    generate_files!(prob_numbers, algo_names, bilevel_options; cons_handle = cons_handle, path = path_jld2, start_point = starter, save = true)
end

# Close the file after logging is complete
#close(io)

## Load already generated historics ##

cd(path_jld2)
N_UL_all_hists = load_object("N_UL_all_hists-cons=$cons_handle-start=$starter-budg_u=$(upper_budget).jld2")
N_LL_all_hists = load_object("N_LL_all_hists-cons=$cons_handle-start=$starter-budg_u=$(upper_budget).jld2")
F_all_hists = load_object("F_all_hists-cons=$cons_handle-start=$starter-budg_u=$(upper_budget).jld2")
f_all_hists = load_object("f_all_hists-cons=$cons_handle-start=$starter-budg_u=$(upper_budget).jld2")
x_all_hists = load_object("x_all_hists-cons=$cons_handle-start=$starter-budg_u=$(upper_budget).jld2")
y_all_hists = load_object("y_all_hists-cons=$cons_handle-start=$starter-budg_u=$(upper_budget).jld2")
t_all_hists = load_object("t_all_hists-cons=$cons_handle-start=$starter-budg_u=$(upper_budget).jld2")
cd("/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark")


## Scale the N using λ ##


N_all_hists = copy(N_UL_all_hists)
λ_choice = "LL" # "UL" or "LL"
if hub_options.ω_toggle
    λ_list = load_object("numerics/time_and_omega_list/omega_list.jld2")
    X = collect(1:length(prob_numbers))
    gr()
    plot = scatter(X, λ_list, xlabel = "Problem index", ylabel = "λ value", title = "Value of λ for each problem", yaxis = :log)
    savefig(plot, "lambdas_scatter.pdf")
    for prob in eachindex(prob_numbers)
        for algo in keys(N_UL_all_hists[prob])
            if λ_choice == "LL"
                N_all_hists[prob][algo] .= λ_list[prob]*N_UL_all_hists[prob][algo] .+ N_LL_all_hists[prob][algo]
            else
                N_all_hists[prob][algo] .= N_UL_all_hists[prob][algo] .+ N_LL_all_hists[prob][algo]./λ_list[prob]
            end
        end
    end
end

## Drawing profiles, convergence plots and, if asked, confirming profiles ##

lim = 1000
αs = collect(1:1:lim)
ks = collect(1:1:lim)
ds = collect(1:1:lim)
y_perf = zeros(Float64, length(αs), length(algo_names))
y_data = zeros(Float64, length(ks), length(algo_names))
y_acc = zeros(Float64, length(ds), length(algo_names))

if hub_options.confirm_profiles
    for τ in [1e-1, 1e-2]
        NapMatrix = Nap_Matrix(F_all_hists, N_all_hists, algo_names, prob_numbers, τ)
        perf_prof = performance_profile(PlotsBackend(), NapMatrix, algo_names, title="Performance Profile τ = $(τ*100)%";) #ylims=(0.35,0.45))

        display(perf_prof)
    end
end

if hub_options.draw_conv
    for p in conv_problems
        draw_convergence!(F_all_hists, N_all_hists, p, algo_names; logscale = log_scaling, type_of_ref = "", cons_handle = cons_handle)
    end
end

if hub_options.draw_profiles
    draw_profiles!([1e-1, 1e-2], αs, ks, algo_names, prob_numbers, F_all_hists, N_all_hists; cons_handle = cons_handle, log_scaling = log_scaling, type_of_ref = "", start_point = starter, ω_toggle = hub_options.ω_toggle, λ_choice = λ_choice)
end

## Referee adjustments ##

if hub_options.referee_please
    if hub_options.typeof_referee == "Single_Final"
        referee_name = "NOMAD"
        algo_blames = Referee_adjust(algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, referee_name)
        F_all_hists_adjusted = copy(F_all_hists)
        N_all_hists_adjusted = copy(N_all_hists)

        if hub_options.draw_profiles_adjusted
            # Change F historics for the accuracy computation
            for algo in algo_names
                for issued_prob in algo_blames[algo]
                    prob_index = findfirst(x->x==issued_prob, prob_numbers)

                    # If the referee found a better final solution, invalidate ALL the historic of this algo for this problem
                    F_all_hists_adjusted[prob_index][algo] .= fill(Inf, length(F_all_hists[prob_index][algo]))
                    N_all_hists_adjusted[prob_index][algo] .= fill(Inf, length(N_all_hists[prob_index][algo]))
                end
            end
            draw_profiles!([1e-1, 1e-2], αs, ks, algo_names, prob_numbers, F_all_hists_adjusted, N_all_hists_adjusted; adjusted = draw_profiles_adjusted, cons_handle = cons_handle, log_scaling = log_scaling, start_point = starter, ω_toggle = hub_options.ω_toggle, λ_choice = λ_choice)
        end

    # Intern Final Referee #
    elseif hub_options.typeof_referee == "All_Final"
        algo_blames = Referee_all_adjust(algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, algo_names)
        F_all_hists_adjusted = copy(F_all_hists)
        N_all_hists_adjusted = copy(N_all_hists)

        if hub_options.draw_profiles_adjusted
            # Change F historics for the accuracy computation
            for algo in algo_names
                for issued_prob in algo_blames[algo]
                    prob_index = findfirst(x->x==issued_prob, prob_numbers)

                    # If the referee found a better final solution, invalidate ALL the historic of this algo for this problem
                    F_all_hists_adjusted[prob_index][algo] = []
                    N_all_hists_adjusted[prob_index][algo] = []
                end
            end
            draw_profiles!([1e-1, 1e-2], αs, ks, algo_names, prob_numbers, F_all_hists_adjusted, N_all_hists_adjusted; adjusted = hub_options.draw_profiles_adjusted, cons_handle = cons_handle, log_scaling = log_scaling, type_of_ref = hub_options.typeof_referee, start_point = starter, ω_toggle = hub_options.ω_toggle, λ_choice = λ_choice)
        end

    # Intern Total Referee #
    elseif hub_options.typeof_referee == "All_All"
        F_all_hists_adjusted = copy(F_all_hists)
        N_all_hists_adjusted = copy(N_all_hists)
        if hub_options.generate_F_adjusted
            Referee_all_historic_adjust!(F_all_hists_adjusted, N_all_hists_adjusted, algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, algo_names)
            cd(path_jld2)
            JLD2.save_object("F_all_hists_adjusted-cons=$cons_handle-start=$starter-$(hub_options.typeof_referee).jld2", F_all_hists_adjusted)
            cd("/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark")
        end
        cd(path_jld2)
        F_all_hists_adjusted = JLD2.load_object("F_all_hists_adjusted-cons=$cons_handle-start=$starter-$(hub_options.typeof_referee).jld2")
        cd("/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark")

        if hub_options.draw_profiles_adjusted
            draw_profiles!([1e-1, 1e-2], αs, ks, algo_names, prob_numbers, F_all_hists_adjusted, N_all_hists; adjusted = hub_options.draw_profiles_adjusted, cons_handle = cons_handle, log_scaling = log_scaling, type_of_ref = hub_options.typeof_referee, start_point = starter, ω_toggle = hub_options.ω_toggle, λ_choice = λ_choice)
        end
        if hub_options.draw_conv_adjusted
            for p in conv_problems
                draw_convergence!(F_all_hists_adjusted, N_all_hists, p, algo_names; logscale = log_scaling, type_of_ref = "adjusted", cons_handle = cons_handle, adjusted = hub_options.draw_conv_adjusted, start_point = starter)
            end
        end
    
    # Intern Backward Referee #
    elseif hub_options.typeof_referee == "All_Backward"
        F_all_hists_adjusted = copy(F_all_hists)
        N_all_hists_adjusted = copy(N_all_hists)
        if hub_options.generate_F_adjusted
            Referee_backward_historic_adjust!(F_all_hists_adjusted, N_all_hists_adjusted, algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, algo_names)
            cd(path_jld2)
            JLD2.save_object("F_all_hists_adjusted-backward-cons=$cons_handle-start=$starter-$(hub_options.typeof_referee).jld2", F_all_hists_adjusted)
            cd("/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark")
        end
        cd(path_jld2)
        F_all_hists_adjusted = JLD2.load_object("F_all_hists_adjusted-backward-cons=$cons_handle-start=$starter-$(hub_options.typeof_referee).jld2")
        cd("/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark")

        if hub_options.draw_profiles_adjusted
            draw_profiles!([1e-1, 1e-2], αs, ks, algo_names, prob_numbers, F_all_hists_adjusted, N_all_hists; adjusted = hub_options.draw_profiles_adjusted, cons_handle = cons_handle, log_scaling = log_scaling, type_of_ref = hub_options.typeof_referee, start_point = starter, ω_toggle = hub_options.ω_toggle, λ_choice = λ_choice)
        end
        if hub_options.draw_conv_adjusted
            for p in conv_problems
                draw_convergence!(F_all_hists_adjusted, N_all_hists, p, algo_names; logscale = log_scaling, type_of_ref = "adjusted", cons_handle = cons_handle, adjusted = hub_options.draw_conv_adjusted, start_point = starter)
            end
        end
    end
end

close(io)