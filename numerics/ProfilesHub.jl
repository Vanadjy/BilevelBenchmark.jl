using BilevelBenchmark
using PGFPlots, JLD2, NOMAD
using Logging

# To check the profiles
using BenchmarkProfiles
using Plots
import Random

## ----------------- BASH ARGUMENTS ------------------------ ##
if length(ARGS) == 0
    log_scaling = false
    num_algos = 4
    referee = "Intern_EndPoint"
    @assert referee ∈ ["Intern_EndPoint", "Intern_Complete", "Intern_Reverse", "Extern_EndPoint", "Extern_Complete", "Extern_Reverse"] "Invalid referee type. Must be one of: Intern_EndPoint, Intern_Complete, Intern_Reverse, Extern_EndPoint, Extern_Complete, Extern_Reverse."
    algo_names = num_algos == 2 ? Union{String, Int}["Algo1", "Algo2"] : Union{String, Int}["Algo1", "Algo2", "Algo3", "Algo4", "Algo5"]
else
    log_scaling = ARGS[1] == "log" ? true : false
    referee = ARGS[2]
    num_algos = parse(Int, ARGS[3])
    @assert referee ∈ ["Intern_EndPoint", "Intern_Complete", "Intern_Reverse", "Extern_EndPoint", "Extern_Complete", "Extern_Reverse"] "Invalid referee type. Must be one of: Intern_EndPoint, Intern_Complete, Intern_Reverse, Extern_EndPoint, Extern_Complete, Extern_Reverse."
    algo_names = num_algos == 2 ? Union{String, Int}["Algo1", "Algo2"] : Union{String, Int}["Algo1", "Algo2", "Algo3", "Algo4", "Algo5"]
end

## ----------------- MAIN CODE ----------------------------- ##

io = open("Output_intern_total_referee.txt", "w+")

include("plot_settings.jl")
include("plot-utils.jl")
include("GenerateFiles.jl")
include("DrawProfiles.jl")
include("Referee_adjust.jl")

seed = 1234
Random.seed!(seed)

all_probs = collect(1:173)
issued_probs = [36, 49, 50, 51, 138, 127, 131, 173]
prob_numbers = filter(x -> !(x in issued_probs), all_probs)
ignored_pbs = JLD2.load_object(joinpath("numerics/logs/", "ignored_pbs.jld2"))
kept_pbs = setdiff(collect(1:length(prob_numbers)), ignored_pbs)
prob_numbers = prob_numbers[kept_pbs] # Keeps only the problems that are refereed
JLD2.save_object("numerics/logs/prob_numbers.jld2", prob_numbers)
#println(prob_numbers)
conv_problems = [10]
conv_problem_indexes = findall(in(conv_problems), prob_numbers)
@assert conv_problems ⊆ prob_numbers "Some convergence problems are not in the list of kept problems."

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

cons_handle = "EB"
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
    typeof_referee = referee, # Possibles : "Single_EndPoint", "Intern_EndPoint", "Intern_Complete", "Intern_Reverse"
    generate_F_adjusted = false, # If it is necessary to rerun the refereeing procedures to generate the adjusted F historics (can be long, especially for the Intern_Complete and Intern_Reverse)
    draw_conv_adjusted = false,
    draw_profiles_adjusted = true,
    confirm_profiles = false,
    save_logs = false,
    λ_toggle = true
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

N_UL_all_hists = load_object(joinpath(path_jld2, "N_UL_all_hists-cons=$cons_handle-start=$starter-budg_u=$(upper_budget).jld2"))
N_LL_all_hists = load_object(joinpath(path_jld2, "N_LL_all_hists-cons=$cons_handle-start=$starter-budg_u=$(upper_budget).jld2"))
F_all_hists = load_object(joinpath(path_jld2, "F_all_hists-cons=$cons_handle-start=$starter-budg_u=$(upper_budget).jld2"))
f_all_hists = load_object(joinpath(path_jld2, "f_all_hists-cons=$cons_handle-start=$starter-budg_u=$(upper_budget).jld2"))
x_all_hists = load_object(joinpath(path_jld2, "x_all_hists-cons=$cons_handle-start=$starter-budg_u=$(upper_budget).jld2"))
y_all_hists = load_object(joinpath(path_jld2, "y_all_hists-cons=$cons_handle-start=$starter-budg_u=$(upper_budget).jld2"))
t_all_hists = load_object(joinpath(path_jld2, "t_all_hists-cons=$cons_handle-start=$starter-budg_u=$(upper_budget).jld2"))

## ------------------- ##
## Scale the N using λ ##
## ------------------- ##

N_all_hists = copy(N_UL_all_hists)
λ_choice = "UL" # "UL" or "LL"
effort_choice = "UL" # "UL", "LL" or "Agregate"
if hub_options.λ_toggle
    λ_list = load_object("numerics/logs/lambda_list.jld2")
    X = collect(1:length(prob_numbers))
    #gr()
    #plot = scatter(X, λ_list, xlabel = "Problem index", ylabel = "λ value", title = "Value of λ for each problem", yaxis = :log)
    #savefig(plot, "lambdas_scatter.pdf")
    for prob in eachindex(prob_numbers)
        for algo in keys(N_UL_all_hists[prob])
            if λ_choice == "LL"
                @assert effort_choice ∈ ["LL", "Agregate"] "Invalid effort choice. Must be one of: LL or Agregate when λ_choice is LL."
                N_all_hists[prob][algo] .= N_UL_all_hists[prob][algo]./λ_list[prob] .+ N_LL_all_hists[prob][algo]
            else
                @assert effort_choice ∈ ["UL", "Agregate"] "Invalid effort choice. Must be one of: UL or Agregate when λ_choice is UL."
                N_all_hists[prob][algo] .= N_UL_all_hists[prob][algo] .+ (λ_list[prob] .* N_LL_all_hists[prob][algo])
            end
        end
    end
else
    for prob in eachindex(prob_numbers)
        for algo in keys(N_UL_all_hists[prob])
            @assert effort_choice ∈ ["UL", "LL"] "Invalid effort choice. Must be one of: UL or LL when upper and lower evaulations are not agregated."
            if effort_choice == "UL"
                N_all_hists[prob][algo] .= N_UL_all_hists[prob][algo]
            elseif effort_choice == "LL"
                N_all_hists[prob][algo] .= N_LL_all_hists[prob][algo]
            end
        end
    end
end

## ---------------------------------------------------------------------- ##
## Drawing profiles, convergence plots and, if asked, confirming profiles ##
## ---------------------------------------------------------------------- ##

τ_values = [1e-1, 1e-4, 1e-8]
lim = 1000

ds = collect(1:lim)
αs = similar(ds)
ks = similar(ds)
if log_scaling # αs and ks need to have the same last digit
    αs = collect(1:0.01:100)
    ks = vcat(collect(0:9), collect(10:10:90), collect(100:100:lim))
    ds = collect(0:0.01:10)
else
    αs = collect(1:0.01:10)
    ks = collect(0:1:300)
    ds = collect(0:0.01:10)
end
y_perf = zeros(Float64, length(αs), length(algo_names))
y_data = zeros(Float64, length(ks), length(algo_names))
y_acc = zeros(Float64, length(ds), length(algo_names))

if hub_options.confirm_profiles
    for τ in τ_values
        NapMatrix = Nap_Matrix(F_all_hists, N_all_hists, algo_names, prob_numbers, τ)
        perf_prof = performance_profile(PlotsBackend(), NapMatrix, ["Algo1", "Algo2"], title="Performance Profile τ = $(τ*100)%";) #ylims=(0.35,0.45))

        display(perf_prof)
    end
end

if hub_options.draw_conv
    for p in conv_problem_indexes
        draw_convergence!(F_all_hists, N_all_hists, p, algo_names; logscale = log_scaling, type_of_ref = "", cons_handle = cons_handle)
    end
end

if hub_options.draw_profiles
    draw_profiles!(τ_values, 
                   αs, 
                   ks, 
                   algo_names, 
                   prob_numbers, 
                   F_all_hists, 
                   N_all_hists; 
                   cons_handle = cons_handle, 
                   log_scaling = log_scaling, 
                   type_of_ref = "No Referee", 
                   start_point = starter, 
                   λ_toggle = hub_options.λ_toggle, 
                   λ_choice = λ_choice, 
                   budget_UL = upper_budget, 
                   budget_LL = lower_budget,
                   effort_choice = effort_choice
    )
end

if hub_options.referee_please
    if hub_options.typeof_referee == "Single_EndPoint"
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
            draw_profiles!(τ_values, 
                           αs, 
                           ks, 
                           algo_names, 
                           prob_numbers, 
                           F_all_hists_adjusted, 
                           N_all_hists_adjusted; 
                           adjusted = draw_profiles_adjusted, 
                           cons_handle = cons_handle, 
                           log_scaling = log_scaling, 
                           start_point = starter, 
                           λ_toggle = hub_options.λ_toggle, 
                           λ_choice = λ_choice, 
                           budget_UL = upper_budget, 
                           budget_LL = lower_budget,
                           effort_choice = effort_choice
            )
        end

    ## --------------------------------------- ##
    ## Referee adjustments : Internal approach ##
    ## --------------------------------------- ##

    # Intern End-Point Referee #
    elseif hub_options.typeof_referee == "Intern_EndPoint"

        if hub_options.generate_F_adjusted
            algo_blames = EndPoint_Referee(algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, algo_names; max_budget = max_neval_lower)
            F_all_hists_adjusted = copy(F_all_hists)
            N_all_hists_adjusted = copy(N_all_hists)
            # Change F historics for the accuracy computation
            for algo in algo_names
                for issued_prob in algo_blames[algo]
                    prob_index = findfirst(x->x==issued_prob, prob_numbers)

                    # If the referee found a better final solution, invalidate ALL the historic of this algo for this problem
                    F_all_hists_adjusted[prob_index][algo] = []
                    N_all_hists_adjusted[prob_index][algo] = []
                end
            end
            JLD2.save_object(joinpath(path_jld2,"F_all_hists_adjusted-intern-endpoint-cons=$cons_handle-start=$starter-$(hub_options.typeof_referee).jld2"), F_all_hists_adjusted)
            JLD2.save_object(joinpath(path_jld2,"N_all_hists_adjusted-intern-endpoint-cons=$cons_handle-start=$starter-$(hub_options.typeof_referee).jld2"), N_all_hists_adjusted)
        end
        F_all_hists_adjusted = JLD2.load_object(joinpath(path_jld2, "F_all_hists_adjusted-intern-endpoint-cons=$cons_handle-start=$starter-$(hub_options.typeof_referee).jld2"))
        N_all_hists_adjusted = JLD2.load_object(joinpath(path_jld2, "N_all_hists_adjusted-intern-endpoint-cons=$cons_handle-start=$starter-$(hub_options.typeof_referee).jld2"))
        if hub_options.draw_profiles_adjusted
            draw_profiles!(τ_values, 
                           αs, 
                           ks, 
                           algo_names, 
                           prob_numbers, 
                           F_all_hists_adjusted, 
                           N_all_hists_adjusted; 
                           adjusted = hub_options.draw_profiles_adjusted, 
                           cons_handle = cons_handle, 
                           log_scaling = log_scaling, 
                           type_of_ref = hub_options.typeof_referee, 
                           start_point = starter, 
                           λ_toggle = hub_options.λ_toggle, 
                           λ_choice = λ_choice, 
                           budget_UL = upper_budget, 
                           budget_LL = lower_budget,
                           effort_choice = effort_choice
            )
        end

    # Intern Complete Referee #
    elseif hub_options.typeof_referee == "Intern_Complete"
        if hub_options.generate_F_adjusted
            F_all_hists_adjusted = copy(F_all_hists)
            N_all_hists_adjusted = copy(N_all_hists)
            Complete_Referee!(F_all_hists_adjusted, N_all_hists_adjusted, algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, algo_names; max_budget = max_neval_lower)
            JLD2.save_object(joinpath(path_jld2,"F_all_hists_adjusted-intern-complete-cons=$cons_handle-start=$starter-$(hub_options.typeof_referee).jld2"), F_all_hists_adjusted)
            JLD2.save_object(joinpath(path_jld2,"N_all_hists_adjusted-intern-complete-cons=$cons_handle-start=$starter-$(hub_options.typeof_referee).jld2"), N_all_hists_adjusted)
        end
        
        F_all_hists_adjusted = JLD2.load_object(joinpath(path_jld2, "F_all_hists_adjusted-intern-complete-cons=$cons_handle-start=$starter-$(hub_options.typeof_referee).jld2"))
        N_all_hists_adjusted = JLD2.load_object(joinpath(path_jld2, "N_all_hists_adjusted-intern-complete-cons=$cons_handle-start=$starter-$(hub_options.typeof_referee).jld2"))

        if hub_options.draw_profiles_adjusted
            draw_profiles!(τ_values, 
                           αs, 
                           ks, 
                           algo_names, 
                           prob_numbers, 
                           F_all_hists_adjusted, 
                           N_all_hists_adjusted; 
                           adjusted = hub_options.draw_profiles_adjusted, 
                           cons_handle = cons_handle, 
                           log_scaling = log_scaling, 
                           type_of_ref = hub_options.typeof_referee, 
                           start_point = starter, 
                           λ_toggle = hub_options.λ_toggle, 
                           λ_choice = λ_choice, 
                           budget_UL = upper_budget, 
                           budget_LL = lower_budget,
                           effort_choice = effort_choice
            )
        end
        if hub_options.draw_conv_adjusted
            for p in conv_problem_indexes
                draw_convergence!(F_all_hists_adjusted, N_all_hists_adjusted, p, algo_names; logscale = log_scaling, type_of_ref = "All_All", cons_handle = cons_handle, adjusted = hub_options.draw_conv_adjusted, start_point = starter)
            end
        end
    
    # Intern Reverse Referee #
    elseif hub_options.typeof_referee == "Intern_Reverse"
        if hub_options.generate_F_adjusted
            F_all_hists_adjusted = copy(F_all_hists)
            N_all_hists_adjusted = copy(N_all_hists)
            Reverse_Referee!(F_all_hists_adjusted, N_all_hists_adjusted, algo_names, prob_numbers, x_all_hists, y_all_hists, f_all_hists, algo_names; max_budget = max_neval_lower)
            JLD2.save_object(joinpath(path_jld2, "F_all_hists_adjusted-intern-reverse-cons=$cons_handle-start=$starter-$(hub_options.typeof_referee).jld2"), F_all_hists_adjusted)
            JLD2.save_object(joinpath(path_jld2, "N_all_hists_adjusted-intern-reverse-cons=$cons_handle-start=$starter-$(hub_options.typeof_referee).jld2"), N_all_hists_adjusted)
        end

        F_all_hists_adjusted = JLD2.load_object(joinpath(path_jld2, "F_all_hists_adjusted-intern-reverse-cons=$cons_handle-start=$starter-$(hub_options.typeof_referee).jld2"))
        N_all_hists_adjusted = JLD2.load_object(joinpath(path_jld2, "N_all_hists_adjusted-intern-reverse-cons=$cons_handle-start=$starter-$(hub_options.typeof_referee).jld2"))
        
        if hub_options.draw_profiles_adjusted
            draw_profiles!(τ_values, 
                           αs, 
                           ks, 
                           algo_names, 
                           prob_numbers, 
                           F_all_hists_adjusted, 
                           N_all_hists_adjusted; 
                           adjusted = hub_options.draw_profiles_adjusted, 
                           cons_handle = cons_handle, 
                           log_scaling = log_scaling, 
                           type_of_ref = hub_options.typeof_referee, 
                           start_point = starter, 
                           λ_toggle = hub_options.λ_toggle, 
                           λ_choice = λ_choice, 
                           budget_UL = upper_budget, 
                           budget_LL = lower_budget,
                           effort_choice = effort_choice
            )
        end
        if hub_options.draw_conv_adjusted
            for p in conv_problem_indexes
                draw_convergence!(F_all_hists_adjusted, N_all_hists_adjusted, p, algo_names; logscale = log_scaling, type_of_ref = "All_Reverse", cons_handle = cons_handle, adjusted = hub_options.draw_conv_adjusted, start_point = starter)
            end
        end

    ## --------------------------------------- ##
    ## Referee adjustments : External approach ##
    ## --------------------------------------- ##

    elseif hub_options.typeof_referee == "Extern_EndPoint"
        @warn "The Extern_EndPoint referee is not implemented yet. Use with caution and check the results."
    end
end



close(io)