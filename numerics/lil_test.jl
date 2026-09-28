using JLD2
using CSV, DataFrames
using BilevelBenchmark
using Statistics

include("plot_settings.jl")
include("plot-utils.jl")
include("GenerateFiles.jl")
include("DrawProfiles.jl")
include("Referee_adjust.jl")

purpose = "Anonymous_Referee"
Possible_algos = Dict("Anonymous_Referee" => Union{String, Int}["Algo1", "Algo2", "Algo5"],
                        "Known_Referee" => Union{String, Int}["Algo1", "Algo2", "Algo3"],
                        "Compare_Lambda" => Union{String, Int}["Algo1", "Algo5"])
Referees = Dict(
    "Known_Referee" => "Intern",
    "Anonymous_Referee" => "Extern"
)
algo_names = Possible_algos[purpose]
typeof_ref = Referees[purpose]

valid_referees = ["intern_endpoint", "intern_complete", "intern_reverse", "extern_endpoint", "extern_complete", "extern_reverse"]

all_probs = collect(1:173)
issued_probs = [36, 49, 50, 51, 138, 127, 131, 173]
prob_numbers = filter(x -> !(x in issued_probs), all_probs)
ignored_pbs = JLD2.load_object(joinpath("numerics/logs/", "ignored_pbs.jld2"))
kept_pbs = setdiff(collect(1:length(prob_numbers)), ignored_pbs)
prob_numbers = prob_numbers[kept_pbs] # Keeps only the problems that are refereed

path_jld2 = "/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/JLD2saves"

λ_list = load_object("numerics/logs/lambda_list.jld2")
λ_list = 60*ones(length(prob_numbers))
λ_choice = "UL" # "LL" or "UL"
effort_choice = "Agregate" # "UL", "LL" or "Agregate"
N_UL_all_hists = load_object(joinpath(path_jld2, "N_UL_all_hists-cons=EB-start=y0-budg_u=300.jld2"))
N_LL_all_hists = load_object(joinpath(path_jld2, "N_LL_all_hists-cons=EB-start=y0-budg_u=300.jld2"))
t_all_hists = load_object(joinpath(path_jld2, "t_all_hists-cons=EB-start=y0-budg_u=300.jld2"))
N_all_hists = copy(N_UL_all_hists)


for prob in eachindex(prob_numbers)
    for algo in keys(N_UL_all_hists[prob])
        if λ_choice == "LL"
            @assert effort_choice ∈ ["LL", "Agregate"] "Invalid effort choice. Must be one of: LL or Agregate when λ_choice is LL."
            N_all_hists[prob][algo] .= (λ_list[prob] .* N_UL_all_hists[prob][algo]) .+ N_LL_all_hists[prob][algo]
        else
            @assert effort_choice ∈ ["UL", "Agregate"] "Invalid effort choice. Must be one of: UL or Agregate when λ_choice is UL."
            N_all_hists[prob][algo] .= N_UL_all_hists[prob][algo] .+ (N_LL_all_hists[prob][algo]./λ_list[prob])
        end
    end
end
for k in [1, 2, 3, 4, 5]
    for prob in eachindex(prob_numbers)
        model = get_bilevel_problem(prob_numbers[prob])

        N_all_hists[prob]["Algo"*string(k)] .= N_all_hists[prob]["Algo"*string(k)]./(model.dim[2] + 1) # Normalization of the number of evaluations by the dimension of the problem, to be able to compare the problems with different dimensions on the same plot
        N_LL_all_hists[prob]["Algo"*string(k)] .= N_LL_all_hists[prob]["Algo"*string(k)]./(model.dim[2] + 1) # Normalization of the number of evaluations by the total number of evaluations (UL + LL), to be able to compare the problems with different budgets on the same plot
        local N_UL_all_hists = load_object(joinpath(path_jld2, "N_UL_all_hists-cons=EB-start=y0-budg_u=300.jld2"))
        N_UL_all_hists[prob]["Algo"*string(k)] .= N_UL_all_hists[prob]["Algo"*string(k)]./(model.dim[1] + 1) # Normalization of the number of evaluations by the total number of evaluations (UL + LL), to be able to compare the problems with different budgets on the same plot
    end
end
F_all_hists = load_object(joinpath(path_jld2, "F_all_hists-cons=EB-start=y0-budg_u=300.jld2"))
F_all_hists_endpoint = JLD2.load_object(joinpath(path_jld2, "F_all_hists_adjusted-purpose=$(purpose)-$(typeof_ref)_EndPoint-lambda=true-effort=Agregate.jld2"))
F_all_hists_complete = JLD2.load_object(joinpath(path_jld2, "F_all_hists_adjusted-purpose=$(purpose)-$(typeof_ref)_Complete-lambda=true-effort=Agregate.jld2"))
F_all_hists_reverse = JLD2.load_object(joinpath(path_jld2, "F_all_hists_adjusted-purpose=$(purpose)-$(typeof_ref)_Reverse-lambda=true-effort=Agregate.jld2"))

global not_solved = Dict{String, Vector{Int}}(algo => Int[] for algo in algo_names)
global not_solved_endpoint = Dict{String, Vector{Int}}(algo => Int[] for algo in algo_names)
global not_solved_complete = Dict{String, Vector{Int}}(algo => Int[] for algo in algo_names)
global not_solved_reverse = Dict{String, Vector{Int}}(algo => Int[] for algo in algo_names)

solved_by_no_one = []
solved_by_no_one_endpoint = []
solved_by_no_one_complete = []
solved_by_no_one_reverse = []

all_dims_nx = [get_bilevel_problem(prob_numbers[p]).dim[1] for p in eachindex(prob_numbers)]
all_dims_ny = [get_bilevel_problem(prob_numbers[p]).dim[2] for p in eachindex(prob_numbers)]

for p in eachindex(prob_numbers)
    max_length = maximum([length(F_all_hists[p][a]) for a in algo_names])
    full_matrix = zeros(max_length, 4*length(algo_names)) # 4 columns for the 4 different data sets
    problem_p_not_solved = []
    for a in algo_names
        if F_all_hists[p][a][1] == F_all_hists[p][a][end] # Problem p not solved by algo a
            push!(not_solved[a], prob_numbers[p])
            push!(problem_p_not_solved, a)
            if length(problem_p_not_solved) == length(algo_names) # Problem p not solved by any algorithm
                push!(solved_by_no_one, prob_numbers[p])
            end
        end
        #=Mat_res_algo = hcat([F_all_hists[p][a]; fill(NaN, max_length-length(F_all_hists[p][a]))],
            [F_all_hists_endpoint[p][a]; fill(NaN, max_length-length(F_all_hists_endpoint[p][a]))], 
            [F_all_hists_complete[p][a]; fill(NaN, max_length-length(F_all_hists_complete[p][a]))], 
            [F_all_hists_reverse[p][a]; fill(NaN, max_length-length(F_all_hists_reverse[p][a]))]
        )
        full_matrix[:, (findfirst(==(a), algo_names)-1)*4 .+ (1:4)] .= Mat_res_algo=#
    end
    #full_df = DataFrame(full_matrix, Symbol.(vec([string(a, "_", suffix) for a in algo_names for suffix in ["hist-problem_$(prob_numbers[p])", "endpoint-problem_$(prob_numbers[p])", "complete-problem_$(prob_numbers[p])", "reverse-problem_$(prob_numbers[p])"]])))
    #CSV.write("problem_$(prob_numbers[p]).csv",full_df)
end
@info "Problems $(solved_by_no_one) are solved by no algorithm."

new_prob_numbers = setdiff(prob_numbers, solved_by_no_one)

for p in eachindex(new_prob_numbers)
    problem_p_not_solved_endpoint = []
    problem_p_not_solved_complete = []
    problem_p_not_solved_reverse = []
    for a in algo_names
        if length(F_all_hists_endpoint[p][a]) == 0  # Problem p not solved by algo a
            push!(not_solved_endpoint[a], prob_numbers[p])
            push!(problem_p_not_solved_endpoint, a)
            if length(problem_p_not_solved_endpoint) == length(algo_names) # Problem p not solved by any algorithm
                push!(solved_by_no_one_endpoint, prob_numbers[p])
            end
        elseif length(F_all_hists_endpoint[p][a]) > 0 &&F_all_hists_endpoint[p][a][1] == F_all_hists_endpoint[p][a][end]
            push!(not_solved_endpoint[a], prob_numbers[p])
            push!(problem_p_not_solved_endpoint, a)
            if length(problem_p_not_solved_endpoint) == length(algo_names) # Problem p not solved by any algorithm
                push!(solved_by_no_one_endpoint, prob_numbers[p])
            end
        end
        if length(F_all_hists_complete[p][a]) == 0 || F_all_hists_complete[p][a][1] == F_all_hists_complete[p][a][end] # Problem p not solved by algo a
            push!(not_solved_complete[a], prob_numbers[p])
            push!(problem_p_not_solved_complete, a)
            if length(problem_p_not_solved_complete) == length(algo_names) # Problem p not solved by any algorithm
                push!(solved_by_no_one_complete, prob_numbers[p])
            end
        elseif length(F_all_hists_complete[p][a]) > 0 && F_all_hists_complete[p][a][1] == F_all_hists_complete[p][a][end]
            push!(not_solved_complete[a], prob_numbers[p])
            push!(problem_p_not_solved_complete, a)
            if length(problem_p_not_solved_complete) == length(algo_names) # Problem p not solved by any algorithm
                push!(solved_by_no_one_complete, prob_numbers[p])
            end
        end
        if length(F_all_hists_reverse[p][a]) == 0 || F_all_hists_reverse[p][a][1] == F_all_hists_reverse[p][a][end] # Problem p not solved by algo a
            push!(not_solved_reverse[a], prob_numbers[p])
            push!(problem_p_not_solved_reverse, a)
            if length(problem_p_not_solved_reverse) == length(algo_names) # Problem p not solved by any algorithm
                push!(solved_by_no_one_reverse, prob_numbers[p])
            end
        elseif length(F_all_hists_reverse[p][a]) > 0 && F_all_hists_reverse[p][a][1] == F_all_hists_reverse[p][a][end]
            push!(not_solved_reverse[a], prob_numbers[p])
            push!(problem_p_not_solved_reverse, a)
            if length(problem_p_not_solved_reverse) == length(algo_names) # Problem p not solved by any algorithm
                push!(solved_by_no_one_reverse, prob_numbers[p])
            end
        end
    end
end

@info "The number of problems solved for $(typeof_ref)al Referee is: $(length(solved_by_no_one_endpoint)) out of $(length(prob_numbers))."
@info "The number of problems solved for $(typeof_ref)al Referee is: $(length(solved_by_no_one_complete)) out of $(length(prob_numbers))."
@info "The number of problems solved for $(typeof_ref)al Referee is: $(length(solved_by_no_one_reverse)) out of $(length(prob_numbers))."


for a in algo_names
    println("Algorithm $(a) did not solve ", length(not_solved[a])/length(prob_numbers)*100, "% of problems: ", not_solved[a])
end

@info "Total CPU time spent to complete the original optimization: $(sum([t_all_hists[p][a][end] for p in eachindex(prob_numbers), a in algo_names]))."
for referee in valid_referees
    time_referee = load_object(joinpath(path_jld2, "elapsed_time_$(referee)-tol_ref=1.0e-9.jld2"))
    @info "Total CPU time spent for the $referee: $(time_referee)"
end

## Additional data profiles

#=lim = 100
n_probs = 50
prob_numbers = collect(1:n_probs)

ds = collect(1:lim)
αs = similar(ds)
ks = similar(ds)

y_perf = zeros(Float64, length(αs), 2)
y_data = zeros(Float64, length(ks), 2)
y_acc = zeros(Float64, length(ds), 2)

## Generating imaginary histories
All_initial_hists_F = [Dict("A" => sort(0.5 .+ rand(50), rev = true), "B" => sort(rand(50), rev = true)) for i in 1:n_probs]
All_initial_hists_N = [Dict("A" => collect(range(1, 100, step=2)), "B" => collect(range(1, 100, step=2))) for i in 1:n_probs]

τ_values = [1e-2]
algo_names = Union{Int64, String}["A", "B"]

for τ in τ_values
    for a in eachindex(algo_names)
        #@views perf_profile!(y_perf[:, a], αs, All_initial_hists_F, N_all_hists, prob_numbers, algo_names[a], τ, algo_names)
        @views data_profile!(y_data[:, a], ks, All_initial_hists_F, All_initial_hists_N, prob_numbers, algo_names[a], τ, algo_names; λ_toggle = false, λ_choice = "UL",  budget_UL = 100, budget_LL = 10, effort_choice = "Agregate")
        #@views accuracy_profile!(y_acc[:, a], ds, All_initial_hists_F, prob_numbers, algo_names[a], algo_names; opt_known = false)
    end
end

display(y_data)=#