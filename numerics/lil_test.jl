using JLD2
using CSV, DataFrames

algo_names = ["Algo1", "Algo2", "Algo3", "Algo4", "Algo5"]

all_probs = collect(1:173)
issued_probs = [36, 49, 50, 51, 138, 127, 131, 173]
prob_numbers = filter(x -> !(x in issued_probs), all_probs)
ignored_pbs = JLD2.load_object(joinpath("numerics/logs/", "ignored_pbs.jld2"))
kept_pbs = setdiff(collect(1:length(prob_numbers)), ignored_pbs)
prob_numbers = prob_numbers[kept_pbs] # Keeps only the problems that are refereed

path_jld2 = "/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/JLD2saves"

F_all_hists = load_object(joinpath(path_jld2, "F_all_hists-cons=EB-start=y0-budg_u=300.jld2"))
F_all_hists_endpoint = JLD2.load_object(joinpath(path_jld2, "F_all_hists_adjusted-endpoint-cons=EB-start=y0-All_Final.jld2"))
display(F_all_hists_endpoint[1])
F_all_hists_complete = JLD2.load_object(joinpath(path_jld2, "F_all_hists_adjusted-complete-cons=EB-start=y0-All_All.jld2"))
F_all_hists_reverse = JLD2.load_object(joinpath(path_jld2, "F_all_hists_adjusted-backward-cons=EB-start=y0-All_Backward.jld2"))

not_solved = Dict{String, Vector{Int}}(algo => Int[] for algo in algo_names)
solved_by_no_one = []

for p in eachindex(prob_numbers)
    max_length = maximum([length(F_all_hists[p][a]) for a in algo_names])
    full_matrix = zeros(max_length, 4*length(algo_names)) # 4 columns for the 4 different data sets
    problem_p_not_solved = []
    for a in algo_names
        #println(F_all_hists[p][a][1])
        #println(F_all_hists[p][a][end])
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
println("Problems ", solved_by_no_one, " are solved by no algorithm.")
for a in algo_names
    println("Algorithm $(a) did not solve ", length(not_solved[a])/length(prob_numbers)*100, "% of problems: ", not_solved[a])
end