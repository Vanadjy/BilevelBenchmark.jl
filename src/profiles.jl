export conv_plot, accuracy, Nap, rap, perf_profile!, data_profile!#, accuracy_profile

using LaTeXStrings, Plots

#=function conv_plot(p::Int; logscale = false)
    graph = plot()
    @inbounds for i in 1:3
        if logscale
            plot!(N[i][2:end], f[i](p)[2:end], linetype=:steppre, xaxis=:log10, yaxis=:log10, label="algorithm $i")
            xlabel!(L"$N$ in $log_{10}$ scale")
            ylabel!(L"$f(x^N)$ in $log_{10}$ scale")
        else
            plot!(N[i], f[i](p), linetype=:steppre, label="algorithm $i")
            xlabel!("Objective evaluations")
            ylabel!("Objective value")
        end
        title!("convergence plot for problem $p")
    end
    display(graph)
end=#

function f_star(f_hists, prob::Int)
    obj_hists = f_hists[prob]
    algos = collect(keys(obj_hists))
    n_algos = length(algos)

    @assert n_algos > 1 "Trying to compare only one algorithm for performace/data profiles"

    best_val = obj_hists[algos[1]][end]
    for a in 2:n_algos
        if best_val > obj_hists[algos[a]][end]
            best_val = obj_hists[algos[a]][end]
        end
    end
    return best_val
end

accuracy(f_hists, k::Int, prob::Union{Int, String}, algo::Union{Int, String}) = ((f_hists[prob][algo][k] - f_hists[prob][algo][1])/(f_star(f_hists, prob) - f_hists[prob][algo][1]))

function Nap(f_hists, N_hists, algo::Union{Int, String}, prob::Int, τ::Real)
    Nap = Inf
    Tap = false
    i = 1
    while !(Tap || i >= length(N_hists[prob][algo]))
        i += 1
        if accuracy(f_hists, i, prob, algo) ≥ 1 - τ
            Nap = N_hists[prob][algo][i]
            Tap = true
        end
    end
    return Nap, Tap
end

function rap(f_hists, N_hists, algo::Union{Int, String}, prob::Union{Int, String}, τ::Real)
    Nap_ref, Tap_ref = Nap(f_hists, N_hists, algo, prob, τ)
    champ_Nap = Nap_ref
    rap = Inf
    if Tap_ref
        all_algos = collect(keys(f_hists[prob]))
        for alg in all_algos
            Nap_algo, Tap_algo = Nap(f_hists, N_hists, alg, prob, τ)
            if Tap_algo && (Nap_algo < champ_Nap)
                champ_Nap = Nap_algo
            end
        end
        rap = Nap_ref / champ_Nap
    end
    return rap
end

function perf_profile!(y, αs, f_hist, N_hist, prob_list::Vector{Int}, algo::Union{Int, String}, τ::Real)
    count = 0
    @inbounds for l in eachindex(αs)
        α = αs[l]
        @inbounds for prob in eachindex(prob_list)
            if (rap(f_hist, N_hist, algo, prob, τ) ≤ α)
                count += 1
            end
        end
        ρ = count / (length(prob_list))
        y[l] = ρ
        count = 0
    end
    #=plot!(αs, y, linetype=:steppre, label = "Algorithm $algo")
    xlabel!("Ratio of number of evaluations")
    ylabel!("Proportion of τ-solved problems")=#
end

function data_profile!(y, ks, f_hist, N_hist, prob_list::Vector{Int}, algo::Union{Int, String}, τ::Real)
    count = 0
    @inbounds for l in eachindex(ks)
        k = ks[l]
        @inbounds for prob in eachindex(prob_list)
            Nap_data, Tap_data = Nap(f_hist, N_hist, algo, prob, τ)
            model = get_bilevel_problem(prob_list[prob])
            dimprob = model.dim[1] + model.dim[2]
            if Nap_data ≤ k * (dimprob + 1) * Tap_data
                count += 1
            end
        end
        dk = count / (length(prob_list))
        y[l] = dk
        count = 0
    end
    return y
    #=plot!(ks, y, linetype=:steppre, label = "Algorithm $a")
    xlabel!("Groups of (p+11) evaluations")
    ylabel!("Proportion of τ-solved problems")=#
end

#=function accuracy_profile(ds, a::Int)
    count = 0
    y = []
    k = length(N[a])
    @inbounds for d in ds
        @inbounds for p in P
            f_acc_tot = accuracy(f, k, p, a)
            if  -log10(1 - f_acc_tot) ≥ d
                count += 1
            end
        end
        ratio = count / (card_P)
        push!(y, ratio)
        count = 0
    end
    plot!(ds, y, linetype=:steppre, label = "Algorithm $a")
    ylabel!(L"Ratio of $f_{acc}^{N_{a,p}^{tot}}$ with more than d decimals")
    xlabel!("Number of decimals d")
end=#