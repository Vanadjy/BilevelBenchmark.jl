export conv_plot, accuracy, Nap, rap, perf_profile!, data_profile!, accuracy_profile!

using LaTeXStrings, Plots

function conv_plot(f_hists, N_hists, prob::Int; logscale::Bool = false)
    graph = plot()
    ns = length(keys(N_hists[prob]))
    @inbounds for algo in keys(N_hists[prob])
        if logscale
            plot!(N_hists[prob][algo], f[prob][algo], linetype=:steppre, xaxis=:log10, yaxis=:log10, label=key(N_hists[prob])[i],
                  xlabel="Number of F evaluations",
                  ylabel="F value")
        else
            plot!(N[i], f[i](p), linetype=:steppre, label="algorithm $i")
            xlabel!("Number of F evaluations")
            ylabel!("F value")
        end
        title!("convergence plot for problem $p")
    end
    display(graph)
end

function f_star(f_hists, prob::Int)
    obj_hists = f_hists[prob]
    algos = collect(keys(obj_hists))
    n_algos = length(algos)

    @assert n_algos > 1 "Trying to compare only one algorithm for performace/data profiles"

    best_val = length(obj_hists[algos[1]]) == 0 ? Inf : obj_hists[algos[1]][end]
    for a in 2:n_algos
        candidate_val = length(obj_hists[algos[a]]) == 0 ? Inf : obj_hists[algos[a]][end]
        if best_val > candidate_val
            best_val = candidate_val
        end
    end
    return best_val
end

function f_0(f_hists, prob::Int, algo::Union{Int, String})
    f0 = length(f_hists[prob][algo]) == 0 ? Inf : f_hists[prob][algo][1]
    return f0
end

function f_0(f_hists, prob::Int)
    algos = collect(keys(f_hists[1]))
    n_algos = length(algos)
    @assert n_algos > 1 "Trying to compare only one algorithm for performace/data profiles"

    f0 = length(f_hists[1]) == 0 ? Inf : f_hists[1][1]
    for a in 2:n_algos
        # Routine to select the highest fist feasible f0
        candidate_f0 = length(f_hists[prob][algos[a]]) == 0 ? Inf : f_hists[prob][algos[a]][1]
        if f0 < candidate_f0
            f0 = candidate_f0
        end
    end  
    return f0
end


function accuracy(f_hists, k::Int, prob::Union{Int, String}, algo::Union{Int, String})
    f_N = length(f_hists[prob][algo]) == 0 ? Inf : f_hists[prob][algo][k]
    return ((f_N - f_0(f_hists, prob))/(f_star(f_hists, prob) - f_0(f_hists, prob)))
end

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
    return y
end

function data_profile!(y, ks, f_hist, N_hist, prob_list::Vector{Int}, algo::Union{Int, String}, τ::Real; ω_toggle::Bool = false, λ_choice::String = "LL")
    count = 0
    @inbounds for l in eachindex(ks)
        k = ks[l]
        @inbounds for prob in eachindex(prob_list)
            Nap_data, Tap_data = Nap(f_hist, N_hist, algo, prob, τ)
            model = get_bilevel_problem(prob_list[prob])
            if ω_toggle # if we scaled the UL evaluations with ω
                dimprob = λ_choice == "LL" ? model.dim[2] + 1 : model.dim[1] + 1
            else
                dimprob = model.dim[2]*(model.dim[1] +1)
            end
            if Nap_data ≤ k * (dimprob) * Tap_data
                count += 1
            end
        end
        dk = count / (length(prob_list))
        y[l] = dk
        count = 0
    end
    return y
end

function accuracy_profile!(y, ds, f_hist, prob_list::Vector{Int}, algo::Union{Int, String})
    count = 0
    @inbounds for i in eachindex(ds)
        d = ds[i]
        @inbounds for prob in eachindex(prob_list)
            k = length(f_hist[prob][algo])
            f_acc_tot = accuracy(f_hist, k, prob, algo)
            if isinf(f_acc_tot)
                f_acc_tot = 0.0
            end
            if  -log10(1 - f_acc_tot) ≥ d
                count += 1
            end
        end
        ratio = count / (length(prob_list))
        y[i] = ratio
        count = 0
    end
    return y
end