function conv_plot(p::Int; logscale = false)
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
end

accuracy(f, k::Int, p::Int, a::Int) = ((f[a](p)[k] - f[a](p)[1])/(f_star(p) - f[a](p)[1]))

function Nap(f, N, a::Int, p::Int, τ::Real)
    Nap = 10^16
    Tap = false
    i = 1
    while !(Tap || i >= length(N[1]))
        i +=1
        if accuracy(f, i, p, a) ≥ 1 - τ
            Nap = N[a][i]
            Tap = true
        end
    end
    return (Int(Nap), Tap)
end

function rap(f, N, a::Int, p::Int, τ::Real)
    tuple = Nap(f, N, a, p, τ)
    min_Nap = tuple[1]
    Tap = tuple[2]
    rap = 1e10
    if Tap
        @inbounds for algo in 1:3
            tuple_algo = Nap(f, N, algo, p, τ)
            if tuple_algo[2] && (tuple_algo[1] < min_Nap)
                min_Nap = tuple_algo[1]
            end
        end
        rap = tuple[1] / min_Nap
    end
    return rap
end

function perf_profile(αs, a::Int, τ::Real)
    count = 0
    y = []
    @inbounds for α in αs
        @inbounds for p in P
            if (rap(f, N, a, p, τ) ≤ α)
                count += 1
            end
        end
        ρ = count / (card_P)
        push!(y, ρ)
        count = 0
    end
    plot!(αs, y, linetype=:steppre, label = "Algorithm $a")
    xlabel!("Ratio of number of evaluations")
    ylabel!("Proportion of τ-solved problems")
end

function data_profile(ks, a::Int, τ::Real)
    count = 0
    y = []
    @inbounds for k in ks
        @inbounds for p in P
            tuple = Nap(f, N, a, p, τ)
            Nap_data, Tap_data = tuple[1], tuple[2]
            if Nap_data ≤ k * (p + 11) * Tap_data
                count += 1
            end
        end
        dk = count / (card_P)
        push!(y, dk)
        count = 0
    end
    plot!(ks, y, linetype=:steppre, label = "Algorithm $a")
    xlabel!("Groups of (p+11) evaluations")
    ylabel!("Proportion of τ-solved problems")
end

function accuracy_profile(ds, a::Int)
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
end


for p in [1, 10]
    conv_plot(p; logscale = true) #on remarque que f(x*) = 0.0
end
x = collect(1:20)

for τ in [1e-1, 1e-2]
    local graph = plot()
    for algo in 1:3
        title!("Performance profile for τ = $τ")
        perf_profile(x, algo, τ)
    end
    display(graph)
end

for τ in [1e-1, 1e-2]
    local graph = plot()
    for algo in 1:3
        title!("Data profile for τ = $τ")
        data_profile(x, algo, τ)
    end
    display(graph)
end

graph = plot()
d = collect(0:0.5:8)
for algo in 1:3
    title!("Accuracy profile")
    accuracy_profile(d, algo)
end
display(graph)