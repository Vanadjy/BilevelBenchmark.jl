function log_scale(n)
  try
    Int(log10(n))
  catch
    error("Input Error: n should be a power of 10")
  end
  log_scale = [k * 10.0^(i) for i in 0:Int(log10(n) - 1) for k in 1.0:9.0]
  return [Int.(log_scale); Int(10.0^(log10(n)))]
end

function Nap_Matrix(f_hists, N_hists, algo_names::Vector{String}, all_probs::Vector{Int}, τ::Real)
    na = length(algo_names)
    np = length(all_probs)

    NapMatrix = zeros(Float64, np, na)
    for a in eachindex(algo_names)
        for p in eachindex(all_probs)
            NapMatrix[p, a] = Nap(f_hists, N_hists, algo_names[a], p, τ)[1]
        end
    end
    return NapMatrix
end

function DataMatrix(F_all_hist, N_all_hist)
    Data = zeros(eltype(F_all_hist[1]["Algo1"][1]), )

end

function ScaleDataDim(all_probs::Vector{Int})
    NpVector = zeros(Int, length(all_probs))
    for p in eachindex(all_probs)
        model = get_bilevel_problem(all_probs[p])
        NpVector[p] = model.dim[1] + model.dim[2] + 1
    end
    return NpVector
end