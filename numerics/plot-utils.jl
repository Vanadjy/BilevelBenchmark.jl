mutable struct HubOPtions{B, S}
    generate_files::B
    draw_conv::B
    draw_profiles::B
    referee_please::B
    typeof_referee::S
    generate_F_adjusted::B
    draw_conv_adjusted::B
    draw_profiles_adjusted::B
    confirm_profiles::B

    function HubOPtions{B, S}(;
                    generate_files::B = true,
                    draw_conv::B = true,
                    draw_profiles::B = false,
                    referee_please::B = false,
                    typeof_referee::S = "All_All",
                    generate_F_adjusted::B = false,
                    draw_conv_adjusted::B = false,
                    draw_profiles_adjusted::B = false,
                    confirm_profiles::B = false
        ) where {B <: Bool, S <: String}

        @assert draw_profiles_adjusted <= referee_please "Cannot adjust the profiles if no referee."
        @assert generate_F_adjusted <= referee_please "Cannot generate F adjusted if no referee."
        @assert draw_conv_adjusted <= referee_please "Cannot draw adjusted convergence plot if no referee."
        @assert typeof_referee ∈ ["All_All", "All_Final", "Single_Final"] "Invalid type of referee. Must be one of: All_All, All_Final, Single_Final."
        typeof_referee = (referee_please ? typeof_referee : "")

        return new{B, S}(
                generate_files,
                draw_conv,
                draw_profiles,
                referee_please,
                typeof_referee,
                generate_F_adjusted,
                draw_conv_adjusted,
                draw_profiles_adjusted,
                confirm_profiles
        )
    end
end

HubOPtions(args...; kwargs...) = HubOPtions{Bool, String}(args...; kwargs...)

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

function DataMatrix(F_all_hist, N_all_hist, algo_names::Vector{String})
    Data = Inf .* ones(eltype(F_all_hist[1]["Algo1"][1]), Int(N_all_hist[1]["Algo1"][end]), length(F_all_hist), length(keys(F_all_hists[1])))
    for p in eachindex(F_all_hist)
        for a in 1:length(keys(F_all_hists[1]))
            Data[Int.(N_all_hist[p][algo_names[a]]), p, a] .= F_all_hist[p][algo_names[a]]
        end
    end
    return Data
end

function ScaleDataDim(all_probs::Vector{Int})
    NpVector = zeros(Int, length(all_probs))
    for p in eachindex(all_probs)
        model = get_bilevel_problem(all_probs[p])
        NpVector[p] = model.dim[1] + model.dim[2] + 1
    end
    return NpVector
end