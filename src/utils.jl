function Householder!(H::Matrix{Float64}, v::Vector{Float64})
    @assert dot(v, v) == 1 "The input vector must be a unit vector."
    H .= I - 2 * (v * v')
    return H
end