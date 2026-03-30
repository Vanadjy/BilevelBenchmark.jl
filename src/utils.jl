export get_x0, get_y0, get_opt_val

function Householder!(H::Matrix{Float64}, v::Vector{Float64})
    @assert dot(v, v) ≈ 1 "The input vector must be a unit vector."
    H .= I - 2 * (v * v')
    return H
end

function get_x0(model::BilevelProblem)
    return model.xy0[1:model.dim[1]]
end

function get_y0(model::BilevelProblem)
    return model.xy0[model.dim[1]+1:end]
end

function get_opt_val(model::BilevelProblem)
    if model.sol !== nothing && !isempty(model.sol)
        if model.sol isa Vector{Float64}
            x_star = model.sol[1:model.dim[1]]
            y_star = model.sol[model.dim[1]+1:model.dim[1]+model.dim[2]]
            return model.F_func(x_star, y_star)
        else
            return model.sol
        end
    else
        return NaN
    end
end