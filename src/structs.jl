export BilevelProblem, NOMADOptions

mutable struct BilevelProblem
    name::String
    dim::Vector{Int}  # [n_x, n_y, n_G, n_g]
    xy0::Vector{Float64}
    Ff::Vector{Float64}  # [F, f, Status]
    F_func::Function
    f_func::Function
    G_func::Function
    g_func::Function
    sol::Union{Vector{Float64}, Float64}  # Solution or optimal value
end

mutable struct NOMADOptions{R, I}
    max_bb_eval::I
    quad_model_search::Bool
    direction_type::String
    eval_queue_sort::String
    max_time::R
    display_stats::Vector{String}
    display_degree::I

    function NOMADOptions{R, I}(;
        max_bb_eval::I                = 1000,
        quad_model_search::Bool       = false,
        direction_type::String        = "ORTHO 2N",
        eval_queue_sort::String       = "DIR_LAST_SUCCESS",
        max_time::R                   = 3600.0,
        display_stats::Vector{String} = ["EVAL", "SOL", "OBJ"],
        display_degree::I             = 0
        ) where {R <: Real, I <: Int}
        DirectionTypes = ["ORTHO 2N" # 2n directions, no quadratic models
                        "ORTHO N+1 NEG" # n directions, the (n+1)th is the negative sum of the n first.
                        "ORTHO N+1 QUAD" # n directions, the (n+1)th is found by solving a quadratic subproblem
                        "ORTHO N+1 QUAD" # n directions, the (n+1)th is found by solving a quadratic subproblem
                        "N+1 UNI" # n+1 uniformly distributed directions
                        "SINGLE" # one direction
                        "DOUBLE" # two opposed direction
                        ]
        @assert max_bb_eval > 0 "Number of black box evaluations must be positive"
        @assert direction_type ∈ DirectionTypes "Direction type indicated not supported. Please use a direction type among $DirectionTypes"
        @assert max_time > 0.0 "Need a positive time budget"
        @assert display_degree ∈ [0, 1, 2, 3] "Display degree indicated not supported in NOMAD. It must be an 0 to 3 integer."

        return new{R, I}(
                max_bb_eval,
                quad_model_search,
                direction_type,
                eval_queue_sort,
                max_time,
                display_stats,
                display_degree
        )
    end
end

NOMADOptions(args...; kwargs...) = NOMADOptions{Float64, Int}(args...; kwargs...)