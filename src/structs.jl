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

mutable struct BilevelOptions
    subsolver_name::String
    γ::Float64
    oppportunistic::Bool
    ordered::Bool
    search::Bool
    orthogonal::Bool
    max_neval_upper::Int
    max_neval_upper_cons::Int
    max_neval_lower::Int
    Δ0::Float64
    tol_upper::Float64
    tol_lower::Float64
    max_time::Float64
    biphase::Bool
    verbose::Bool

    function BilevelOptions(;
        subsolver_name::String = "NOMAD",
        γ::Float64 = 1/2,
        oppportunistic::Bool = true,
        ordered::Bool = false,
        search::Bool = false,
        orthogonal::Bool = true,
        max_neval_upper::Int = 1000,
        max_neval_upper_cons::Int = 1000,
        max_neval_lower::Int = 100,
        Δ0::Float64 = 1.0,
        tol_upper::Float64 = 1e-6,
        tol_lower::Float64 = 1e-6,
        max_time::Float64 = 3600.0,
        biphase::Bool = true,
        verbose::Bool = true
    )
        @assert Δ0 > 0.0 "Value Error: Initial mesh size Δ0 must be positive"
        @assert ordered ≤ oppportunistic "Logical Error: Ordering scheme cannot be selected without an oppportunistic one"
        @assert γ > 0.0 "Value Error: γ must be positive"
        @assert γ < 1.0 "Value Error: γ must be lower than 1.0"
        return new(subsolver_name, γ, oppportunistic, ordered, search, orthogonal, max_neval_upper, max_neval_upper_cons, max_neval_lower, Δ0, tol_upper, tol_lower, max_time, biphase, verbose)
    end
end

mutable struct NOMADOptions{R, I}
    max_bb_eval::I
    quad_model_search::Bool
    direction_type::String
    eval_queue_sort::String
    max_time::R
    display_stats::Vector{String}
    display_degree::I
    cons_handle::String
    start_points::String

    function NOMADOptions{R, I}(;
        max_bb_eval::I                = 1000,
        quad_model_search::Bool       = false,
        direction_type::String        = "ORTHO 2N",
        eval_queue_sort::String       = "DIR_LAST_SUCCESS",
        max_time::R                   = 3600.0,
        display_stats::Vector{String} = ["EVAL", "SOL", "OBJ"],
        display_degree::I             = 0,
        cons_handle::String           = "PB",
        start_points::String          = "y0"
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
        @assert cons_handle ∈ ["PB", "EB"] "Constraint handling error: choose a supported way to handle constraints i.e. either EB (extreme barrier) or PB (progressive barrier)"
        @assert start_points ∈ ["y0", "yk-1"] "Start points must be either 'y0' or 'yk-1'. Other start points are not supported yet."

        return new{R, I}(
                max_bb_eval,
                quad_model_search,
                direction_type,
                eval_queue_sort,
                max_time,
                display_stats,
                display_degree,
                cons_handle,
                start_points
        )
    end
end

NOMADOptions(args...; kwargs...) = NOMADOptions{Float64, Int}(args...; kwargs...)