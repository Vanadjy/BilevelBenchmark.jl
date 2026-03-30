using JLD2, BilevelBenchmark

function generate_lambda_list(pb_bank; nb_runs::Int = 100)
    t_UL_list = zeros(length(pb_bank))
    t_LL_list = similar(t_UL_list)
    λ_list = similar(t_UL_list)
    @inbounds for i in eachindex(pb_bank)
        # Get the model and initial points
        prob = pb_bank[i]
        model = get_bilevel_problem(prob)
        nx = model.dim[1]
        ny = model.dim[2]
        x0y0 = model.xy0
        
        xk = x0y0[1:nx]
        yk = x0y0[nx+1:nx+ny]
        F = model.F_func
        f = model.f_func
        G = model.G_func
        g = model.g_func

        t_UL = 0
        t_LL = 0

        # computes the mean of the times on 100 runs
        for i in 1:nb_runs
            # Begin procedure to compute time of UL
            t_UL_start = time()
            F(xk, yk)
            G(xk, yk)
            t_UL_i = time() - t_UL_start
            t_UL += t_UL_i

            # Begin procedure to compute time of LL
            t_LL_start = time()
            f(xk, yk)
            g(xk, yk)
            t_LL_i = time() - t_LL_start
            t_LL += t_LL_i
        end

        # Compute the ratio for λ
        t_UL_list[prob] = t_UL/nb_runs
        t_LL_list[prob] = t_LL/nb_runs
        λ_list[prob] = t_LL / t_UL
    end

    save_object(joinpath("numerics/logs", "lambda_list.jld2"), λ_list)
    save_object(joinpath("numerics/logs", "t_UL_list.jld2"), t_UL_list)
    save_object(joinpath("numerics/logs", "t_LL_list.jld2"), t_LL_list)
    return λ_list
end

all_probs = collect(1:173)
issued_probs = [36, 49, 50, 51, 138, 127, 131, 173]
prob_numbers = filter(x -> !(x in issued_probs), all_probs)
generate_lambda_list(prob_numbers)