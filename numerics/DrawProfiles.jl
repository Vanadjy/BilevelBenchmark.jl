using PGFPlots

function draw_profiles!(τs, αs, ks, algo_names, prob_numbers, F_all_hists, N_all_hists; adjusted::Bool = false, cons_handle = "PB", log_scaling::Bool = true, type_of_ref::String = "All_Final")

    for τ in τs
        for a in eachindex(algo_names)
            @views perf_profile!(y_perf[:, a], αs, F_all_hists, N_all_hists, prob_numbers, algo_names[a], τ)
            @views data_profile!(y_data[:, a], ks, F_all_hists, N_all_hists, prob_numbers, algo_names[a], τ)
            @assert y_perf[:, a][end] ≈ y_data[:, a][end] "Profiles Error: Last value of performace and date profiles should be the same. Check the code to draw profiles or increase their horizon."
        end
        perf_prof_plots = Union{PGFPlots.Plots.Linear, PGFPlots.Plots.Scatter}[]
        data_prof_plots = Union{PGFPlots.Plots.Linear, PGFPlots.Plots.Scatter}[]

        xlim_perf = 0
        xlim_data = 0

        for i in eachindex(algo_names)
            x_ind_perf = findfirst(y_perf[:, i] .≥ maximum(y_perf[:, i]))
            if x_ind_perf > xlim_perf # Aims for giving to xlim the largest value so that the plot is not cut too early
                xlim_perf = min(x_ind_perf + 2 * Int(αs[end] / 10), length(αs))
            end

            x_ind_data = findfirst(y_data[:, i] .≥ maximum(y_data[:, i]))
            if x_ind_data > xlim_data # Aims for giving to xlim the largest value so that the plot is not cut too early
                xlim_data = min(x_ind_data + Int(ks[end] / 10), length(ks))
            end
        end

        if log_scaling
            scatter_log = log_scale(αs[end])

            # For legend display
            for i in eachindex(algo_names)
                legend_perf_prof = PGFPlots.Plots.Linear(αs[1:2], [y_perf[:, i][j] for j in 1:2], mark = "$(marks[i])", style="$(line_color[i]), const plot", legendentry = "$(algo_names[i])")
                push!(perf_prof_plots, legend_perf_prof)

                legend_data_prof = PGFPlots.Plots.Linear(ks[1:2], [y_data[:, i][j] for j in 1:2], mark = "$(marks[i])", style="$(line_color[i]), const plot", legendentry = "$(algo_names[i])")
                push!(data_prof_plots, legend_data_prof)
            end


            for i in eachindex(algo_names)
                perf_profile_data = PGFPlots.Plots.Linear(αs[1:xlim_perf], y_perf[:, i][1:xlim_perf], style="$(line_color[i]), const plot, solid", mark = "none")
                filtered_αs_log = filter(x -> x < αs[xlim_perf], scatter_log)
                marker_indexes = findall(x -> x in filtered_αs_log, αs)
                @assert filtered_αs_log[end] <= xlim_perf "Filtered αs log should not exceed xlim_perf"
                @assert length(marker_indexes) == length(filtered_αs_log) "Marker indexes should match the filtered αs log length"
                perf_profile_markers = PGFPlots.Plots.Scatter(filtered_αs_log, y_perf[:, i][marker_indexes], style="$(line_color[i])", mark = "$(marks[i])")

                push!(perf_prof_plots, perf_profile_data, perf_profile_markers)
            end

            for i in eachindex(algo_names)
                data_profile_data = PGFPlots.Plots.Linear(ks[1:xlim_data], y_data[:, i][1:xlim_data], style="$(line_color[i]), const plot, solid", mark = "none")
                filtered_ks_log = filter(x -> x < ks[xlim_data], scatter_log)
                marker_indexes = findall(x -> x in filtered_ks_log, ks)
                data_profile_markers = PGFPlots.Plots.Scatter(filtered_ks_log, y_data[:, i][marker_indexes], style="$(line_color[i])", mark = "$(marks[i])")

                push!(data_prof_plots, data_profile_data, data_profile_markers)
            end

            plt_perf = PGFPlots.Axis(
                perf_prof_plots,
                xlabel = "Ratio of function evaluation \$ \\alpha \$",
                xmode = "log",
                ylabel = "Porportion of problem solved",
                title = "Performance profile \$\\tau = $(Int(τ*100))\\%\$",
                legendPos= "south east",
                ymin = 0.0,
                ymax = 1.0,
                xmax = αs[xlim_perf] + (αs[xlim_perf] % 10 == 0 ? Int(10^(ceil(log10(αs[xlim_perf])))) : Int(10^(ceil(log10(αs[xlim_perf])) - 1))) # Leave an additionnal blank space
            )

            plt_data = PGFPlots.Axis(
                data_prof_plots,
                xlabel = "Groups of \$ n_p + 1\$ evaluations \$k\$",
                xmode="log",
                ylabel = "Porportion of problem solved",
                title = "Data profile \$\\tau = $(Int(τ*100))\\%\$",
                legendPos= "south east",
                ymin = 0.0,
                ymax = 1.0,
                xmax = ks[xlim_data] + (ks[xlim_data] % 10 == 0 ? Int(10^(ceil(log10(ks[xlim_data])))) : Int(10^(ceil(log10(ks[xlim_data])) - 1))) # Leave an additionnal blank space
            )
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/Plots/PerfProfiles")
            PGFPlots.save("PerformanceProfile-tau=$(Int(τ*100))-n_probs=$(length(all_probs))-cons_handle=$(cons_handle)-referee=$adjusted-logscale-$type_of_ref.tikz", plt_perf)
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/Plots/DataProfiles")
            PGFPlots.save("DataProfile-tau=$(Int(τ*100))-n_probs=$(length(all_probs))-cons_handle=$(cons_handle)-referee=$adjusted-logscale-$type_of_ref.tikz", plt_data)
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark")
        else
            for i in eachindex(algo_names)
                perf_profile_data = PGFPlots.Plots.Linear(αs, y_perf[:, i], style="$(line_color[i]), const plot", mark = "$(marks[i])", legendentry = "$(algo_names[i])")
                push!(perf_prof_plots, perf_profile_data)
            end

            for i in eachindex(algo_names)
                data_profile_data = PGFPlots.Plots.Linear(ks, y_data[:, i], style="$(line_color[i]), const plot", mark = "$(marks[i])", legendentry = "$(algo_names[i])")
                push!(data_prof_plots, data_profile_data)
            end

            plt_perf = PGFPlots.Axis(
                perf_prof_plots,
                xlabel = "Ratio of function evaluation \$ \\alpha \$",
                ylabel = "Porportion of problem solved",
                title = "Performance profile \$\\tau = $(Int(τ*100))\\%\$",
                legendPos= "south east",
                ymin = 0.0,
                ymax = 1.0
            )

            plt_data = PGFPlots.Axis(
                data_prof_plots,
                xlabel = "Groups of \$ n_p + 1\$ evaluations \$k\$",
                ylabel = "Porportion of problem solved",
                title = "Data profile \$\\tau = $(Int(τ*100))\\%\$",
                legendPos= "south east",
                ymin = 0.0,
                ymax = 1.0
            )
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/Plots/PerfProfiles")
            PGFPlots.save("PerformanceProfile-tau=$(Int(τ*100))-n_probs=$(length(all_probs))-cons_handle=$(cons_handle)-referee=$adjusted-$type_of_ref.tikz", plt_perf)
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/Plots/DataProfiles")
            PGFPlots.save("DataProfile-tau=$(Int(τ*100))-n_probs=$(length(all_probs))-cons_handle=$(cons_handle)-referee=$adjusted-$type_of_ref.tikz", plt_data)
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark")
        end
    end
end