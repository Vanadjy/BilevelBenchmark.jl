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
                perf_profile_data = PGFPlots.Plots.Linear(αs, y_perf[:, i], style="$(line_color[i]), const plot, solid", mark = "none")
                marker_indexes = findall(x -> x in scatter_log, αs)
                perf_profile_markers = PGFPlots.Plots.Scatter(scatter_log, y_perf[:, i][marker_indexes], style="$(line_color[i])", mark = "$(marks[i])")

                push!(perf_prof_plots, perf_profile_data, perf_profile_markers)
            end

            for i in eachindex(algo_names)
                data_profile_data = PGFPlots.Plots.Linear(ks, y_data[:, i], style="$(line_color[i]), const plot, solid", mark = "none")
                marker_indexes = findall(x -> x in scatter_log, ks)
                data_profile_markers = PGFPlots.Plots.Scatter(scatter_log, y_data[:, i][marker_indexes], style="$(line_color[i])", mark = "$(marks[i])")

                push!(data_prof_plots, data_profile_data, data_profile_markers)
            end

            plt_perf = PGFPlots.Axis(
                perf_prof_plots,
                xlabel = "Ratio of function evaluation \$ \\alpha \$",
                xmode = "log",
                ylabel = "Porportion of problem solved",
                title = "Performance profile \$\\tau = $(Int(τ*100))\\%\$",
                legendPos= "south east",
                #ymin = 0.0,
                #ymax = 1.0
            )

            plt_data = PGFPlots.Axis(
                data_prof_plots,
                xlabel = "Groups of \$ n_p + 1\$ evaluations \$k\$",
                xmode="log",
                ylabel = "Porportion of problem solved",
                title = "Data profile \$\\tau = $(Int(τ*100))\\%\$",
                legendPos= "south east",
                #ymin = 0.0,
                #ymax = 1.0
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
                #ymin = 0.0,
                #ymax = 1.0
            )

            plt_data = PGFPlots.Axis(
                data_prof_plots,
                xlabel = "Groups of \$ n_p + 1\$ evaluations \$k\$",
                ylabel = "Porportion of problem solved",
                title = "Data profile \$\\tau = $(Int(τ*100))\\%\$",
                legendPos= "south east",
                #ymin = 0.0,
                #ymax = 1.0
            )
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/Plots/PerfProfiles")
            PGFPlots.save("PerformanceProfile-tau=$(Int(τ*100))-n_probs=$(length(all_probs))-cons_handle=$(cons_handle)-referee=$adjusted-$type_of_ref.tikz", plt_perf)
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/Plots/DataProfiles")
            PGFPlots.save("DataProfile-tau=$(Int(τ*100))-n_probs=$(length(all_probs))-cons_handle=$(cons_handle)-referee=$adjusted-$type_of_ref.tikz", plt_data)
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark")
        end
    end
end