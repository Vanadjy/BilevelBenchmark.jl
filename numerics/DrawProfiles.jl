using PGFPlots

function draw_convergence!(F_all_hists, N_all_hists, prob::Int, algo_names; logscale::Bool = false, max_budget::Int = 1000, type_of_ref::String = "All_Final", adjusted::Bool = false, cons_handle = "PB", start_point::String = "y0")
    F_hist = F_all_hists[prob]
    N_hist = N_all_hists[prob]
    #@assert N_hist[algo_names[1]][end] ≤ max_budget "The last value of N_hist exceeds the max budget."
    conv_plot = Union{PGFPlots.Plots.Linear, PGFPlots.Plots.Scatter}[]

    xlim_conv = 0
    for i in eachindex(algo_names)
        # Fixes xlim for plots
        x_ind_conv = findfirst(F_hist[algo_names[i]] .≤ minimum(F_hist[algo_names[i]]))
        if x_ind_conv > xlim_conv # Aims for giving to xlim the largest value so that the plot is not cut too early
            xlim_conv = min(x_ind_conv + 2 * Int(N_hist[algo_names[i]][end] / 10), length(N_hist[algo_names[i]]))
        end

        # Filters F historic if referee
        if adjusted
            filter_indexes = findall(x -> !isinf(x), F_hist[algo_names[i]])
            F_hist[algo_names[i]] = F_hist[algo_names[i]][filter_indexes]
            N_hist[algo_names[i]] = N_hist[algo_names[i]][filter_indexes]
        end
    end

    if logscale
        scatter_log = log_scale(max_budget)
        for i in eachindex(algo_names)
            # For legend display
            if !isempty(F_hist[algo_names[i]])
                legend_conv_plot = PGFPlots.Plots.Linear(N_hist[algo_names[i]][1:2], F_hist[algo_names[i]][1:2], mark = "$(marks[i])", style="$(line_color[i]), const plot", legendentry = "$(algo_names[i])")
                push!(conv_plot, legend_conv_plot)
            else
                @warn "No convergence data for algorithm $(algo_names[i]) on problem $(prob). Skipping legend entry."
            end
        end

        for i in eachindex(algo_names)
            if !isempty(F_hist[algo_names[i]])
                # Plot the convergence data
                filtered_Ns_log = filter(x -> x < N_hist[algo_names[i]][end], scatter_log)
                marker_indexes = findall(x -> x in filtered_Ns_log, N_hist[algo_names[i]])
                
                conv_plot_data = PGFPlots.Plots.Linear(N_hist[algo_names[i]], F_hist[algo_names[i]], style="$(line_color[i]), const plot, solid", mark = "none")
                conv_plot_markers = PGFPlots.Plots.Scatter(intersect(filtered_Ns_log, N_hist[algo_names[i]][marker_indexes]), F_hist[algo_names[i]][marker_indexes], style="$(line_color[i])", mark = "$(marks[i])")
                push!(conv_plot, conv_plot_data, conv_plot_markers)
            else
                @warn "No convergence data for algorithm $(algo_names[i]) on problem $(prob). Skipping plot."
            end
        end
        plt_conv = PGFPlots.Axis(
                conv_plot,
                xlabel = "Number \$N\$ of \$ F \$ evaluations",
                xmode = "log",
                ylabel = "\$ F(x_N) \$",
                title = "Convergence plot for problem \$ $(prob) \$",
                legendPos= "north east",
                #xmax = N_hist[algo_names[i]][xlim_conv] + (N_hist[algo_names[i]][xlim_conv] % 10 == 0 ? Int(10^(ceil(log10(N_hist[algo_names[i]][xlim_pxlim_converf])))) : Int(10^(ceil(log10(N_hist[algo_names[i]][xlim_conv])) - 1))) # Leave an additionnal blank space
            )
        cd("/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/Plots/ConvPlots")
        PGFPlots.save("ConvPlot-p=$prob-cons_handle=$(cons_handle)-start=$start_point-referee=$adjusted-logscale-$type_of_ref.tex", plt_conv)
        cd("/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark")

    else
        scatter_log = dec_scale(max_budget)
        for i in eachindex(algo_names)
            # For legend display
            if !isempty(F_hist[algo_names[i]])
                legend_conv_plot = PGFPlots.Plots.Linear(N_hist[algo_names[i]][1:2], F_hist[algo_names[i]][1:2], mark = "$(marks[i])", style="$(line_color[i]), const plot", legendentry = "$(algo_names[i])")
                push!(conv_plot, legend_conv_plot)
            else
                @warn "No convergence data for algorithm $(algo_names[i]) on problem $(prob). Skipping legend entry."
            end
        end
        for i in eachindex(algo_names)
            if !isempty(F_hist[algo_names[i]])
                # Plot the convergence data
                filtered_Ns_log = filter(x -> x < N_hist[algo_names[i]][end], scatter_log)
                marker_indexes = findall(x -> x in filtered_Ns_log, N_hist[algo_names[i]])
                
                conv_plot_data = PGFPlots.Plots.Linear(N_hist[algo_names[i]], F_hist[algo_names[i]], style="$(line_color[i]), const plot, solid", mark = "none")
                conv_plot_markers = PGFPlots.Plots.Scatter(intersect(filtered_Ns_log, N_hist[algo_names[i]][marker_indexes]), F_hist[algo_names[i]][marker_indexes], style="$(line_color[i])", mark = "$(marks[i])")
                push!(conv_plot, conv_plot_data, conv_plot_markers)
            else
                @warn "No convergence data for algorithm $(algo_names[i]) on problem $(prob). Skipping plot."
            end
        end
        plt_conv = PGFPlots.Axis(
            conv_plot,
            xlabel = "Number \$N\$ of \$ F \$ evaluations",
            ylabel = "\$ F(x_N) \$",
            title = "Convergence plot for problem \$ $(prob) \$",
            legendPos= "north east"
        )
        cd("/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/Plots/ConvPlots")
        PGFPlots.save("ConvPlot-p=$prob-cons_handle=$(cons_handle)-start=$start_point-referee=$adjusted-$type_of_ref.tex", plt_conv)
        cd("/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark")

    end
end

function draw_profiles!(τs, αs, ks, algo_names, prob_numbers, F_all_hists, N_all_hists; adjusted::Bool = false, cons_handle = "PB", log_scaling::Bool = true, type_of_ref::String = "All_Final", start_point::String = "y0", λ_toggle::Bool = false, λ_choice::String = "LL")
    if type_of_ref == "All_All"
        type_ref = "Total Referee"
    elseif type_of_ref == "All_Final"
        type_ref = "Final Referee"
    elseif type_of_ref == "All_Backward"
        type_ref = "Backward Referee"
    else
        type_ref = "No Referee"
    end
    for τ in τs
        for a in eachindex(algo_names)
            @views perf_profile!(y_perf[:, a], αs, F_all_hists, N_all_hists, prob_numbers, algo_names[a], τ, algo_names)
            @views data_profile!(y_data[:, a], ks, F_all_hists, N_all_hists, prob_numbers, algo_names[a], τ, algo_names; λ_toggle = λ_toggle, λ_choice = λ_choice)
            @views accuracy_profile!(y_acc[:, a], ds, F_all_hists, prob_numbers, algo_names[a], algo_names)
            #@assert y_perf[:, a][end] ≈ y_data[:, a][end] "Profiles Error: Last value of performace and date profiles should be the same. Check the code to draw profiles or increase their horizon."
        end
        perf_prof_plots = Union{PGFPlots.Plots.Linear, PGFPlots.Plots.Scatter}[]
        data_prof_plots = Union{PGFPlots.Plots.Linear, PGFPlots.Plots.Scatter}[]
        acc_prof_plots = Union{PGFPlots.Plots.Linear, PGFPlots.Plots.Scatter}[]

        xlim_perf = 0
        xlim_data = 0
        xlim_acc = 0

        for i in eachindex(algo_names)
            x_ind_perf = findfirst(y_perf[:, i] .≥ maximum(y_perf[:, i]))
            if x_ind_perf > xlim_perf # Aims for giving to xlim the largest value so that the plot is not cut too early
                xlim_perf = min(x_ind_perf + 2 * Int(αs[end] / 10), length(αs))
            end

            x_ind_data = findfirst(y_data[:, i] .≥ maximum(y_data[:, i]))
            if x_ind_data > xlim_data # Aims for giving to xlim the largest value so that the plot is not cut too early
                xlim_data = min(x_ind_data + Int(ks[end] / 10), length(ks))
            end

            x_ind_acc = findfirst(y_acc[:, i] .≥ maximum(y_acc[:, i]))
            if x_ind_acc > xlim_acc # Aims for giving to xlim the largest value so that the plot is not cut too early
                xlim_acc = min(x_ind_acc + 2 * Int(ds[end] / 10), length(ds))
            end
        end

        if log_scaling
            scatter_log_perf = log_scale(Int(αs[end]))
            scatter_log_data = log_scale(ks[end])
            scatter_log_acc = log_scale(ds[end])

            # For legend display - Only for right sided figures
            if type_ref in ["Final Referee", "Total Referee"]
                for i in eachindex(algo_names)
                    legend_perf_prof = PGFPlots.Plots.Linear(αs[1:2], [y_perf[:, i][j] for j in 1:2], mark = "$(marks[i])", style="$(line_color[i]), const plot", legendentry = "$(algo_names[i])")
                    push!(perf_prof_plots, legend_perf_prof)

                    legend_data_prof = PGFPlots.Plots.Linear(ks[1:2], [y_data[:, i][j] for j in 1:2], mark = "$(marks[i])", style="$(line_color[i]), const plot", legendentry = "$(algo_names[i])")
                    push!(data_prof_plots, legend_data_prof)

                    legend_acc_prof = PGFPlots.Plots.Linear(ds[1:2], [y_acc[:, i][j] for j in 1:2], mark = "$(marks[i])", style="$(line_color[i]), const plot", legendentry = "$(algo_names[i])")
                    push!(acc_prof_plots, legend_acc_prof)
                end
            end


            for i in eachindex(algo_names)
                perf_profile_data = PGFPlots.Plots.Linear(αs[1:xlim_perf], y_perf[:, i][1:xlim_perf], style="$(line_color[i]), const plot, solid", mark = "none")
                filtered_αs_log = filter(x -> x < αs[xlim_perf], scatter_log_perf)
                marker_indexes = findall(x -> x in filtered_αs_log, αs)
                @assert length(marker_indexes) == length(filtered_αs_log) "Marker indexes should match the filtered αs log length"
                perf_profile_markers = PGFPlots.Plots.Scatter(filtered_αs_log, y_perf[:, i][marker_indexes], style="$(line_color[i])", mark = "$(marks[i])")

                push!(perf_prof_plots, perf_profile_data, perf_profile_markers)
            end

            for i in eachindex(algo_names)
                data_profile_data = PGFPlots.Plots.Linear(ks[1:xlim_data], y_data[:, i][1:xlim_data], style="$(line_color[i]), const plot, solid", mark = "none")
                filtered_ks_log = filter(x -> x < ks[xlim_data], scatter_log_data)
                marker_indexes = findall(x -> x in filtered_ks_log, ks)
                data_profile_markers = PGFPlots.Plots.Scatter(filtered_ks_log, y_data[:, i][marker_indexes], style="$(line_color[i])", mark = "$(marks[i])")

                push!(data_prof_plots, data_profile_data, data_profile_markers)
            end

            for i in eachindex(algo_names)
                acc_profile_data = PGFPlots.Plots.Linear(ds[1:xlim_acc], y_acc[:, i][1:xlim_acc], style="$(line_color[i]), const plot, solid", mark = "none")
                filtered_ds_log = filter(x -> x < ds[xlim_acc], scatter_log_acc)
                marker_indexes = findall(x -> x in filtered_ds_log, ds)
                acc_profile_markers = PGFPlots.Plots.Scatter(filtered_ds_log, y_acc[:, i][marker_indexes], style="$(line_color[i])", mark = "$(marks[i])")

                push!(acc_prof_plots, acc_profile_data, acc_profile_markers)
            end

            plt_perf = PGFPlots.Axis(
                perf_prof_plots,
                xlabel =  "Ratio of function evaluation \$ \\alpha \$",
                xmode = "log",
                ylabel = (type_ref in ["No Referee", "Backward Referee"]) ? "Portion of \$ \\tau \$-solved instances \$ \\rho_a(\\alpha) \$" : "",
                title = "Performance profile - $type_ref",
                legendPos= "south east",
                xmin = 1.0,
                xmax = αs[xlim_perf] + (αs[xlim_perf] % 10 == 0 ? Int(10^(ceil(log10(αs[xlim_perf])))) : Int(10^(ceil(log10(αs[xlim_perf])) - 1))), # Leave an additionnal blank space
                ymin = 0.0,
                ymax = 1.0,
            )

            plt_data = PGFPlots.Axis(
                data_prof_plots,
                xlabel = λ_toggle ? (λ_choice == "LL" ? "Groups of \$ n_y + 1\$ evaluations \$k\$" : "Groups of \$ n_x + 1\$ evaluations \$k\$") : "Groups of \$ n_y(n_x + 1) \$ evaluations \$k\$",
                xmode="log",
                ylabel = (type_ref in ["No Referee", "Backward Referee"]) ? "Portion of \$ \\tau \$-solved instances \$ d_a(k) \$" : "",
                title = "Data profile - $type_ref",
                legendPos= "south east",
                xmin = 0.0,
                xmax = ks[xlim_data] + (ks[xlim_data] % 10 == 0 ? Int(10^(ceil(log10(ks[xlim_data])))) : Int(10^(ceil(log10(ks[xlim_data])) - 1))), # Leave an additionnal blank space
                ymin = 0.0,
                ymax = 1.0,
            )

            plt_acc = PGFPlots.Axis(
                acc_prof_plots,
                xlabel = "Relative accuracy \$ d \$",
                xmode="log",
                ylabel = (type_ref in ["No Referee", "Backward Referee"]) ? "Portion of \$ \\tau \$-solved instances \$ r_a(d) \$" : "",
                title = "Accuracy profile - $type_ref",
                legendPos= "north east",
                ymin = 0.0,
                ymax = 1.0,
                xmax = ds[xlim_acc] + (ds[xlim_acc] % 10 == 0 ? Int(10^(ceil(log10(ds[xlim_acc])))) : Int(10^(ceil(log10(ds[xlim_acc])) - 1))) # Leave an additionnal blank space
            )
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/Plots/PerfProfiles")
            PGFPlots.save("PerformanceProfile-tau=$(τ*100)-n_algos=$(length(algo_names))-cons_handle=$(cons_handle)-start=$start_point-referee=$adjusted-logscale-$type_of_ref-costly=$λ_choice.tex", plt_perf)
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/Plots/DataProfiles")
            PGFPlots.save("DataProfile-tau=$(τ*100)-n_algos=$(length(algo_names))-cons_handle=$(cons_handle)-start=$start_point-referee=$adjusted-logscale-$type_of_ref-costly=$λ_choice.tex", plt_data)
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/Plots/AccProfiles")
            PGFPlots.save("AccuracyProfile-n_algos=$(length(algo_names))-cons_handle=$(cons_handle)-start=$start_point-referee=$adjusted-logscale-$type_of_ref-costly=$λ_choice.tex", plt_acc)
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark")
        else
            scatter_dec_perf = dec_scale(Int(αs[end]))
            scatter_dec_data = dec_scale(ks[end])
            scatter_dec_acc = dec_scale(ds[end])

            # For legend display
            if type_ref in ["Final Referee", "Total Referee"]
                for i in eachindex(algo_names)
                    legend_perf_prof = PGFPlots.Plots.Linear(αs[1:2], [y_perf[:, i][j] for j in 1:2], mark = "$(marks[i])", style="$(line_color[i]), const plot", legendentry = "$(algo_names[i])")
                    push!(perf_prof_plots, legend_perf_prof)

                    legend_data_prof = PGFPlots.Plots.Linear(ks[1:2], [y_data[:, i][j] for j in 1:2], mark = "$(marks[i])", style="$(line_color[i]), const plot", legendentry = "$(algo_names[i])")
                    push!(data_prof_plots, legend_data_prof)

                    legend_acc_prof = PGFPlots.Plots.Linear(ds[1:2], [y_acc[:, i][j] for j in 1:2], mark = "$(marks[i])", style="$(line_color[i]), const plot", legendentry = "$(algo_names[i])")
                    push!(acc_prof_plots, legend_acc_prof)
                end
            end


            for i in eachindex(algo_names)
                perf_profile_data = PGFPlots.Plots.Linear(αs, y_perf[:, i], style="$(line_color[i]), const plot, solid", mark = "none")
                marker_indexes = findall(x -> x in scatter_dec_perf, αs)
                perf_profile_markers = PGFPlots.Plots.Scatter(αs[marker_indexes], y_perf[:, i][marker_indexes], style="$(line_color[i])", mark = "$(marks[i])")

                push!(perf_prof_plots, perf_profile_data, perf_profile_markers)
            end

            for i in eachindex(algo_names)
                data_profile_data = PGFPlots.Plots.Linear(ks, y_data[:, i], style="$(line_color[i]), const plot, solid", mark = "none")
                marker_indexes = findall(x -> x in scatter_dec_data, ks)
                data_profile_markers = PGFPlots.Plots.Scatter(ks[marker_indexes], y_data[:, i][marker_indexes], style="$(line_color[i])", mark = "$(marks[i])")

                push!(data_prof_plots, data_profile_data, data_profile_markers)
            end

            for i in eachindex(algo_names)
                acc_profile_data = PGFPlots.Plots.Linear(ds, y_acc[:, i], style="$(line_color[i]), const plot, solid", mark = "none")
                marker_indexes = findall(x -> x in scatter_dec_acc, ds)
                acc_profile_markers = PGFPlots.Plots.Scatter(ds[marker_indexes], y_acc[:, i][marker_indexes], style="$(line_color[i])", mark = "$(marks[i])")

                push!(acc_prof_plots, acc_profile_data, acc_profile_markers)
            end

            plt_perf = PGFPlots.Axis(
                perf_prof_plots,
                xlabel = "Ratio of function evaluation \$ \\alpha \$",
                ylabel = (type_ref in ["No Referee", "Backward Referee"]) ? "Portion of \$ \\tau \$-solved instances \$ \\rho_a(\\alpha) \$" : "",
                title = "Performance profile - $type_ref",
                legendPos= "south east",
                xmin = 1.0,
                ymin = 0.0,
                ymax = 1.0
            )

            plt_data = PGFPlots.Axis(
                data_prof_plots,
                xlabel = λ_toggle ? (λ_choice == "LL" ? "Groups of \$ n_y + 1\$ evaluations \$k\$" : "Groups of \$ n_x + 1\$ evaluations \$k\$") : "Groups of \$ n_y(n_x + 1) \$ evaluations \$k\$",
                ylabel = (type_ref in ["No Referee", "Backward Referee"]) ? "Portion of \$ \\tau \$-solved instances \$ d_a(k) \$" : "",
                title = "Data profile - $type_ref",
                legendPos= "south east",
                xmin = 0.0,
                ymin = 0.0,
                ymax = 1.0
            )

            plt_acc = PGFPlots.Axis(
                acc_prof_plots,
                xlabel = "Relative accuracy \$ d \$",
                ylabel = (type_ref in ["No Referee", "Backward Referee"]) ? "Portion of \$ \\tau \$-solved instances \$ r_a(d) \$" : "",
                title = "Accuracy profile - $type_ref",
                legendPos= "north east",
                xmin = 0.0,
                ymin = 0.0,
                ymax = 1.0
            )
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/Plots/PerfProfiles")
            PGFPlots.save("PerformanceProfile-tau=$(τ*100)-n_algos=$(length(algo_names))-cons_handle=$(cons_handle)-start=$start_point-referee=$adjusted-$type_of_ref-costly=$λ_choice.tex", plt_perf)
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/Plots/DataProfiles")
            PGFPlots.save("DataProfile-tau=$(τ*100)-n_algos=$(length(algo_names))-cons_handle=$(cons_handle)-start=$start_point-referee=$adjusted-$type_of_ref-costly=$λ_choice.tex", plt_data)
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/Plots/AccProfiles")
            PGFPlots.save("AccuracyProfile-n_algos=$(length(algo_names))-cons_handle=$(cons_handle)-start=$start_point-referee=$adjusted-$type_of_ref-costly=$λ_choice.tex", plt_acc)
            cd(raw"/home/dijovale/Documents/Dijon_PhD/P1-BiObjBenchmarking/BilevelBenchmark")
        end
    end
end