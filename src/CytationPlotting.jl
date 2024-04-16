using JSON
using Plots, LaTeXStrings, StatsPlots
using DataFrames, CSV
using Statistics
using HypothesisTests: UnequalVarianceTTest, pvalue
using LsqFit
using Colors: JULIA_LOGO_COLORS
    
pgfplotsx()

function dose_response(conditions, plot_conditions, 
                      plot_normalization, normalization_method, plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, dose_concs, plot_size, plots_directory)
    data = [combine(df, names(df) .=> maximum .=> names(df)) for df in data]
    x = dose_concs
    y = []
    error = []

    for condition in plot_conditions
        condition_data = [] 
        for i in eachindex(conditions)
            if condition in keys(conditions[i])
                subset = data[i][!, conditions[i][condition]]
                for entry in subset[1, :]
                    push!(condition_data, entry)
                end
            end
        end
        means = mean(condition_data)
        push!(y, means[1])
        stds = std(condition_data)
        push!(error, stds[1])
    end
    
    if plot_normalization != ""
        norm_data = [] 
        for i in eachindex(conditions)
            if plot_normalization in keys(condition[i])
                subset = data[i][!, conditions[i][plot_normalization]]
                for entry in subset[1, :]
                    push!(norm_data, entry)
                end
            end
        end
        norm_mean = mean(norm_data)
        norm_std = std(norm_data)
        if normalization_method == "percent"
            y = (y ./ norm_mean .- 1) .* 100
            error = (sqrt.((error ./ norm_mean).^2 .+ (norm_std / norm_mean)^2) .- 1) .* 100
        elseif normalization_method == "fold-change"
            y = y ./ norm_mean
            error = sqrt.((error ./ norm_mean).^2 .+ (norm_std / norm_mean)^2)
        else
            error("Normalization method not implemented yet!")
        end
    end

    error .= ifelse.(isnan.(error), 0, error)

    @. model(x, p) = p[1] + (x^p[4])*(p[2]-p[1])/(x^p[4]+p[3]^p[4]) 

    p0 = [0.0, maximum(y), 0.01, 2]
    lb = [0.0, 0.0, 0.0, 0.0]
    ub = [Inf, Inf, maximum(x), 20.0]

    fit = curve_fit(model, x, y, p0, lower=lb, upper=ub)
    pstar = coef(fit)
    plt = scatter(x, y, yerror=error, color="black", mc=:white, label="Data", size=plot_size)
    xbase = collect(range(minimum(x), maximum(x), 100))
    plot!(plt, xbase, model.(xbase, (pstar,)), color="black", label="Fit")
    if occursin("\$", plot_xlabel) 
        plot_xlabel = latexstring(plot_xlabel)
    end
    if occursin("\$", plot_ylabel[1])
        plot_ylabel[1] = latexstring(plot_ylabel[1])
    end
    if occursin("\$", plot_title)
        plot_title = latexstring(plot_title)
    end
    xlabel!(plt, plot_xlabel)
    ylabel!(plt, plot_ylabel[1])
    title!(plt, plot_title)
    savefig(plt, "$plots_directory/$plot_filename"*".svg")
end

function well_name_to_position(well_name::String)
    row = Int(well_name[1]) - Int('A') + 1
    col = parse(Int, well_name[2:end])
    return row, col
end

function find_block_boundaries(well_names::Array{String, 1})
    min_row, max_row = 8, 1
    min_col, max_col = 12, 1

    for well_name in well_names
        row, col = well_name_to_position(well_name)
        min_row = min(min_row, row)
        max_row = max(max_row, row)
        min_col = min(min_col, col)
        max_col = max(max_col, col)
    end

    return min_row, max_row, min_col, max_col
end

function heatplot(conditions, plot_conditions, 
                      plot_normalization, normalization_method, plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, plot_clab, plot_size, plots_directory)

    data = [combine(df, names(df) .=> maximum .=> names(df)) for df in data]
    plot_condition = [key for key in plot_conditions][1] 
    replicates = 0
    for i in eachindex(conditions)
        if plot_condition in keys(conditions[i])
            subset = conditions[i][plot_condition]
            for well in subset 
                replicates += 1
            end
        end
    end

    block_wells = Array{String, 2}(undef, length(plot_conditions), replicates)
    plate = Array{Int, 1}(undef, replicates)
    plate_matrix = zeros(8, 12, length(conditions))
    for (j, condition) in enumerate(plot_conditions)
        replicate = 0
        for i in eachindex(conditions)
            if condition in keys(conditions[i])
                for well in conditions[i][condition]
                    replicate += 1
                    row, col = well_name_to_position(well)
                    plate_matrix[row, col, i] = data[i][1, Symbol(well)]
                    block_wells[j, replicate] = well
                    if j == 1
                        plate[replicate] = i
                    end
                end
            end
        end
    end

    min_row, max_row, min_col, max_col = find_block_boundaries(block_wells[:,1])
    nrows = max_row - min_row + 1
    ncols = max_col - min_col + 1
    block_matrices = Array{Float64, 3}(undef, nrows, ncols, replicates)
    for i in 1:replicates 
        min_row, max_row, min_col, max_col = find_block_boundaries(block_wells[:,i])
        block_matrices[:,:,i] = plate_matrix[min_row:max_row, min_col:max_col, plate[i]]
    end

    block_mean = mean(block_matrices, dims=3)
    block_std = std(block_matrices, dims=3)

    if plot_normalization != ""
        norms = []
        for i in eachindex(conditions)
            if plot_normalization in keys(conditions[i])
                subset = data[i][!, conditions[i][plot_normalization]]
                for entry in subset[1, :]
                    push!(norms, entry)
                end
            end
        end
        norm_mean = mean(norms)
        norm_std = std(norms)
        if normalization_method == "percent"
            block_std = (sqrt.((block_std ./ norm_mean).^2 .+ (norm_std / norm_mean)^2) .- 1) .* 100
            block_mean = (block_mean ./ norm_mean .- 1) .* 100
            clab_title = "%"
        elseif normalization_method == "fold-change"
            block_mean = block_mean ./ norm_mean
            block_std = sqrt.((block_std ./ norm_mean).^2 .+ (norm_std / norm_mean)^2)
            clab_title = "Fold-change"
        else
            error("Normalization method not implemented yet!")
        end
    end
    if occursin("\$", plot_xlabel)
        plot_xlabel = latexstring(plot_xlabel)
    end
    if occursin("\$", plot_ylabel[1])
        plot_ylabel[1] = latexstring(plot_ylabel[1])
    end
    if occursin("\$", clab_title)
        clab_title = latexstring(clab_title)
    end
    if occursin("\$", plot_title)
        plot_title = latexstring(plot_title)
    end
    plt = heatmap(plot_xticks, reverse(plot_yticks), block_mean[:,:,1], xrotation = 45, yflip=true, color=:cool,
                  colorbar_title=clab_title, size=plot_size)
    xlabel!(plt, plot_xlabel)
    ylabel!(plt, plot_ylabel[1])
    title!(plt, plot_title)
    savefig("$plots_directory/$plot_filename"*"_mean.svg")
    plt = heatmap(plot_xticks, reverse(plot_yticks), block_std[:,:,1], xrotation = 45, yflip=true, color=:cool,
                  colorbar_title=clab_title, size=plot_size)
    xlabel!(plt, plot_xlabel)
    ylabel!(plt, plot_ylabel[1])
    title!(plt, plot_title)
    savefig(plt, "$plots_directory/$plot_filename"*"_std.svg")
end

function twin_y(conditions, plot_conditions,
                      plot_normalization, normalization_method, plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, xaxis_data, plot_xaxis, plot_size, plots_directory)
    colors = [JULIA_LOGO_COLORS.green, JULIA_LOGO_COLORS.purple]
    p = plot(size=plot_size)
    p_twin = twinx(p)
    
    for j in 1:length(data) 
        if plot_normalization != ""
            condition_data = DataFrame() 
            for i in eachindex(conditions)
                if plot_normalization in keys(conditions[i])
                    subset = data[j][i][!, conditions[i][plot_normalization]]
                    subset = DataFrame(collect.(eachrow(subset)), :auto)
                    append!(condition_data, subset)
                end
            end
            norm_means = maximum.(eachrow(condition_data))
            max_norm = mean(norm_means)
            max_std = std(norm_means)
        end

        for condition in plot_conditions
            condition_data = DataFrame() 
            if plot_xaxis == "OD"
                OD = DataFrame()
            end
            for i in eachindex(conditions)
                if condition in keys(conditions[i])
                    subset = data[j][i][!, conditions[i][condition]]
                    subset = DataFrame(collect.(eachrow(subset)), :auto)
                    append!(condition_data, subset)
                    if plot_xaxis == "OD"
                        subset = xaxis_data[i][!, conditions[i][condition]]
                        subset = DataFrame(collect.(eachrow(subset)), :auto)
                        append!(OD, subset)
                    end
                end
            end
            means = mean.(eachcol(condition_data))
            stds = std.(eachcol(condition_data))
            if plot_normalization != ""
                norm_method_idx = min(length(normalization_method), j)
                if normalization_method[norm_method_idx] == "percent"
                    means = (means ./ max_norm .- 1) .* 100
                    stds = (sqrt.((stds ./ max_norm).^2 .+ (max_std / max_norm)^2) .- 1) .* 100
                elseif normalization_method[norm_method_idx] == "fold-change"
                    means = means ./ max_norm
                    stds = sqrt.((stds ./ max_norm).^2 .+ (max_std / max_norm)^2)
                end
            end
            if occursin("\$", plot_ylabel[j])
                plot_ylabel[j] = latexstring(plot_ylabel[j])
            end
            if plot_xaxis == "Time"
                xaxis=xaxis_data
            elseif plot_xaxis == "OD"
                xaxis = mean.(eachcol(OD))
            else
                error("Can only plot time or OD on the x-axis.")
            end
            if j == 1
                stds .= ifelse.(isnan.(stds), 0, stds)
                plot!(p, xaxis, means, marker=:circle, ribbon=stds, label=latexstring(condition), color=colors[j], 
                      ytickfontcolor=colors[j], legend=false)
                ylabel!(p, plot_ylabel[j], yguidefontcolor=colors[j])
            else
                stds .= ifelse.(isnan.(stds), 0, stds)
                plot!(p_twin, xaxis, means, marker=:circle, ribbon=stds, 
                      label=latexstring(condition), color=colors[j], 
                      ytickfontcolor=colors[j], legend=false)
                ylabel!(p_twin, plot_ylabel[j], yguidefontcolor=colors[j])
            end
        end
    end
    if occursin("\$", plot_xlabel)
        plot_xlabel = latexstring(plot_xlabel)
    end
    if occursin("\$", plot_title)
        plot_title = latexstring(plot_title)
    end
    xlabel!(p, plot_xlabel)
    xlabel!(p_twin, "")
    title!(p, plot_title)
    savefig(p, "$plots_directory/$plot_filename"*".svg")
end

function line_plot(conditions, plot_conditions, 
                      plot_normalization, normalization_method, plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, xaxis_data, plot_xaxis, plot_size, plots_directory)

    if occursin("\$", plot_xlabel)
        plot_xlabel = latexstring(plot_xlabel)
    end
    if occursin("\$", plot_ylabel[1])
        plot_ylabel[1] = latexstring(plot_ylabel[1])
    end
    if occursin("\$", plot_title)
        plot_title = latexstring(plot_title)
    end
    p = plot(size=plot_size)

    if plot_normalization != ""
        condition_data = DataFrame() 
        for i in eachindex(conditions)
            if plot_normalization in keys(conditions[i])
                subset = data[i][!, conditions[i][plot_normalization]]
                subset = DataFrame(collect.(eachrow(subset)), :auto)
                append!(condition_data, subset)
            end
        end
        norm_means = maximum.(eachrow(condition_data))
        max_norm = mean(norm_means)
        max_std = std(norm_means)
    end

    for condition in plot_conditions
        if plot_xaxis == "OD"
            OD = DataFrame()
        end
        condition_data = DataFrame() 
        for i in eachindex(conditions)
            if condition in keys(conditions[i])
                subset = data[i][!, conditions[i][condition]]
                subset = DataFrame(collect.(eachrow(subset)), :auto)
                append!(condition_data, subset)
                if plot_xaxis == "OD"
                    subset = xaxis_data[i][!, conditions[i][condition]]
                    subset = DataFrame(collect.(eachrow(subset)), :auto)
                    append!(OD, subset)
                end
            end
        end
        means = mean.(eachcol(condition_data))
        stds = std.(eachcol(condition_data))
        if plot_normalization != ""
            if normalization_method == "percent"
                means = (means ./ max_norm .- 1) .* 100
                stds = (sqrt.((stds ./ max_norm).^2 .+ (max_std / max_norm)^2) .- 1) .* 100
            elseif normalization_method == "fold-change"
                means = means ./ max_norm
                stds = sqrt.((stds ./ max_norm).^2 .+ (max_std / max_norm)^2)
            else
                error("Normalization method not implemented yet!")
            end
        end
        if plot_xaxis == "Time"
            xaxis=xaxis_data
        elseif plot_xaxis == "OD"
            xaxis = mean.(eachcol(OD))
        else
            error("Can only plot time or OD on the x-axis.")
        end
        stds .= ifelse.(isnan.(stds), 0, stds)
        plot!(p, xaxis, means, marker=:circle, ribbon=stds, label=latexstring(condition))
    end
    xlabel!(p, plot_xlabel)
    ylabel!(p, plot_ylabel[1])
    title!(p, plot_title)
    savefig(p, "$plots_directory/$plot_filename"*".svg")
end

function jitter_vals(values; width=0.05)
    return values .+ width .* (rand(length(values)) .- 0.5)
end

function pval_annotation(pval)
    if pval < 0.0001
        return "****"
    elseif pval < 0.001
        return "***"
    elseif pval < 0.01
        return "**"
    elseif pval < 0.05
        return "*"
    else
        return "n.s."
    end
end

function sig_annot!(p, categories, unique_cats, values)
	ymin, ymax = ylims(p)
	dy = (ymax - ymin)/10
    ymax += dy * length(unique_cats)
	xt = xticks(p[1])[1]

	plot!(ylim=(ymin, ymax + dy))

	for (i, category) in enumerate(unique_cats[2:end])
		data1 = values[categories .== unique_cats[1]]
		data2 = values[categories .== category]

		test_result = UnequalVarianceTTest(data1, data2)
		pval = pvalue(test_result)

		xpos = [xt[1], xt[i+1]]
        ypos = ymax - i*dy

		plot!(xpos, [ypos, ypos], color=:black)
		annotate!(mean(xpos), ypos + dy/4, text(pval_annotation(pval), 10))
	end
end

function jitter_plot(conditions, plot_conditions, 
                      plot_normalization, normalization_method, plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, default_color, plot_size, plots_directory)
    data = [combine(df, names(df) .=> maximum .=> names(df)) for df in data]

    if plot_normalization != ""
        norm_data = [] 
        for i in eachindex(conditions)
            if plot_normalization in keys(conditions[i])
                subset = data[i][!, conditions[i][plot_normalization]]
                for entry in subset[1, :]
                    push!(norm_data, entry)
                end
            end
        end
        norm_mean = mean(norm_data)
    end

    categories = String[]
    values = Float64[]
    cat_labels = String[]

    for condition in plot_conditions
        for i in eachindex(conditions)
            if condition in keys(conditions[i])
                for col in conditions[i][condition]
                    append!(categories, repeat([condition], length(data[i][!, col])))
                    if normalization_method == "percent"
                        append!(values, (data[i][!, col] ./ norm_mean .- 1) .* 100)
                   elseif normalization_method == "fold-change"
                        append!(values, data[i][!, col] ./ norm_mean)
                    else
                        error("Normalization method not implemented yet!")
                    end
                    append!(cat_labels, repeat([col], length(data[i][!, col])))
                end
            end
        end
    end

    unique_cats = unique(categories)
    cat_indices = Dict(cat => i for (i, cat) in enumerate(unique_cats))
    x_vals = [cat_indices[cat] + jitter_vals([0]; width=0.1)[1] for cat in categories]
    x_min = minimum(x_vals) - 0.5
    x_max = maximum(x_vals) + 0.5
    box_x = [cat_indices[cat] for cat in categories]
    if occursin("\$", plot_xlabel)
        plot_xlabel = latexstring(plot_xlabel)
    end
    if occursin("\$", plot_ylabel[1])
        plot_ylabel[1] = latexstring(plot_ylabel[1])
    end
    if occursin("\$", plot_title)
        plot_title = latexstring(plot_title)
    end
    p = scatter(x_vals, values, group=categories, color=default_color, 
                markerstrokecolor=default_color, alpha=0.4, 
                xticks=(1:length(unique_cats), unique_cats), xrotation=45, 
                xlims=(x_min, x_max), size=plot_size)
    boxplot!(p, box_x, values, color=default_color, linecolor=default_color, 
             markerstrokecolor=default_color, leg=false, outliers=false, 
             fillalpha=0.1, linewidth=1.5)
	sig_annot!(p, categories, unique_cats, values)
    xlabel!(p, plot_xlabel)
    ylabel!(p, plot_ylabel[1])
    title!(p, plot_title)
    savefig(p, "$plots_directory/$plot_filename"*".svg")
end

function select_data(plot_dtype, lum, OD, lux, BF_imaging, CFP_imaging, 
                    YFP_imaging, texas_red_imaging, CY5_imaging, YFP, CY5)
    if plot_dtype == "lum"
        data = lum
    elseif plot_dtype == "OD"
        data = OD
    elseif plot_dtype == "BF_imaging"
        data = BF_imaging
    elseif plot_dtype == "CFP_imaging"
        data = CFP_imaging
    elseif plot_dtype == "YFP_imaging"
        data = YFP_imaging
    elseif plot_dtype == "texas_red_imaging"
        data = texas_red_imaging
    elseif plot_dtype == "CY5_imaging"
        data = CY5_imaging
    elseif plot_dtype == "YFP"
        data = YFP
    elseif plot_dtype == "CY5"
        data = CY5
    elseif plot_dtype == "RLU" 
        data = lux
    else
        error("Can't select a data type that doesn't exist.")
    end
    return data
end

function generate_plot(conditions, acquisition_frequency, plot_num, plot_type, 
                      plot_dtypes, plot_numerators, plot_denominators, 
                      plot_xaxis, plot_conditions, 
                      plot_normalization, normalization_method, plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_clab, plot_size, plot_filename, 
                      lum, OD, lux, BF_imaging, CFP_imaging, YFP_imaging, 
                      texas_red_imaging, CY5_imaging, YFP, CY5, 
                      default_color, dose_concs, plots_directory)
    plot_size = Tuple(plot_size)
    if length(plot_dtypes) == 1
        data = select_data(plot_dtypes[1], lum, OD, lux, BF_imaging, CFP_imaging, 
                           YFP_imaging, texas_red_imaging, CY5_imaging, YFP, CY5)
    elseif length(plot_dtypes) == 2 && plot_type != "two-axis" 
        if length(plot_numerators) == 0
            error("Passed too many data types without specifying numerator and denominator.")
        elseif length(plot_numerators) > 1
            error("Passed too many numerators/denominators.")
        else
            numerator = select_data(plot_numerators[1], lum, OD, lux, BF_imaging, 
                                  CFP_imaging, YFP_imaging, texas_red_imaging, 
                                  CY5_imaging, YFP, CY5)
            denominator = select_data(plot_denominators[1], lum, OD, lux, BF_imaging, 
                                  CFP_imaging, YFP_imaging, texas_red_imaging, 
                                  CY5_imaging, YFP, CY5)
            quotient_df = DataFrame()
            for col_name in names(numerator[1])
                quotient_df[!, col_name] = numerator[1][!, col_name] ./ denominator[1][!, col_name]
            end
            quotient_df .= ifelse.(isnan.(quotient_df), 0, quotient_df)
            data = quotient_df 
        end
    elseif length(plot_dtypes) == 2 && plot_type == "two-axis"
        data = Array{Vector{Union{Nothing, DataFrame}}, 1}(undef, 2)
        for i in 1:2
            data[i] = select_data(plot_dtypes[i], lum, OD, lux, BF_imaging, 
                                  CFP_imaging, YFP_imaging, texas_red_imaging, 
                                  CY5_imaging, YFP, CY5)
        end
    elseif length(plot_dtypes) > 2 && plot_type == "two-axis"
        if length(plot_numerators) == 0
            error("Passed too many data types without specifying numerator and denominator.")
        elseif length(plot_numerators) > 2
            error("Passed too many numerators/denominators.")
        else
            data = Array{Vector{Union{Nothing, DataFrame}}, 1}(undef, 2)
            plot_dtype1 = plot_dtypes[1]
            if in(plot_dtype1, plot_denominators) || in(plot_dtype1, plot_denominators)  
                numerator = select_data(plot_numerators[1], lum, OD, lux, BF_imaging, 
                                      CFP_imaging, YFP_imaging, texas_red_imaging, 
                                      CY5_imaging, YFP, CY5)
                denominator = select_data(plot_denominators[1], lum, OD, lux, BF_imaging, 
                                      CFP_imaging, YFP_imaging, texas_red_imaging, 
                                      CY5_imaging, YFP, CY5)
                quotient_df = DataFrame()
                for col_name in names(numerator[1])
                    quotient_df[!, col_name] = numerator[1][!, col_name] ./ denominator[1][!, col_name]
                end
                quotient_df .= ifelse.(isnan.(quotient_df), 0, quotient_df)
                data[1] = [quotient_df] 
                if length(plot_dtypes) == 3
                    plot_dtype = filter(x -> x != plot_numerators[1] && x != plot_denominators[1], plot_dtypes)[1]
                    data[2] = select_data(plot_dtype, lum, OD, lux, BF_imaging, 
                                      CFP_imaging, YFP_imaging, texas_red_imaging, 
                                      CY5_imaging, YFP, CY5) 
                else
                    numerator = select_data(plot_numerators[2], lum, OD, lux, BF_imaging, 
                                          CFP_imaging, YFP_imaging, texas_red_imaging, 
                                          CY5_imaging, YFP, CY5)
                    denominator = select_data(plot_denominators[2], lum, OD, lux, BF_imaging, 
                                          CFP_imaging, YFP_imaging, texas_red_imaging, 
                                          CY5_imaging, YFP, CY5)
                    quotient_df = DataFrame()
                    for col_name in names(numerator[1])
                        quotient_df[!, col_name] = numerator[1][!, col_name] ./ denominator[1][!, col_name]
                    end
                    quotient_df .= ifelse.(isnan.(quotient_df), 0, quotient_df)
                    data[2] = [quotient_df] 
                end
            else
                data[1] = select_data(plot_dtype1, lum, OD, lux, BF_imaging, 
                                      CFP_imaging, YFP_imaging, texas_red_imaging, 
                                      CY5_imaging, YFP, CY5)
                numerator = select_data(plot_numerators[1], lum, OD, lux, BF_imaging, 
                                      CFP_imaging, YFP_imaging, texas_red_imaging, 
                                      CY5_imaging, YFP, CY5)
                denominator = select_data(plot_denominators[1], lum, OD, lux, BF_imaging, 
                                      CFP_imaging, YFP_imaging, texas_red_imaging, 
                                      CY5_imaging, YFP, CY5)
                quotient_df = DataFrame()
                for col_name in names(numerator[1])
                    quotient_df[!, col_name] = numerator[1][!, col_name] ./ denominator[1][!, col_name]
                end
                quotient_df .= ifelse.(isnan.(quotient_df), 0, quotient_df)
                data[2] = [quotient_df] 
            end
        end
    else
        error("Passed too many data types.")
    end
    if plot_conditions[1] == "all"
        plot_conditions = keys(conditions[1])
    end
    if plot_type == "jitter"
        jitter_plot(conditions, plot_conditions, 
                      plot_normalization, normalization_method[1], plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, default_color, plot_size, 
                      plots_directory)
    elseif plot_type == "line"
        if plot_xaxis == "Time"
            t = range(0,stop=nrow(data[1]) - 1,length=nrow(data[1])) ./ acquisition_frequency 
            line_plot(conditions, plot_conditions, 
                          plot_normalization, normalization_method[1], plot_title, 
                          plot_ylabel, plot_xlabel, 
                          plot_yticks, plot_xticks, 
                          plot_filename, data, t, plot_xaxis, plot_size,
                          plots_directory)
        elseif plot_xaxis == "OD"
            line_plot(conditions, plot_conditions, 
                          plot_normalization, normalization_method[1], plot_title, 
                          plot_ylabel, plot_xlabel, 
                          plot_yticks, plot_xticks, 
                          plot_filename, data, OD, plot_xaxis, plot_size,
                          plots_directory)
        else
            error("Can only plot time or OD on the x-axis of a lineplot.")
        end
    elseif plot_type == "two-axis"
        if plot_xaxis == "Time"
            t = range(0,stop=nrow(data[1][1]) - 1,length=nrow(data[1][1])) ./ acquisition_frequency 
            twin_y(conditions, plot_conditions, 
                          plot_normalization, normalization_method, plot_title, 
                          plot_ylabel, plot_xlabel, 
                          plot_yticks, plot_xticks, 
                          plot_filename, data, t, plot_xaxis, plot_size,
                          plots_directory)
        elseif plot_xaxis == "OD"
            twin_y(conditions, plot_conditions, 
                          plot_normalization, normalization_method, plot_title, 
                          plot_ylabel, plot_xlabel, 
                          plot_yticks, plot_xticks, 
                          plot_filename, data, OD, plot_xaxis, plot_size,
                          plots_directory)
        else
            error("Can only plot time or OD on the x-axis of a two-axis plot.")
        end
    elseif plot_type == "heatmap"
        heatplot(conditions, plot_conditions, 
                      plot_normalization, normalization_method[1], plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, plot_clab, plot_size,
                      plots_directory)
    elseif plot_type == "dose-response"
        dose_response(conditions, plot_conditions, 
                      plot_normalization, normalization_method[1], plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, dose_concs, plot_size,
                      plots_directory)
    else
        error("Plot type not implemented yet!")
    end
end

function main()
    #push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\usepackage{sfmath}\n\\renewcommand{\\familydefault}{\\sfdefault}")
    push!(PGFPlotsX.CUSTOM_PREAMBLE,
          """
    \\usepackage[scaled]{helvet}
    \\renewcommand\\familydefault{\\sfdefault} 
    \\usepackage[T1]{fontenc}
    \\usepackage{helvet, sansmath}
    \\sansmath
    """)
    default(titlefont = (20, "computer modern"), legendfontsize = 15, 
            guidefont = (20, :black), colorbar_tickfontsize=15, colorbar_titlefontsize=15, tickfont = (15, :black), 
            guide = L"x", linewidth=2, grid=false, formatter = :plain)
    default_color = RGB{Float64}(0.8888735002725198,0.43564919034818994,0.2781229361419438) 
    
    config = JSON.parsefile("experiment_config.json")
    conditions = config["conditions"]
    nplates = length(conditions)
    acquisition_frequency = config["acquisition_frequency"]
    dose_concs = config["dose_concs"]
    images_directories = config["images_directory"]
    bulk_data = config["bulk_data"]
    plot_types = config["plot_types"] 
    plot_dtypes = config["plot_dtypes"] 
    plot_numerators = config["plot_numerators"] 
    plot_denominators = config["plot_denominators"] 
    plot_xaxis = config["plot_xaxis"] 
    plot_conditions = config["plot_conditions"] 
    plot_normalization = config["plot_normalization"] 
    normalization_method = config["normalization_method"] 
    plot_titles = config["plot_titles"]
    plot_ylabels = config["plot_ylabs"] 
    plot_xlabels = config["plot_xlabs"] 
    plot_xticks = config["plot_xticks"] 
    plot_yticks = config["plot_xticks"] 
    plot_clabs = config["color_label"] 
    plot_size = config["plot_size"] 
    plot_filenames = config["plot_filenames"] 
    parent_directory = length(images_directories) > 0 ? images_directories[1] : bulk_data[1][1:end-4] 
    plots_directory = "$parent_directory/Plots"
    if isdir(plots_directory)
        rm(plots_directory; recursive = true)
    end
    mkdir(plots_directory)

    data_types = vcat(values(plot_dtypes)...)
    lum = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    OD = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    lux = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    BF_imaging = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    CFP_imaging = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    YFP_imaging = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    texas_red_imaging = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    CY5_imaging = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    YFP = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    CY5 = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    for i in 1:nplates
        if images_directories[i] == "nothing"
            images_directories[i] = bulk_data[i][1:end-4]
        end
        data_directory = images_directories[i]
        if "lum" in data_types || "RLU" in data_types && isfile("$data_directory/lum.csv")
            lum[i] = CSV.read("$data_directory/lum.csv", DataFrame)
        else
            lum[i] = nothing
        end
        if "OD" in data_types || "RLU" in data_types && isfile("$data_directory/OD600.csv")
            OD[i] = CSV.read("$data_directory/OD600.csv", DataFrame)
        else
            OD[i] = nothing
        end
        if "BF_imaging" in data_types && isfile("$data_directory/results_data/BF_imaging.csv")
            BF_imaging[i] = CSV.read("$data_directory/results_data/BF_imaging.csv", DataFrame)
        else  
            BF_imaging[i] = nothing
        end
        if "CFP_imaging" in data_types && isfile("$data_directory/results_data/CFP_imaging.csv")
            CFP_imaging[i] = CSV.read("$data_directory/results_data/CFP_imaging.csv", DataFrame)
        else  
            CFP_imaging[i] = nothing
        end
        if "YFP_imaging" in data_types && isfile("$data_directory/results_data/YFP_imaging.csv")
            YFP_imaging[i] = CSV.read("$data_directory/results_data/YFP_imaging.csv", DataFrame)
        else  
            YFP_imaging[i] = nothing
        end
        if "texas_red_imaging" in data_types && isfile("$data_directory/results_data/texas_red_imaging.csv")
            texas_red_imaging[i] = CSV.read("$data_directory/results_data/texas_red_imaging.csv", DataFrame)
        else  
            texas_red_imaging[i] = nothing
        end
        if "CY5_imaging" in data_types && isfile("$data_directory/results_data/CY5_imaging.csv")
            CY5_imaging[i] = CSV.read("$data_directory/results_data/CY5_imaging.csv", DataFrame)
        else  
            CY5_imaging[i] = nothing
        end
        if "YFP" in data_types && isfile("$data_directory/YFP.csv")
            YFP[i] = CSV.read("$data_directory/YFP.csv", DataFrame)
        else  
            YFP[i] = nothing
        end
        if "CY5" in data_types && isfile("$data_directory/CY5.csv")
            CY5[i] = CSV.read("$data_directory/CY5.csv", DataFrame)
        else  
            CY5[i] = nothing
        end
        if "RLU" in data_types && lum[i] != nothing && OD[i] != nothing 
            lux[i] = lum[i] ./ OD[i]
        else
            lux[i] = nothing
        end
    end
    println("Plotting...")
    for (plot_num, plot_type) in plot_types
        generate_plot(conditions, acquisition_frequency, plot_num, plot_type, 
                      plot_dtypes[plot_num], plot_numerators[plot_num], plot_denominators[plot_num], 
                      plot_xaxis[plot_num], plot_conditions[plot_num], 
                      plot_normalization[plot_num], normalization_method[plot_num],
                      plot_titles[plot_num], 
                      plot_ylabels[plot_num], plot_xlabels[plot_num], 
                      plot_yticks[plot_num], plot_xticks[plot_num], 
                      plot_clabs[plot_num], plot_size[plot_num], plot_filenames[plot_num], 
                      lum, OD, lux, BF_imaging, CFP_imaging, YFP_imaging, 
                      texas_red_imaging, CY5_imaging, YFP, CY5, 
                      default_color, dose_concs[plot_num], plots_directory)
    end
end

main()
