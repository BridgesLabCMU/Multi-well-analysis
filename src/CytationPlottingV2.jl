using JSON
using LaTeXStrings
using PythonCall
using PythonPlot
using DataFrames, CSV
using Statistics
using HypothesisTests: UnequalVarianceTTest, pvalue
using LsqFit
using Colors: JULIA_LOGO_COLORS
using CategoricalArrays

sns = pyimport("seaborn")
    
propdiv(a,b,c,d) = sqrt.((a./b).^2 .+ (c./d).^2)

function linux_path(path)
	if occursin("C:", path)
		return replace(path, "C:" => "/mnt/c")
	elseif occursin("//bridgeslab.bio.cmu.edu/data/", path)
		return replace(path, "//bridgeslab.bio.cmu.edu/data" => "/mnt/b")
	elseif occursin("B:", path)
		return replace(path, "B:" => "/mnt/b")
	end
end

function dose_response(conditions, plot_conditions, 
                      plot_normalization, normalization_method, plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, nums, denoms, 
                      dose_concs, plots_directory)
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
            error = (sqrt.((error ./ y).^2 .+ (norm_std / norm_mean)^2) .- 1) .* 100
        elseif normalization_method == "fold-change"
            y = y ./ norm_mean
            error = sqrt.((error ./ y).^2 .+ (norm_std / norm_mean)^2)
        end
    end

    error .= ifelse.(isnan.(error), 0, error)

    @. model(x, p) = p[1] + (x^p[4])*(p[2]-p[1])/(x^p[4]+p[3]^p[4]) 

    p0 = [0.0, maximum(y), 0.01, 2]
    lb = [0.0, 0.0, 0.0, 0.0]
    ub = [Inf, Inf, maximum(x), 20.0]

    fit = curve_fit(model, x, y, p0, lower=lb, upper=ub)
    pstar = coef(fit)
    fig, ax = pyplot.subplots()
	ax.errorbar(x, y, yerr=error, fmt="o", color="black", mfc="white", label="Data")
	xbase = collect(range(minimum(x), maximum(x), 100))
	ax.plot(xbase, model.(xbase, (pstar,)), color="black", label="Fit")

    if occursin("\$", plot_xlabel) 
        plot_xlabel = latexstring(plot_xlabel)
    end
    if occursin("\$", plot_ylabel[1])
        plot_ylabel[1] = latexstring(plot_ylabel[1])
    end
    if occursin("\$", plot_title)
        plot_title = latexstring(plot_title)
    end
    ax.set_xlabel(plot_xlabel)
	ax.set_ylabel(plot_ylabel[1])
	ax.set_title(plot_title)
    ax.set_yscale(yscale)
    pyplot.locator_params(axis='x', min_n_ticks=5)
    pyplot.locator_params(axis='y', min_n_ticks=5)
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    ax.legend(bbox_to_anchor=(1, 1), frameon=false)
    if length(plot_filename) == 0
        savefig("$plots_directory/plot.svg")
    else 
        savefig("$plots_directory/$plot_filename"*".svg")
    end
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
                      plot_filename, default_color, data, nums, denoms, 
                      plot_clab, plots_directory)
    
    default_color = default_color[1] # In this case, it is either "" or a colorscheme
    if default_color == ""
        default_color = "cool"
    end

    data = [combine(df, names(df) .=> maximum .=> names(df)) for df in data]
    if nums != nothing
        nums = [combine(df, names(df) .=> maximum .=> names(df)) for df in nums]
        denoms = [combine(df, names(df) .=> maximum .=> names(df)) for df in denoms]
    end
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
    if nums != nothing
        nums_matrix = zeros(8, 12, length(conditions))
        denoms_matrix = zeros(8, 12, length(conditions))
    end
    for (j, condition) in enumerate(plot_conditions)
        replicate = 0
        for i in eachindex(conditions)
            if condition in keys(conditions[i])
                for well in conditions[i][condition]
                    replicate += 1
                    row, col = well_name_to_position(well)
                    plate_matrix[row, col, i] = data[i][1, Symbol(well)]
                    if nums != nothing
                        nums_matrix[row, col, i] = nums[i][1, Symbol(well)]
                        denoms_matrix[row, col, i] = denoms[i][1, Symbol(well)]
                    end
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
    if nums != nothing
        block_nums = Array{Float64, 3}(undef, nrows, ncols, replicates)
        block_denoms = Array{Float64, 3}(undef, nrows, ncols, replicates)
    end
    for i in 1:replicates 
        min_row, max_row, min_col, max_col = find_block_boundaries(block_wells[:,i])
        block_matrices[:,:,i] = plate_matrix[min_row:max_row, min_col:max_col, plate[i]]
        if nums != nothing
            block_nums[:,:,i] = nums_matrix[min_row:max_row, min_col:max_col, plate[i]]
            block_denoms[:,:,i] = denoms_matrix[min_row:max_row, min_col:max_col, plate[i]]
        end
    end

    block_mean = mean(block_matrices, dims=3)
    if nums != nothing
        block_nums_mean = mean(block_nums, dims=3)
        block_denoms_mean = mean(block_denoms, dims=3)
        block_nums_std = std(block_nums, dims=3)
        block_denoms_std = std(block_denoms, dims=3)
        block_std = block_mean.*propdiv(block_nums_std, block_nums_mean, block_denoms_std, block_denoms_mean)
    else
        block_std = std(block_matrices, dims=3)
    end

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
            block_mean = (block_mean ./ norm_mean .- 1) .* 100
            block_std = (sqrt.((block_std ./ block_mean).^2 .+ (norm_std / norm_mean)^2) .- 1) .* 100
            clab_title = "%"
        elseif normalization_method == "fold-change"
            block_mean = block_mean ./ norm_mean
            block_std = sqrt.((block_std ./ block_mean).^2 .+ (norm_std / norm_mean)^2)
            clab_title = "Fold-change"
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
    fig, ax = pyplot.subplots()
	im = ax.imshow(block_mean[:,:,1], cmap=default_color, origin="upper", aspect="auto")
	ax.set_xticks(0:length(plot_xticks)-1)
	ax.set_xticklabels(plot_xticks, rotation=45, ha="right")
	ax.set_yticks(0:length(plot_yticks)-1)
	ax.set_yticklabels(reverse(plot_yticks))
	ax.set_xlabel(plot_xlabel)
	ax.set_ylabel(plot_ylabel[1])
	ax.set_title(plot_title)
	cbar = fig.colorbar(im)
	cbar.ax.set_ylabel(clab_title)
	pyplot.tight_layout()
	savefig(plots_directory*"/"*"$plot_filename" * "_mean.svg")
	fig, ax = pyplot.subplots()
	im = ax.imshow(block_std[:,:,1], cmap=default_color, origin="upper", aspect="auto")
	ax.set_xticks(0:length(plot_xticks)-1)
	ax.set_xticklabels(plot_xticks, rotation=45, ha="right")
	ax.set_yticks(0:length(plot_yticks)-1)
	ax.set_yticklabels(reverse(plot_yticks))
	ax.set_xlabel(plot_xlabel)
	ax.set_ylabel(plot_ylabel[1])
	ax.set_title(plot_title)
	cbar = fig.colorbar(im)
	cbar.ax.set_ylabel(clab_title)
	pyplot.tight_layout()
	savefig(plots_directory*"/"*"$plot_filename" * "_std.svg")
end

function twin_y(conditions, plot_conditions,
                      plot_normalization, normalization_method, plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, default_color, data, nums, denoms, 
                      xaxis_data, plot_xaxis, plots_directory, yscale)
    if default_color[1] == ""
        colors = ["#de8466", "#7dc6e3"]
    else
        colors = default_color
    end
    fig, ax1 = pyplot.subplots()
    ax2 = ax1.twinx()
    
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
            if nums[j] != nothing
                nums_data = DataFrame()
                denoms_data = DataFrame()
            end
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
                    if nums[j] != nothing
                        subset = nums[j][i][!, conditions[i][condition]]
                        subset = DataFrame(collect.(eachrow(subset)), :auto)
                        append!(nums_data, subset)
                        subset = denoms[j][i][!, conditions[i][condition]]
                        subset = DataFrame(collect.(eachrow(subset)), :auto)
                        append!(denoms_data, subset)
                    end
                end
            end
            means = mean.(eachcol(condition_data))
            if nums[j] != nothing
                stds = means.*propdiv(std.(eachcol(nums_data)), mean.(eachcol(nums_data)), 
                                     std.(eachcol(denoms_data)), mean.(eachcol(denoms_data)))
            else
                stds = std.(eachcol(condition_data))
            end
            if plot_normalization != ""
                norm_method_idx = min(length(normalization_method), j)
                if normalization_method[norm_method_idx] == "percent"
                    means = (means ./ max_norm .- 1) .* 100
                    stds = (sqrt.((stds ./ means).^2 .+ (max_std / max_norm)^2) .- 1) .* 100
                elseif normalization_method[norm_method_idx] == "fold-change"
                    means = means ./ max_norm
                    stds = sqrt.((stds ./ means).^2 .+ (max_std / max_norm)^2)
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
            if occursin("\$", condition)
                condition = latexstring(condition)
            end
            if j == 1
                stds .= ifelse.(isnan.(stds), 0, stds)
                ax1.plot(xaxis, means, marker="o", markeredgewidth=1, markeredgecolor="black", color=colors[j])
                ax1.tick_params(axis="y", colors=colors[j])
                ax1.set_yscale(yscale[1])
                ax1.set_ylabel(plot_ylabel[j], color=colors[j])
                ax1.fill_between(xaxis, means .- stds, means .+ stds, color=colors[j], alpha=0.3)
                ax1.locator_params(axis='y', min_n_ticks=5)
            else
                stds .= ifelse.(isnan.(stds), 0, stds)
                ax2.plot(xaxis, means, marker="o", markeredgewidth=1, markeredgecolor="black", color=colors[j])
                ax2.tick_params(axis="y", colors=colors[j])
                ax2.set_yscale(yscale[2])
                ax2.set_ylabel(plot_ylabel[j], color=colors[j])
                ax2.fill_between(xaxis, means .- stds, means .+ stds, color=colors[j], alpha=0.3)
                ax2.locator_params(axis='y', min_n_ticks=5)
            end
        end
    end
    if occursin("\$", plot_xlabel)
        plot_xlabel = latexstring(plot_xlabel)
    end
    if occursin("\$", plot_title)
        plot_title = latexstring(plot_title)
    end
    ax1.set_xlabel(plot_xlabel)
    ax1.set_title(plot_title)
    pyplot.locator_params(axis='x', min_n_ticks=5)
    ax1.spines["top"].set_visible(false)
    if length(plot_filename) == 0
        savefig("$plots_directory/plot.svg")
    else 
        savefig("$plots_directory/$plot_filename"*".svg")
    end
end

function line_plot(conditions, plot_conditions, 
                      plot_normalization, normalization_method, plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, default_color, data, nums, denoms, 
                      xaxis_data, plot_xaxis, plots_directory, yscale)
    
    fig, ax = pyplot.subplots()
    
    if length(default_color) == 1
        default_color = default_color[1] # In this case, it is either "" or a colorscheme or a color
        if default_color == ""
            ax.set_prop_cycle(color=pyplot.get_cmap("Set2").colors)
        elseif !(occursin("#", default_color))
            ax.set_prop_cycle(color=default_color)
        end
    end

    if occursin("\$", plot_xlabel)
        plot_xlabel = latexstring(plot_xlabel)
    end
    if occursin("\$", plot_ylabel[1])
        plot_ylabel[1] = latexstring(plot_ylabel[1])
    end
    if occursin("\$", plot_title)
        plot_title = latexstring(plot_title)
    end

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

    for (k, condition) in enumerate(plot_conditions)
        if plot_xaxis == "OD"
            OD = DataFrame()
        end
        condition_data = DataFrame() 
        numerators = DataFrame() 
        denominators = DataFrame() 
        for i in eachindex(conditions)
            if condition in keys(conditions[i])
                subset = data[i][!, conditions[i][condition]]
                subset = DataFrame(collect.(eachrow(subset)), :auto)
                if ncol(condition_data) > ncol(subset)
                    condition_data = condition_data[:, 1:ncol(subset)]
                elseif (ncol(subset) > ncol(condition_data)) && (ncol(condition_data) > 0)
                    subset = subset[:, 1:ncol(condition_data)]
                end
                append!(condition_data, subset)
                if plot_xaxis == "OD"
                    subset = xaxis_data[i][!, conditions[i][condition]]
                    subset = DataFrame(collect.(eachrow(subset)), :auto)
                    append!(OD, subset)
                end
                if nums != nothing
                    subset = nums[i][!, conditions[i][condition]]
                    subset = DataFrame(collect.(eachrow(subset)), :auto)
                    append!(numerators, subset)
                    subset = denoms[i][!, conditions[i][condition]]
                    subset = DataFrame(collect.(eachrow(subset)), :auto)
                    append!(denominators, subset)
                end
            end
        end
        means = mean.(eachcol(condition_data))
        stds = nums==nothing ? std.(eachcol(condition_data)) : means.*propdiv(std.(eachcol(numerators)), mean.(eachcol(numerators)), 
                                                                     std.(eachcol(denominators)), mean.(eachcol(denominators)))
        if plot_normalization != ""
            if normalization_method == "percent"
                means = (means ./ max_norm .- 1) .* 100
                stds = (sqrt.((stds ./ means).^2 .+ (max_std / max_norm)^2) .- 1) .* 100
            elseif normalization_method == "fold-change"
                means = means ./ max_norm
                stds = sqrt.((stds ./ means).^2 .+ (max_std / max_norm)^2)
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
        if occursin("\$", condition)
            condition = latexstring(condition)
        end
        stds .= ifelse.(isnan.(stds), 0, stds)
        if isa(default_color, Array)
            ax.plot(xaxis, means, marker="o", markeredgewidth=1, markeredgecolor="black", color=default_color[k], label=condition)
            ax.fill_between(xaxis, means .- stds, means .+ stds, color=default_color[k], alpha=0.3)
        elseif occursin("#", default_color)
            ax.plot(xaxis, means, marker="o", markeredgewidth=1, markeredgecolor="black", color=default_color, label=condition)
            ax.fill_between(xaxis, means .- stds, means .+ stds, color=default_color, alpha=0.3)
        else
            ax.plot(xaxis, means, marker="o", markeredgewidth=1, markeredgecolor="black", label=condition)
            ax.fill_between(xaxis, means .- stds, means .+ stds, alpha=0.3)
        end
    end
    ax.set_ylabel(plot_ylabel[1])
    ax.set_xlabel(plot_xlabel)
    ax.set_title(plot_title)
    ax.set_yscale(yscale)
    pyplot.locator_params(axis='x', min_n_ticks=5)
    pyplot.locator_params(axis='y', min_n_ticks=5)
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    ax.legend(bbox_to_anchor=(1, 1), frameon=false)
    if length(plot_filename) == 0
        savefig("$plots_directory/plot.svg")
    else 
        savefig("$plots_directory/$plot_filename"*".svg")
    end
end

function jitter_plot(conditions, plot_conditions, 
                      plot_normalization, normalization_method, plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, default_color, plots_directory, ylims, yscale)
    if length(default_color) == 1
        default_color = default_color[1] # In this case, it is either "" or a colorscheme or a color
        if default_color == ""
            default_color = "#589fc4"
        end
    end

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

    for condition in plot_conditions
        for i in eachindex(conditions)
            if condition in keys(conditions[i])
                for col in conditions[i][condition]
                    append!(categories, repeat([condition], length(data[i][!, col])))
                    if normalization_method == "percent"
                        append!(values, (data[i][!, col] ./ norm_mean .- 1) .* 100)
                   elseif normalization_method == "fold-change"
                        append!(values, data[i][!, col] ./ norm_mean)
                   elseif normalization_method == ""
                        append!(values, data[i][!, col])
                    end
                end
            end
        end
    end

    unique_cats = unique(categories)

    if occursin("\$", plot_xlabel)
        plot_xlabel = latexstring(plot_xlabel)
    end
    if occursin("\$", plot_ylabel[1])
        plot_ylabel[1] = latexstring(plot_ylabel[1])
    end
    if occursin("\$", plot_title)
        plot_title = latexstring(plot_title)
    end
    for (k,e) in enumerate(unique_cats)
        if occursin("\$", e)
            unique_cats[k] = latexstring(e)
        end
    end

    fig, ax1 = pyplot.subplots()

    if ylims == "default"
        if isa(default_color, Array)
            sns.boxplot(x=categories, y=values, ax=ax1, palette=default_color, showfliers=false)
            sns.stripplot(x=categories, y=values, ax=ax1, palette=default_color, dodge=false, linewidth=1, alpha=0.8, jitter=false)
        elseif occursin("#", default_color)
            sns.boxplot(x=categories, y=values, ax=ax1, palette=(default_color,), showfliers=false)
            sns.stripplot(x=categories, y=values, ax=ax1, palette=(default_color,), dodge=false, linewidth=1, alpha=0.8, jitter=false)
        elseif default_color == ""
            sns.boxplot(x=categories, y=values, ax=ax1, palette=(default_color,), showfliers=false)
            sns.stripplot(x=categories, y=values, ax=ax1, palette=(default_color,), dodge=false, linewidth=1, alpha=0.8, jitter=false)
        else
            sns.boxplot(x=categories, y=values, ax=ax1, palette=default_color, showfliers=false)
            sns.stripplot(x=categories, y=values, ax=ax1, palette=default_color, dodge=false, linewidth=1, alpha=0.8, jitter=false)
        end
    else
        for (k, e) in enumerate(ylims)
            if e == "nothing"
                if k == 1
                    ylims[k] = nothing 
                else
                    ylims[k] = nothing 
                end
            end
        end
        ylims = Tuple(ylims)
        if isa(default_color, Array)
            sns.boxplot(x=categories, y=values, ax=ax1, palette=default_color, showfliers=false)
            sns.stripplot(x=categories, y=values, ax=ax1, palette=default_color, dodge=false, linewidth=1, alpha=0.8, jitter=false)
        elseif occursin("#", default_color)
            sns.boxplot(x=categories, y=values, ax=ax1, palette=(default_color,), showfliers=false)
            sns.stripplot(x=categories, y=values, ax=ax1, palette=(default_color,), dodge=false, linewidth=1, alpha=0.8, jitter=false)
        elseif default_color == ""
            sns.boxplot(x=categories, y=values, ax=ax1, palette=(default_color,), showfliers=false)
            sns.stripplot(x=categories, y=values, ax=ax1, palette=(default_color,), dodge=false, linewidth=1, alpha=0.8, jitter=false)
        else
            sns.boxplot(x=categories, y=values, ax=ax1, palette=default_color, showfliers=false)
            sns.stripplot(x=categories, y=values, ax=ax1, palette=default_color, dodge=false, linewidth=1, alpha=0.8, jitter=false)
        end
        ax1.set_ylim(ylims)
    end

    if plot_xlabel != ""
        pyplot.xticks(0:length(unique_cats)-1, string.(plot_xticks))
    end

	ax1.set_xticklabels(unique_cats, rotation=45)
    ax1.set_yscale(yscale)
    ax1.set_xlabel(plot_xlabel)
    ax1.set_ylabel(plot_ylabel[1])
    ax1.set_title(plot_title)
    pyplot.locator_params(axis='y', min_n_ticks=5)
    ax1.spines["right"].set_visible(false)
    ax1.spines["top"].set_visible(false)
    if length(plot_filename) == 0
        savefig("$plots_directory/plot.svg")
    else 
        savefig("$plots_directory/$plot_filename"*".svg")
    end
end

function tex_split(str::String)
    m = match(r"^((?:[^$\s]|\$[^$]*\$)*)\s", str)
    return m !== nothing ? m.captures[1] : str
end

function extract_groups(arr)
    groups = []
    group_dict = Dict{String, String}()
    for str in arr
        group = tex_split(str) 
        if occursin("\$", group)
            group = latexstring(group)
        end
        push!(groups, group)
        group_dict[str] = group
    end
	groups = unique(groups) 	
    return group_dict, collect(groups)
end

function Base.unique(ctg::CategoricalArray)
    l = levels(ctg)
    newctg = CategoricalArray(l)
    levels!(newctg, l)
end

function grouped_jitter_plot(conditions, plot_conditions, 
                      plot_normalization, normalization_method, plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, default_color, plots_directory, ylims, yscale)
    
    if length(default_color) == 1
        default_color = default_color[1] # In this case, it is either "" or a colorscheme
    end

    if plot_xticks == []
        error("Can't plot a grouped jitter plot without specifying x-ticks.")
    end

    group_dict, groups = extract_groups(plot_conditions)

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

    values = Float64[]
    box_x = []
    groups_plot = []

    x_val = 1 
    last_group = nothing
    for condition in plot_conditions
        for i in eachindex(conditions)
            if condition in keys(conditions[i])
                group = group_dict[condition]
                if last_group == group 
                    x_val += 1
                else
                    x_val = 1
                end
                for col in conditions[i][condition]
                    push!(groups_plot, group)
                    if normalization_method == "percent"
                        append!(values, (data[i][!, col] ./ norm_mean .- 1) .* 100)
                    elseif normalization_method == "fold-change"
                        append!(values, data[i][!, col] ./ norm_mean)
                    elseif normalization_method == ""
                        append!(values, data[i][!, col])
                    end
                    push!(box_x, x_val)
                end
                last_group = group
            end
        end
    end
    
	unique_vals = unique(box_x)
	val_to_str = Dict(zip(unique_vals, plot_xticks))
    box_x = [val_to_str[val] for val in box_x]
    if occursin("\$", plot_xlabel)
        plot_xlabel = latexstring(plot_xlabel)
    end
    if occursin("\$", plot_ylabel[1])
        plot_ylabel[1] = latexstring(plot_ylabel[1])
    end
    if occursin("\$", plot_title)
        plot_title = latexstring(plot_title)
    end
    fig, ax1 = pyplot.subplots()
    if isa(default_color, Array)
        sns.boxplot(x=box_x, y=values, hue=groups_plot, showfliers=false, palette=default_color)
        sns.stripplot(x=box_x, y=values, hue=groups_plot, dodge=true, alpha=0.8, palette=default_color, linewidth=1, jitter=true, legend=nothing)
    elseif default_color == ""
        sns.boxplot(x=box_x, y=values, hue=groups_plot, showfliers=false, palette="Set2")
        sns.stripplot(x=box_x, y=values, hue=groups_plot, dodge=true, alpha=0.8, palette="Set2", linewidth=1, jitter=true, legend=nothing)
    else
        sns.boxplot(x=box_x, y=values, hue=groups_plot, showfliers=false, palette=default_color)
        sns.stripplot(x=box_x, y=values, hue=groups_plot, dodge=true, alpha=0.8, palette=default_color, linewidth=1, jitter=true, legend=nothing)
    end
    if ylims != "default"
        for (k, e) in enumerate(ylims)
            if e == "nothing"
                if k == 1
                    ylims[k] = nothing 
                else
                    ylims[k] = nothing 
                end
            end
        end
        ylims = Tuple(ylims)
        ax1.set_ylim(ylims)
    end
	
    pyplot.xticks(0:length(unique_vals)-1, string.(plot_xticks),rotation=45)
    ax1.set_xlabel(plot_xlabel)
    ax1.set_ylabel(plot_ylabel[1])
    ax1.set_title(plot_title)
    pyplot.locator_params(axis='y', min_n_ticks=5)
    ax1.spines["right"].set_visible(false)
    ax1.spines["top"].set_visible(false)
    ax1.legend(bbox_to_anchor=(1, 1), frameon=false)
    if length(plot_filename) == 0
        savefig("$plots_directory/plot.svg")
    else 
        savefig("$plots_directory/$plot_filename"*".svg")
    end
end

function select_data(plot_dtype, lum, OD, BF_imaging, CFP_imaging, 
                    YFP_imaging, texas_red_imaging, CY5_imaging, GFP, CY5)
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
    elseif plot_dtype == "GFP"
        data = GFP
    elseif plot_dtype == "CY5"
        data = CY5
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
                      lum, OD, BF_imaging, CFP_imaging, YFP_imaging, 
                      texas_red_imaging, CY5_imaging, GFP, CY5, 
                      default_color, dose_concs, plots_directory, ylims, yscale)
    PythonPlot.matplotlib.rcParams["figure.figsize"] = plot_size 
    if length(plot_dtypes) == 1
        data = select_data(plot_dtypes[1], lum, OD, BF_imaging, CFP_imaging, 
                           YFP_imaging, texas_red_imaging, CY5_imaging, GFP, CY5)
        nums = nothing
        denoms = nothing
    elseif length(plot_dtypes) == 2 && plot_type != "two-axis" 
        if length(plot_numerators) == 0
            error("Passed too many data types without specifying numerator and denominator.")
        elseif length(plot_numerators) > 1
            error("Passed too many numerators/denominators.")
        else
            numerator = select_data(plot_numerators[1], lum, OD, BF_imaging, 
                                  CFP_imaging, YFP_imaging, texas_red_imaging, 
                                  CY5_imaging, GFP, CY5)
            denominator = select_data(plot_denominators[1], lum, OD, BF_imaging, 
                                  CFP_imaging, YFP_imaging, texas_red_imaging, 
                                  CY5_imaging, GFP, CY5)
            data = []
            quotient_df = DataFrame()
            for j in 1:length(numerator)
                for col_name in names(numerator[j])
                    quotient_df[!, col_name] = numerator[j][!, col_name] ./ denominator[j][!, col_name]
                end
                quotient_df .= ifelse.(isnan.(quotient_df), 0, quotient_df)
                push!(data, quotient_df)
            end
            nums = numerator
            denoms = denominator 
        end
    elseif length(plot_dtypes) == 2 && plot_type == "two-axis"
        data = Array{Vector{Union{Nothing, DataFrame}}, 1}(undef, 2)
        nums = []
        denoms = []
        for i in 1:2
            data[i] = select_data(plot_dtypes[i], lum, OD, BF_imaging, 
                                  CFP_imaging, YFP_imaging, texas_red_imaging, 
                                  CY5_imaging, GFP, CY5)
            push!(nums, nothing)
            push!(denoms, nothing)
        end
    elseif length(plot_dtypes) > 2 && plot_type == "two-axis"
        if length(plot_numerators) == 0
            error("Passed too many data types without specifying numerator and denominator.")
        elseif length(plot_numerators) > 2
            error("Passed too many numerators/denominators.")
        else
            nums = []
            denoms = [] 
            data = Array{Vector{Union{Nothing, DataFrame}}, 1}(undef, 2)
            plot_dtype1 = plot_dtypes[1]
            if in(plot_dtype1, plot_numerators) || in(plot_dtype1, plot_denominators)  
                numerator = select_data(plot_numerators[1], lum, OD, BF_imaging, 
                                      CFP_imaging, YFP_imaging, texas_red_imaging, 
                                      CY5_imaging, GFP, CY5)
                denominator = select_data(plot_denominators[1], lum, OD, BF_imaging, 
                                      CFP_imaging, YFP_imaging, texas_red_imaging, 
                                      CY5_imaging, GFP, CY5)
                quotient_df = DataFrame()
                quotient_dfs = []
                for j in 1:length(numerator)
                    for col_name in names(numerator[j])
                        quotient_df[!, col_name] = numerator[j][!, col_name] ./ denominator[j][!, col_name]
                    end
                    quotient_df .= ifelse.(isnan.(quotient_df), 0, quotient_df)
                    push!(quotient_dfs, quotient_df)
                end
                data[1] = quotient_dfs
                push!(nums, numerator)
                push!(denoms, denominator)
                if length(plot_dtypes) == 3
                    plot_dtype = filter(x -> x != plot_numerators[1] && x != plot_denominators[1], plot_dtypes)[1]
                    data[2] = select_data(plot_dtype, lum, OD, BF_imaging, 
                                      CFP_imaging, YFP_imaging, texas_red_imaging, 
                                      CY5_imaging, GFP, CY5) 
                    push!(nums, nothing)
                    push!(denoms, nothing)
                else
                    numerator = select_data(plot_numerators[2], lum, OD, BF_imaging, 
                                          CFP_imaging, YFP_imaging, texas_red_imaging, 
                                          CY5_imaging, GFP, CY5)
                    denominator = select_data(plot_denominators[2], lum, OD, BF_imaging, 
                                          CFP_imaging, YFP_imaging, texas_red_imaging, 
                                          CY5_imaging, GFP, CY5)
                    quotient_df = DataFrame()
                    for j in 1:length(numerator)
                        for col_name in names(numerator[j])
                            quotient_df[!, col_name] = numerator[j][!, col_name] ./ denominator[j][!, col_name]
                        end
                        quotient_df .= ifelse.(isnan.(quotient_df), 0, quotient_df)
                        push!(data[2], quotient_df) 
                    end
                    push!(nums, numerator)
                    push!(denoms, denominator)
                end
            else
                data[1] = select_data(plot_dtype1, lum, OD, BF_imaging, 
                                      CFP_imaging, YFP_imaging, texas_red_imaging, 
                                      CY5_imaging, GFP, CY5)
                push!(nums, nothing)
                push!(denoms, nothing)
                numerator = select_data(plot_numerators[1], lum, OD, BF_imaging, 
                                      CFP_imaging, YFP_imaging, texas_red_imaging, 
                                      CY5_imaging, GFP, CY5)
                denominator = select_data(plot_denominators[1], lum, OD, BF_imaging, 
                                      CFP_imaging, YFP_imaging, texas_red_imaging, 
                                      CY5_imaging, GFP, CY5)
                quotient_df = DataFrame()
                for col_name in names(numerator[1])
                    quotient_df[!, col_name] = numerator[1][!, col_name] ./ denominator[1][!, col_name]
                end
                quotient_df .= ifelse.(isnan.(quotient_df), 0, quotient_df)
                data[2] = [quotient_df] 
                push!(nums, numerator)
                push!(denoms, denominator)
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
                      plot_filename, data, default_color, 
                      plots_directory, ylims, yscale[1])
    elseif plot_type == "grouped jitter"
        grouped_jitter_plot(conditions, plot_conditions, 
                      plot_normalization, normalization_method[1], plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, default_color, 
                      plots_directory, ylims, yscale[1])
    elseif plot_type == "line"
        if plot_xaxis == "Time"
            t = range(0,stop=nrow(data[1]) - 1,length=nrow(data[1])) ./ acquisition_frequency 
            line_plot(conditions, plot_conditions, 
                          plot_normalization, normalization_method[1], plot_title, 
                          plot_ylabel, plot_xlabel, 
                          plot_yticks, plot_xticks, 
                          plot_filename, default_color, data, nums, denoms, t, plot_xaxis,
                          plots_directory, yscale[1])
        elseif plot_xaxis == "OD"
            line_plot(conditions, plot_conditions, 
                          plot_normalization, normalization_method[1], plot_title, 
                          plot_ylabel, plot_xlabel, 
                          plot_yticks, plot_xticks, 
                          plot_filename, default_color, data, nums, denoms, OD, plot_xaxis,
                          plots_directory, yscale[1])
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
                          plot_filename, default_color, data, nums, denoms, t, plot_xaxis,
                          plots_directory, yscale)
        elseif plot_xaxis == "OD"
            twin_y(conditions, plot_conditions, 
                          plot_normalization, normalization_method, plot_title, 
                          plot_ylabel, plot_xlabel, 
                          plot_yticks, plot_xticks, 
                          plot_filename, default_color, data, nums, denoms, OD, plot_xaxis,
                          plots_directory, yscale)
        else
            error("Can only plot time or OD on the x-axis of a two-axis plot.")
        end
    elseif plot_type == "heatmap"
        heatplot(conditions, plot_conditions, 
                      plot_normalization, normalization_method[1], plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, default_color, data, nums, denoms, plot_clab,
                      plots_directory)
    elseif plot_type == "dose-response"
        dose_response(conditions, plot_conditions, 
                      plot_normalization, normalization_method[1], plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, nums, denoms, dose_concs,
                      plots_directory, yscale[1])
    else
        error("Plot type not implemented yet!")
    end
end

function main()
    config = JSON.parsefile("experiment_config.json")
    color = config["plot_color"] 
    font = config["plot_font"]["plot1"][1]
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
    plot_yticks = config["plot_yticks"] 
    plot_clabs = config["color_label"] 
    plot_size = config["plot_size"] 
    plot_filenames = config["plot_filenames"] 
    yscale = config["plot_scale"] 
    ylims = config["ylims"]
    images_directories = linux_path.(images_directories)
    bulk_data = linux_path.(bulk_data)
    parent_directory = length(images_directories) > 0 ? images_directories[1] : bulk_data[1][1:end-4] 
    plots_directory = "$parent_directory/Plots"
    if isdir(plots_directory)
        rm(plots_directory; recursive = true)
    end
    mkdir(plots_directory)

    # Font
    PythonPlot.matplotlib.rcParams["lines.linewidth"] = 2.5
    PythonPlot.matplotlib.rcParams["axes.titlesize"] = 20    
    PythonPlot.matplotlib.rcParams["axes.labelsize"] = 20   
    PythonPlot.matplotlib.rcParams["xtick.labelsize"] = 14 
    PythonPlot.matplotlib.rcParams["ytick.labelsize"] = 14
    PythonPlot.matplotlib.rcParams["legend.fontsize"] = 12    
    if font == "Arial"
        PythonPlot.matplotlib.rcParams["font.family"] = "Arial"
        PythonPlot.matplotlib.rcParams["text.usetex"] = false
    elseif font == "Computer modern"
        PythonPlot.matplotlib.rcParams["text.latex.preamble"] = raw"\usepackage{siunitx} \sffamily"
        PythonPlot.matplotlib.rcParams["text.usetex"] = true 
        PythonPlot.matplotlib.rcParams["text.latex.preamble"] = raw"\usepackage{sfmath}" 
    end

    data_types = vcat(values(plot_dtypes)...)
    lum = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    OD = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    BF_imaging = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    CFP_imaging = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    YFP_imaging = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    texas_red_imaging = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    CY5_imaging = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    GFP = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    CY5 = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    for i in 1:nplates
        if i > length(images_directories)
            push!(images_directories, bulk_data[i][1:end-4])
        end
        data_directory = images_directories[i]
        if "lum" in data_types  && isfile("$data_directory/lum.csv")
            lum[i] = CSV.read("$data_directory/lum.csv", DataFrame)
        else
            lum[i] = nothing
        end
        if "OD" in data_types && isfile("$data_directory/OD600.csv")
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
        if "GFP" in data_types && isfile("$data_directory/GFP.csv")
            GFP[i] = CSV.read("$data_directory/GFP.csv", DataFrame)
        else  
            GFP[i] = nothing
        end
        if "CY5" in data_types && isfile("$data_directory/CY5.csv")
            CY5[i] = CSV.read("$data_directory/CY5.csv", DataFrame)
        else  
            CY5[i] = nothing
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
                      lum, OD, BF_imaging, CFP_imaging, YFP_imaging, 
                      texas_red_imaging, CY5_imaging, GFP, CY5, 
                      color[plot_num], dose_concs[plot_num], plots_directory, ylims[plot_num], yscale[plot_num])
    end
end

main()
