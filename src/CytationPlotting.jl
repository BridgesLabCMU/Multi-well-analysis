using JSON
using Plots, LaTeXStrings, StatsPlots
using DataFrames, CSV
using Statistics
using HypothesisTests: UnequalVarianceTTest, pvalue
using LsqFit
using Colors: JULIA_LOGO_COLORS
using CategoricalArrays
    
gr()

function convert_latex(string)
    latex_to_unicode = Dict(
        "\$\\alpha\$" => "α",
        "\$\\beta\$" => "β",
        "\$\\gamma\$" => "γ",
        "\$\\delta\$" => "δ",
        "\$\\Delta\$" => "Δ",
        "\$\\epsilon\$" => "ε",
        "\$\\zeta\$" => "ζ",
        "\$\\eta\$" => "η",
        "\$\\theta\$" => "θ",
        "\$\\iota\$" => "ι",
        "\$\\kappa\$" => "κ",
        "\$\\lambda\$" => "λ",
        "\$\\mu\$" => "μ",
        "\$\\nu\$" => "ν",
        "\$\\xi\$" => "ξ",
        "\$\\omicron\$" => "ο",
        "\$\\pi\$" => "π",
        "\$\\rho\$" => "ρ",
        "\$\\sigma\$" => "σ",
        "\$\\tau\$" => "τ",
        "\$\\upsilon\$" => "υ",
        "\$\\phi\$" => "φ",
        "\$\\chi\$" => "χ",
        "\$\\psi\$" => "ψ",
        "\$\\omega\$" => "ω"
    )
    for (latex, unicode) in latex_to_unicode
        string = replace(string, latex => unicode)
    end
    return string
end

propdiv(a,b,c,d) = sqrt.((a./b).^2 .+ (c./d).^2)

function dose_response(conditions, plot_conditions, 
                      plot_normalization, normalization_method, plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, nums, denoms, 
                      dose_concs, plot_size, plots_directory)
    data = [combine(df, names(df) .=> maximum .=> names(df)) for df in data]
    x = dose_concs
    y = []
    error = []
    
    if yscale == "linear"
        yscale = :identity
    elseif yscale == "log"
        yscale = :log10
    end

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
    plt = scatter(x, y, yerror=error, yscale=yscale, color="black", mc=:white, label="Data", size=plot_size)
    xbase = collect(range(minimum(x), maximum(x), 100))
    plot!(plt, xbase, model.(xbase, (pstar,)), yscale=yscale, color="black", label="Fit")
    if occursin("\$", plot_xlabel) 
        plot_xlabel = convert_latex(plot_xlabel)
    end
    if occursin("\$", plot_ylabel[1])
        plot_ylabel[1] = convert_latex(plot_ylabel[1])
    end
    if occursin("\$", plot_title)
        plot_title = convert_latex(plot_title)
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
                      plot_filename, data, nums, denoms, 
                      plot_clab, plot_size, plots_directory)

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
        plot_xlabel = convert_latex(plot_xlabel)
    end
    if occursin("\$", plot_ylabel[1])
        plot_ylabel[1] = convert_latex(plot_ylabel[1])
    end
    if occursin("\$", clab_title)
        clab_title = convert_latex(clab_title)
    end
    if occursin("\$", plot_title)
        plot_title = convert_latex(plot_title)
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
                      plot_filename, data, nums, denoms, 
                      xaxis_data, plot_xaxis, plot_size, plots_directory, yscale)
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
                plot_ylabel[j] = convert_latex(plot_ylabel[j])
            end
            if plot_xaxis == "Time"
                xaxis=xaxis_data
            elseif plot_xaxis == "OD"
                xaxis = mean.(eachcol(OD))
            else
                error("Can only plot time or OD on the x-axis.")
            end
            if occursin("\$", condition)
                condition = convert_latex(condition)
            end
            if j == 1
                if yscale[1] == "linear"
                    yaxisscale = :identity
                elseif yscale[1] == "log"
                    yaxisscale = :log10
                end
                stds .= ifelse.(isnan.(stds), 0, stds)
                plot!(p, xaxis, means, yscale=yaxisscale, marker=:circle, ribbon=stds, label=condition, color=colors[j], 
                      ytickfontcolor=colors[j], legend=false)
                ylabel!(p, plot_ylabel[j], yguidefontcolor=colors[j])
            else
                if yscale[1] == "linear"
                    yaxisscale = :identity
                elseif yscale[1] == "log"
                    yaxisscale = :log10
                end
                stds .= ifelse.(isnan.(stds), 0, stds)
                plot!(p_twin, xaxis, means, yscale=yaxisscale, marker=:circle, ribbon=stds, 
                      label=condition, color=colors[j], 
                      ytickfontcolor=colors[j], legend=false)
                ylabel!(p_twin, plot_ylabel[j], yguidefontcolor=colors[j])
            end
        end
    end
    if occursin("\$", plot_xlabel)
        plot_xlabel = convert_latex(plot_xlabel)
    end
    if occursin("\$", plot_title)
        plot_title = convert_latex(plot_title)
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
                      plot_filename, data, nums, denoms, 
                      xaxis_data, plot_xaxis, plot_size, plots_directory, yscale)

    if occursin("\$", plot_xlabel)
        plot_xlabel = convert_latex(plot_xlabel)
    end
    if occursin("\$", plot_ylabel[1])
        plot_ylabel[1] = convert_latex(plot_ylabel[1])
    end
    if occursin("\$", plot_title)
        plot_title = convert_latex(plot_title)
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
    if yscale == "linear"
        yscale = :identity
    elseif yscale == "log"
        yscale = :log10
    end

    for condition in plot_conditions
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
            condition = convert_latex(condition)
        end
        plot!(p, xaxis, means, yscale=yscale, marker=:circle, ribbon=stds, label=condition)
    end
    xlabel!(p, plot_xlabel)
    ylabel!(p, plot_ylabel[1])
    title!(p, plot_title)
    savefig(p, "$plots_directory/$plot_filename"*".svg")
end

function jitter_vals(values; width=0.05)
    return values .+ width .* (rand(length(values)) .- 0.5)
end

function jitter_plot(conditions, plot_conditions, 
                      plot_normalization, normalization_method, plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, default_color, plot_size, plots_directory, ylims, yscale)
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
                   elseif normalization_method == ""
                        append!(values, data[i][!, col])
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
        plot_xlabel = convert_latex(plot_xlabel)
    end
    if occursin("\$", plot_ylabel[1])
        plot_ylabel[1] = convert_latex(plot_ylabel[1])
    end
    if occursin("\$", plot_title)
        plot_title = convert_latex(plot_title)
    end
    for (k,e) in enumerate(unique_cats)
        if occursin("\$", e)
            unique_cats[k] = convert_latex(e)
        end
    end
    if yscale == "linear"
        yscale = :identity
    elseif yscale == "log"
        yscale = :log10
    end
    if ylims == "default"
        p = scatter(x_vals, values, group=categories, color=default_color, 
                    markerstrokecolor=default_color, alpha=0.4, 
                    xrotation=45, yscale=yscale, 
                    xlims=(x_min, x_max), size=plot_size)
        boxplot!(p, box_x, values, color=default_color, yscale=yscale, linecolor=default_color, 
                 markerstrokecolor=default_color, leg=false, outliers=false, 
                 fillalpha=0.1, linewidth=1.5)
    else
        for (k, e) in enumerate(ylims)
            if e == "nothing"
                if k == 1
                    ylims[k] = -Inf
                else
                    ylims[k] = Inf
                end
            end
        end
        ylims = Tuple(ylims)
        p = scatter(x_vals, values, group=categories, color=default_color, 
                    markerstrokecolor=default_color, yscale=yscale, alpha=0.4, 
                    xrotation=45, 
                    xlims=(x_min, x_max), ylims=ylims, size=plot_size)
        boxplot!(p, box_x, values, color=default_color, yscale=yscale, linecolor=default_color, 
                 markerstrokecolor=default_color, leg=false, outliers=false, 
                 fillalpha=0.1, linewidth=1.5)
    end
    if plot_xlabel == ""
        xticks!(1:length(unique_cats), unique_cats)
    else
        xticks!(1:length(unique_cats), string.(plot_xticks))
    end
    xlabel!(p, plot_xlabel)
    ylabel!(p, plot_ylabel[1])
    title!(p, plot_title)
    savefig(p, "$plots_directory/$plot_filename"*".svg")
end

function extract_groups(arr)
    groups = []
    group_dict = Dict{String, String}()
    for str in arr
        group = split(str, " ", limit=2)[1] 
        if occursin("\$", group)
            group = latexstring(convert_latex(group))
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
                      plot_filename, data, default_color, plot_size, plots_directory, ylims, yscale)

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
	ctg = CategoricalArray(groups_plot)
	levels!(ctg, unique(groups_plot))

    if occursin("\$", plot_xlabel)
        plot_xlabel = latexstring(convert_latex(plot_xlabel))
    end
    if occursin("\$", plot_ylabel[1])
        plot_ylabel[1] = latexstring(convert_latex(plot_ylabel[1]))
    end
    if occursin("\$", plot_title)
        plot_title = latexstring(convert_latex(plot_title))
    end
    if yscale == "linear"
        yscale = :identity
    elseif yscale == "log"
        yscale = :log10
    end
    if ylims == "default"
        p = plot()
        groupedboxplot!(p, box_x, values, group=ctg, yscale=yscale, 
                 leg=true, outliers=true,  
                 linewidth=1.5)
    else
        for (k, e) in enumerate(ylims)
            if e == "nothing"
                if k == 1
                    ylims[k] = -Inf
                else
                    ylims[k] = Inf
                end
            end
        end
        ylims = Tuple(ylims)
        p = plot()
        groupedboxplot!(p, box_x, values, group=ctg, 
                 leg=true, outliers=true,  
                 linewidth=1.5)
    end
    xlabel!(p, plot_xlabel)
    ylabel!(p, plot_ylabel[1])
    title!(p, plot_title)
    savefig(p, "$plots_directory/$plot_filename"*".svg")
end

function select_data(plot_dtype, lum, OD, BF_imaging, CFP_imaging, 
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
                      texas_red_imaging, CY5_imaging, YFP, CY5, 
                      default_color, dose_concs, plots_directory, ylims, yscale)
    plot_size = Tuple(plot_size)
    if length(plot_dtypes) == 1
        data = select_data(plot_dtypes[1], lum, OD, BF_imaging, CFP_imaging, 
                           YFP_imaging, texas_red_imaging, CY5_imaging, YFP, CY5)
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
                                  CY5_imaging, YFP, CY5)
            denominator = select_data(plot_denominators[1], lum, OD, BF_imaging, 
                                  CFP_imaging, YFP_imaging, texas_red_imaging, 
                                  CY5_imaging, YFP, CY5)
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
                                  CY5_imaging, YFP, CY5)
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
                                      CY5_imaging, YFP, CY5)
                denominator = select_data(plot_denominators[1], lum, OD, BF_imaging, 
                                      CFP_imaging, YFP_imaging, texas_red_imaging, 
                                      CY5_imaging, YFP, CY5)
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
                                      CY5_imaging, YFP, CY5) 
                    push!(nums, nothing)
                    push!(denoms, nothing)
                else
                    numerator = select_data(plot_numerators[2], lum, OD, BF_imaging, 
                                          CFP_imaging, YFP_imaging, texas_red_imaging, 
                                          CY5_imaging, YFP, CY5)
                    denominator = select_data(plot_denominators[2], lum, OD, BF_imaging, 
                                          CFP_imaging, YFP_imaging, texas_red_imaging, 
                                          CY5_imaging, YFP, CY5)
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
                                      CY5_imaging, YFP, CY5)
                push!(nums, nothing)
                push!(denoms, nothing)
                numerator = select_data(plot_numerators[1], lum, OD, BF_imaging, 
                                      CFP_imaging, YFP_imaging, texas_red_imaging, 
                                      CY5_imaging, YFP, CY5)
                denominator = select_data(plot_denominators[1], lum, OD, BF_imaging, 
                                      CFP_imaging, YFP_imaging, texas_red_imaging, 
                                      CY5_imaging, YFP, CY5)
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
                      plot_filename, data, default_color, plot_size, 
                      plots_directory, ylims, yscale[1])
    elseif plot_type == "grouped jitter"
        grouped_jitter_plot(conditions, plot_conditions, 
                      plot_normalization, normalization_method[1], plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, default_color, plot_size, 
                      plots_directory, ylims, yscale[1])
    elseif plot_type == "line"
        if plot_xaxis == "Time"
            t = range(0,stop=nrow(data[1]) - 1,length=nrow(data[1])) ./ acquisition_frequency 
            line_plot(conditions, plot_conditions, 
                          plot_normalization, normalization_method[1], plot_title, 
                          plot_ylabel, plot_xlabel, 
                          plot_yticks, plot_xticks, 
                          plot_filename, data, nums, denoms, t, plot_xaxis, plot_size,
                          plots_directory, yscale[1])
        elseif plot_xaxis == "OD"
            line_plot(conditions, plot_conditions, 
                          plot_normalization, normalization_method[1], plot_title, 
                          plot_ylabel, plot_xlabel, 
                          plot_yticks, plot_xticks, 
                          plot_filename, data, nums, denoms, OD, plot_xaxis, plot_size,
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
                          plot_filename, data, nums, denoms, t, plot_xaxis, plot_size,
                          plots_directory, yscale)
        elseif plot_xaxis == "OD"
            twin_y(conditions, plot_conditions, 
                          plot_normalization, normalization_method, plot_title, 
                          plot_ylabel, plot_xlabel, 
                          plot_yticks, plot_xticks, 
                          plot_filename, data, nums, denoms, OD, plot_xaxis, plot_size,
                          plots_directory, yscale)
        else
            error("Can only plot time or OD on the x-axis of a two-axis plot.")
        end
    elseif plot_type == "heatmap"
        heatplot(conditions, plot_conditions, 
                      plot_normalization, normalization_method[1], plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, nums, denoms, plot_clab, plot_size,
                      plots_directory)
    elseif plot_type == "dose-response"
        dose_response(conditions, plot_conditions, 
                      plot_normalization, normalization_method[1], plot_title, 
                      plot_ylabel, plot_xlabel, 
                      plot_yticks, plot_xticks, 
                      plot_filename, data, nums, denoms, dose_concs, plot_size,
                      plots_directory, yscale[1])
    else
        error("Plot type not implemented yet!")
    end
end

function main()
    config = JSON.parsefile("experiment_config.json")
    font = config["plot_font"]["plot1"]
    default(fontfamily=font, titlefont = (15, "computer modern"), legendfontsize = 15, 
            guidefont = (15, :black), colorbar_tickfontsize=12, colorbar_titlefontsize=15, tickfont = (12, :black), 
            guide = L"x", linewidth=2, grid=false, formatter = :plain)

    color = config["plot_color"] 
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
    yscale = config["plot_scale"] 
    ylims = config["ylims"]
    parent_directory = length(images_directories) > 0 ? images_directories[1] : bulk_data[1][1:end-4] 
    plots_directory = "$parent_directory/Plots"
    if isdir(plots_directory)
        rm(plots_directory; recursive = true)
    end
    mkdir(plots_directory)

    data_types = vcat(values(plot_dtypes)...)
    lum = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    OD = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    BF_imaging = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    CFP_imaging = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    YFP_imaging = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    texas_red_imaging = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    CY5_imaging = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
    YFP = Array{Union{Nothing, DataFrame}, 1}(undef, nplates)
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
                      texas_red_imaging, CY5_imaging, YFP, CY5, 
                      color[plot_num], dose_concs[plot_num], plots_directory, ylims[plot_num], yscale[plot_num])
    end
end

main()
