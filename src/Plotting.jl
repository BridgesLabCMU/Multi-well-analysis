function plot_screen_line(upper_FC, lower_FC, upper_final, images_directories, protocol, interval)
    bf_reads = filter(r -> 
        r.action    == "Imaging Read" &&
        r.channel   == "Bright Field", protocol)
    for bf in eachrow(bf_reads)
        bf = bf[1]
        mag = bf.magnification
        data_file = "$(mag)_BF_biomass.csv"
        plot_name = "$(mag)_screen_lines.html"
        results_name = "$(mag)_screen_results_time.csv"
        color1       = "#66c2a5" 
        color2       = "#fc8d62"
        color3       = "#8da0cb"
        color4       = "#e78ac3"
        default_color= "black"

        legend_desc = Dict(
          :cond1 => "High peak biomass",
          :cond2 => "High peak & high final biomass",
          :cond3 => "High final biomass",
          :cond4 => "Low peak biomass"
        )

        category_label = Dict(
          :cond1 => "high peak biomass",
          :cond2 => "high peak and high final biomass",
          :cond3 => "high final biomass",
          :cond4 => "low peak biomass"
        )

        all_series = Vector{Tuple{Vector{Float64}, String}}()
        for (i, folder) in enumerate(images_directories)
            if !isfile(joinpath(folder, "Numerical data", data_file))
                continue
            end
          df = CSV.read(joinpath(folder, "Numerical data", data_file), DataFrame)
          old = names(df)
          new = [
	      let fn = basename(col)
              m = match(r"([^_]+)_.*$", fn)
              m !== nothing ? "Plate $(i): $(m.captures[1])" : col
              end for col in old
          ]
          rename!(df, new)
          for col in new
            push!(all_series, (Float64.(df[!, col]), col))
          end
        end

        peak_vals = [maximum(y) for (y, _) in all_series]
        mean_peak = mean(peak_vals)
	final_vals = [last(y) for (y, _) in all_series]

        traces = GenericTrace[]
        seen   = Dict(k => false for k in keys(legend_desc))
        cats   = Symbol[]

        for (y_vals, label) in all_series
          pmax   = maximum(y_vals)
          pfinal = last(y_vals)

          cat = :none
          if   pmax >  upper_FC*mean_peak
            cat = (pfinal > upper_final*mean_peak ? :cond2 : :cond1)
          elseif pmax <  lower_FC*mean_peak
            cat = :cond4
          elseif pfinal > upper_final*mean_peak
            cat = :cond3
          end
          push!(cats, cat)

          color = cat == :cond1 ? color1 :
                  cat == :cond2 ? color2 :
                  cat == :cond3 ? color3 :
                  cat == :cond4 ? color4 :
                  default_color

          showlegend = false
          name       = ""
          if cat != :none && !seen[cat]
            seen[cat]    = true
            showlegend   = true
            name          = legend_desc[cat]
          end

          push!(traces, PlotlyJS.scatter(
			x = (0:length(y_vals)-1) .* interval ./ 60,
            y          = y_vals,
            mode       = "lines",
            line       = attr(color=color),
            opacity    = 0.5,
            name       = name,
            showlegend = showlegend,
            hovertext  = label,
            hoverinfo  = "text"
          ))
        end

        abn_idx     = findall(c -> c != :none, cats)
        abn_wells   = [all_series[i][2] for i in abn_idx]
        abn_cats    = [category_label[cats[i]] for i in abn_idx]
	abn_peak_norms = peak_vals[abn_idx] ./ mean_peak
	abn_final_norms = final_vals[abn_idx] ./ mean_peak
	df_abnormal = DataFrame(
	  Well      = abn_wells,
	  Category  = abn_cats,
	  NormPeak  = abn_peak_norms,
	  NormFinal = abn_final_norms
	)
        CSV.write(results_name, df_abnormal)

        layout = Layout(
          showlegend = true,
          hovermode  = "closest",
          font       = attr(family="Helvetica, Arial, sans-serif", size=25),
          xaxis      = attr(title="Time (hr)", titlefont=attr(size=25), tickfont=attr(size=25)),
          yaxis      = attr(title="Biofilm biomass (a.u.)", titlefont=attr(size=25), tickfont=attr(size=25))
        )

        fig = PlotlyJS.plot(traces, layout)

        open(plot_name, "w") do io
          PlotlyBase.to_html(io, fig.plot)
        end
    end
end

function plot_screen_scatter(upper_FC, lower_FC, images_directories, protocol)
    bf_reads = filter(r -> 
        r.action    == "Imaging Read" &&
        r.channel   == "Bright Field", protocol)
    for bf in eachrow(bf_reads)
        bf = bf[1]
        mag = bf.magnification
        data_file = "$(mag)_BF_biomass.csv"
        plot_name = "$(mag)_screen_scatter.html"
        results_name = "$(mag)_screen_results_peak.csv"
        high_color     = "#66c2a5"
        low_color      = "#e78ac3"
        default_color  = "black" 
        category_desc = Dict(
            :high => "High peak biomass",
            :low  => "Low peak biomass"
        )
        
        category_label = Dict(
            :high => "high peak biomass",
            :low  => "low peak biomass",
            :none => "normal"
        )

        peaks  = Float64[]
	finals = Float64[]
        labels = String[]

        for (i, folder) in enumerate(images_directories)
            if !isfile(joinpath(folder, "Numerical data", data_file))
                continue
            end
            df = CSV.read(joinpath(folder, "Numerical data", data_file), DataFrame)
            old = names(df)
            new = [ 
	      let fn = basename(col)
              m = match(r"([^_]+)_.*$", fn)
              m !== nothing ? "Plate $(i): $(m.captures[1])" : col
              end for col in old
            ]
            rename!(df, new)
            for col in new
		vals = Float64.(df[!, col])
                push!(peaks, maximum(vals))
		push!(finals, last(vals))
                push!(labels, col)
            end
        end

        mean_peak = mean(peaks)
        cats = Symbol[]
        for p in peaks
            if   p > upper_FC*mean_peak
                push!(cats, :high)
            elseif p < lower_FC*mean_peak
                push!(cats, :low)
            else
                push!(cats, :none)
            end
        end
        abn_idx    = findall(c -> c != :none, cats)
        abn_wells  = labels[abn_idx]
        abn_cats   = cats[abn_idx]
        abn_labels = [ category_label[c] for c in abn_cats ]
	abn_peak_norms  = peaks[abn_idx]   ./ mean_peak
	abn_final_norms = finals[abn_idx]  ./ mean_peak

        df_abnormal = DataFrame(
            Well     = abn_wells,
            Category = abn_labels,
	    NormPeak = abn_peak_norms,
	    NormFinal = abn_final_norms
        )

        CSV.write(results_name, df_abnormal)

        idx           = sortperm(peaks; rev=true)
        ranked_peaks  = peaks[idx]
        ranked_labels = labels[idx]
        ranked_cats   = cats[idx]
        ranks         = 1:length(ranked_peaks)

        hx, hy, htext = Int[], Float64[], String[]
        lx, ly, ltext = Int[], Float64[], String[]
        nx, ny, ntext = Int[], Float64[], String[]

        for (r, p, lbl, cat) in zip(ranks, ranked_peaks, ranked_labels, ranked_cats)
            if cat == :high
                push!(hx, r); push!(hy, p); push!(htext, lbl)
            elseif cat == :low
                push!(lx, r); push!(ly, p); push!(ltext, lbl)
            else
                push!(nx, r); push!(ny, p); push!(ntext, lbl)
            end
        end

        traces = GenericTrace[]

        push!(traces, PlotlyJS.scatter(
            x         = hx, y = hy,
            mode      = "markers",
            marker    = attr(color=high_color),
            opacity   = 0.5,
            name      = category_desc[:high],
            hovertext = htext,
            hoverinfo = "text"
        ))

        push!(traces, PlotlyJS.scatter(
            x         = lx, y = ly,
            mode      = "markers",
            marker    = attr(color=low_color),
            opacity   = 0.5,
            name      = category_desc[:low],
            hovertext = ltext,
            hoverinfo = "text"
        ))

        push!(traces, PlotlyJS.scatter(
            x         = nx, y = ny,
            mode      = "markers",
            marker    = attr(color=default_color),
            opacity   = 0.5,
            showlegend= false,
            hovertext = ntext,
            hoverinfo = "text"
        ))

        layout = Layout(
            showlegend = true,
            hovermode  = "closest",
            font       = attr(family="Helvetica, Arial, sans-serif", size=25),
            xaxis      = attr(
              title     = "Rank",
              titlefont = attr(size=25),
              tickfont  = attr(size=25),
            ),
            yaxis      = attr(
              title     = "Peak biofilm biomass (a.u.)",
              titlefont = attr(size=25),
              tickfont  = attr(size=25),
            )
        )

        fig = PlotlyJS.plot(traces, layout)

        open(plot_name, "w") do io
            PlotlyBase.to_html(io, fig.plot)
        end
    end
end

function plot_all_line(images_directories, protocol, interval)
    bf_reads = filter(r -> 
        r.action    == "Imaging Read" &&
        r.channel   == "Bright Field", protocol)
    for bf in eachrow(bf_reads)
        bf = bf[1]
        mag = bf.magnification
        data_file = "$(mag)_BF_biomass.csv"
        plot_name = "$(mag)_all_lines.html"
        all_series = Vector{Tuple{Vector{<:Real}, String}}()

        for (i, folder) in enumerate(images_directories)
            csv_path = joinpath(folder, "Numerical data", data_file)
            if !isfile(csv_path)
                continue
            end
            df = CSV.read(csv_path, DataFrame)

            old_names = names(df)
            new_names = [
	      let fn = basename(col)
              m = match(r"([^_]+)_.*$", fn)
              m !== nothing ? "Plate $(i): $(m.captures[1])" : col
              end for col in old_names
            ]
            rename!(df, new_names)

            for col in names(df)
                push!(all_series, (df[!, col], col))
            end
        end

        traces = GenericTrace[]
        for (y_vals, label) in all_series
            push!(traces, PlotlyJS.scatter(
				x = (0:length(y_vals)-1) .* interval ./ 60,
                y       = y_vals,
                mode    = "lines",
                name    = label,
                line    = attr(color="black"),
                opacity = 0.5,
            ))
        end

        layout = Layout(
          showlegend = false,
          hovermode  = "closest",
          font       = attr(family="Helvetica, Arial, sans-serif", size=25),
          xaxis = attr(
            title     = "Time (hr)",
            titlefont = attr(family="Helvetica, Arial, sans-serif", size=25),
            tickfont  = attr(family="Helvetica, Arial, sans-serif", size=25),
          ),
          yaxis = attr(
            title     = "Biofilm biomass (a.u.)",
            titlefont = attr(family="Helvetica, Arial, sans-serif", size=25),
            tickfont  = attr(family="Helvetica, Arial, sans-serif", size=25),
          ),
        )

        fig = PlotlyJS.plot(traces, layout)

        open(plot_name, "w") do io
            PlotlyBase.to_html(io, fig.plot)
        end
    end
end

function plot_all_scatter(images_directories, protocol)
    bf_reads = filter(r -> 
        r.action    == "Imaging Read" &&
        r.channel   == "Bright Field", protocol)
    for bf in eachrow(bf_reads)
        bf = bf[1]
        mag = bf.magnification
        data_file = "$(mag)_BF_biomass.csv"
        plot_name = "$(mag)_all_scatter.html"
        all_peaks = Vector{Tuple{Float64,String}}()

        for (i, folder) in enumerate(images_directories)
            csv_path = joinpath(folder, "Numerical data", data_file)
            if !isfile(csv_path)
                continue
            end
            df = CSV.read(csv_path, DataFrame)

            old_names = names(df)
            new_names = [
	      let fn = basename(col)
              m = match(r"([^_]+)_.*$", fn)
              m !== nothing ? "Plate $(i): $(m.captures[1])" : col
              end for col in old_names
            ]
            rename!(df, new_names)

            for col in new_names
                peak = maximum(df[!, col])
                push!(all_peaks, (peak, col))
            end
        end

        values = first.(all_peaks)
        labels = last.(all_peaks)
        idx = sortperm(values; rev=true)
        ranked_peaks  = values[idx]
        ranked_labels = labels[idx]
        ranks = 1:length(ranked_peaks)

        trace = PlotlyJS.scatter(
            x         = ranks,
            y         = ranked_peaks,
            mode      = "markers",
            marker    = attr(color="black"),
            opacity   = 0.5,
            hovertext = ranked_labels,
            hoverinfo = "text"
        )

        layout = Layout(
          showlegend = false,
          hovermode  = "closest",
          font       = attr(family="Helvetica, Arial, sans-serif", size=25),
          xaxis = attr(
            title     = "Rank",
            titlefont = attr(family="Helvetica, Arial, sans-serif", size=25),
            tickfont  = attr(family="Helvetica, Arial, sans-serif", size=25),
          ),
          yaxis = attr(
            title     = "Peak biofilm biomass (a.u.)",
            titlefont = attr(family="Helvetica, Arial, sans-serif", size=25),
            tickfont  = attr(family="Helvetica, Arial, sans-serif", size=25),
          ),
        )

        fig = PlotlyJS.plot(trace, layout)
        open(plot_name, "w") do io
            PlotlyBase.to_html(io, fig.plot)
        end
    end
end

function plot_wells_line(selected_wells, images_directories, protocol, interval)
    bf_reads = filter(r -> 
        r.action    == "Imaging Read" &&
        r.channel   == "Bright Field", protocol)
    for bf in eachrow(bf_reads)
        bf = bf[1]
        mag = bf.magnification
        data_file = "$(mag)_BF_biomass.csv"
        plot_name = "$(mag)_well_lines.html"
        all_series = Vector{Tuple{Vector{<:Real},String}}()

        for (i, folder) in enumerate(images_directories)
            csv_path = joinpath(folder, "Numerical data", data_file)
            if !isfile(csv_path)
                continue
            end
            df = CSV.read(csv_path, DataFrame)

            old_names = names(df)
            new_names = [
	      let fn = basename(col)
              m = match(r"([^_]+)_.*$", fn)
              m !== nothing ? "Plate $(i): $(m.captures[1])" : col
              end for col in old_names
            ]
            rename!(df, new_names)

            for well in selected_wells[i]
                label = "Plate $(i): $(well)"
                if label in names(df)
                    push!(all_series, (df[!, label], label))
                end
            end
        end

        traces = GenericTrace[]  
        for (y_vals, label) in all_series
            push!(traces, PlotlyJS.scatter(
                x = (0:length(y_vals)-1) .* interval ./ 60,
                y       = y_vals,
                mode    = "lines",
                name    = label,
                line    = attr(color="black"),
                opacity = 0.5,
            ))
        end

        layout = Layout(
          showlegend = false,
          hovermode  = "closest",
          font       = attr(family="Helvetica, Arial, sans-serif", size=25),
          xaxis = attr(
            title     = "Time (hr)",
            titlefont = attr(size=25),
            tickfont  = attr(size=25),
          ),
          yaxis = attr(
            title     = "Biofilm biomass (a.u.)",
            titlefont = attr(size=25),
            tickfont  = attr(size=25),
          ),
        )

        fig = PlotlyJS.plot(traces, layout)

        open(plot_name, "w") do io
            PlotlyBase.to_html(io, fig.plot)
        end
    end
end

function plot_wells_scatter(selected_wells, images_directories, protocol)
    bf_reads = filter(r -> 
        r.action    == "Imaging Read" &&
        r.channel   == "Bright Field", protocol)
    for bf in eachrow(bf_reads)
        bf = bf[1]
        mag = bf.magnification
        data_file = "$(mag)_BF_biomass.csv"
        plot_name = "$(mag)_well_scatter.html"
        all_peaks = Vector{Tuple{Float64,String}}()

        for (i, folder) in enumerate(images_directories)
            csv_path = joinpath(folder, "Numerical data", data_file)
            if !isfile(csv_path)
                continue
            end
            df = CSV.read(csv_path, DataFrame)

            old_names = names(df)
            new_names = [
	      let fn = basename(col)
              m = match(r"([^_]+)_.*$", fn)
              m !== nothing ? "Plate $(i): $(m.captures[1])" : col
              end for col in old_names
            ]
            rename!(df, new_names)

            for well in selected_wells[i]
                label = "Plate $(i): $(well)"
                if label in names(df)
                    peak = maximum(df[!, label])
                    push!(all_peaks, (peak, label))
                end
            end
        end

        values = first.(all_peaks)
        labels = last.(all_peaks)
        idx = sortperm(values; rev=true)
        ranked_peaks  = values[idx]
        ranked_labels = labels[idx]
        ranks = 1:length(ranked_peaks)

        trace = PlotlyJS.scatter(
            x         = ranks,
            y         = ranked_peaks,
            mode      = "markers",
            marker    = attr(color="black"),
            opacity   = 0.5,
            hovertext = ranked_labels,
            hoverinfo = "text"
        )

        layout = Layout(
          showlegend = false,
          hovermode  = "closest",
          font       = attr(family="Helvetica, Arial, sans-serif", size=25),
          xaxis = attr(
            title     = "Rank",
            titlefont = attr(family="Helvetica, Arial, sans-serif", size=25),
            tickfont  = attr(family="Helvetica, Arial, sans-serif", size=25),
          ),
          yaxis = attr(
            title     = "Peak biofilm biomass (a.u.)",
            titlefont = attr(family="Helvetica, Arial, sans-serif", size=25),
            tickfont  = attr(family="Helvetica, Arial, sans-serif", size=25),
          ),
        )

        fig = PlotlyJS.plot(trace, layout)
        open(plot_name, "w") do io
            PlotlyBase.to_html(io, fig.plot)
        end
    end
end

function plot_conditions_line(conditions, conditions_dict,
        images_directories, protocol, interval)
    bf_reads = filter(r -> 
        r.action    == "Imaging Read" &&
        r.channel   == "Bright Field", protocol)
    for bf in eachrow(bf_reads)
        bf = bf[1]
        mag = bf.magnification
        data_file = "$(mag)_BF_biomass.csv"
        plot_name = "$(mag)_conditions_line.pdf"
        maxlen = maximum(length.(values(conditions_dict)))
        for v in values(conditions_dict)
            isempty(v) && continue
            for _ in 1:(maxlen - length(v))
                push!(v, similar(v[1], 0))
            end
        end
        dfs = Vector{DataFrame}(undef, length(images_directories))
        for (i, folder) in enumerate(images_directories)
            if !isfile(joinpath(folder, "Numerical data", data_file))
                continue
            end
            df = CSV.read(joinpath(folder, "Numerical data", data_file), DataFrame)
            old = names(df)
            new = [
	      let fn = basename(col)
              m = match(r"([^_]+)_.*$", fn)
              m !== nothing ? "Plate $(i): $(m.captures[1])" : col
              end for col in old
            ]
            rename!(df, new)
            dfs[i] = df
        end

        conds_stats = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
        N = nrow(dfs[1])  

        for key in conditions 
            series_list = Vector{Vector{Float64}}()
            for (i, df) in enumerate(dfs)
                for w in conditions_dict[key][i]
                    lbl = "Plate $(i): $(w)"
                    if lbl ∈ names(df)
                        push!(series_list, Float64.(df[!, lbl]))
                    end
                end
            end
            mat = reduce(hcat, series_list)
            μ = vec(mean(mat; dims=2))                  
            σ = vec(std(mat;  dims=2))                 
            conds_stats[key] = (μ, σ)
        end

        f = Figure(size = (6*72, 3*72))
        ax = CairoMakie.Axis(f[1,1];
            xlabel       = "Time (hr)",
            ylabel       = "Biofilm biomass (a.u.)",
        )

        for (j, key) in enumerate(conditions)
            μ, σ = conds_stats[key]
            x = (0:length(μ)-1) .* interval ./ 60

            poly!(ax,
                vcat(x,  reverse(x)),
                vcat(μ .+ σ, reverse(μ .- σ)),
                color       = (Makie.to_colormap(:Dark2_8)[j], 0.3)
            )

            lines!(ax, x, μ;
                color     = Makie.to_colormap(:Dark2_8)[j],
                linewidth = 2,
                label     = key
            )
        end
        f[1,2] = Legend(f, ax, merge = true, unique = true, framevisible=false, rowgap=0)
        ax.rightspinevisible = false
        ax.topspinevisible = false
        CairoMakie.save(plot_name, f)
    end
end

function plot_conditions_jitter(conditions, conditions_dict,
        images_directories, protocol)
    bf_reads = filter(r -> 
        r.action    == "Imaging Read" &&
        r.channel   == "Bright Field", protocol)
    for bf in eachrow(bf_reads)
        bf = bf[1]
        mag = bf.magnification
        data_file = "$(mag)_BF_biomass.csv"
        plot_name = "$(mag)_conditions_jitter.pdf"
        maxlen = maximum(length.(values(conditions_dict)))
        for v in values(conditions_dict)
            isempty(v) && continue
            for _ in 1:(maxlen - length(v))
                push!(v, similar(v[1], 0))
            end
        end
        dfs = Vector{DataFrame}(undef, length(images_directories))
        for (i, folder) in enumerate(images_directories)
            if !isfile(joinpath(folder, "Numerical data", data_file))
                continue
            end
            df = CSV.read(joinpath(folder, "Numerical data", data_file), DataFrame)
            old = names(df)
            new = [ 
	      let fn = basename(col)
              m = match(r"([^_]+)_.*$", fn)
              m !== nothing ? "Plate $(i): $(m.captures[1])" : col
              end for col in old
		   ]
            rename!(df, new)
            dfs[i] = df
        end

        cond_maxima = Dict{String, Vector{Float64}}()
        for key in conditions 
            vals = Float64[]
            for (i, df) in enumerate(dfs)
                for w in conditions_dict[key][i]
                    lbl = "Plate $(i): $(w)"
                    if lbl ∈ names(df)
                        push!(vals, maximum(Float64.(df[!, lbl])))
                    end
                end
            end
            cond_maxima[key] = vals
        end

        categories = conditions 
        cat_to_idx = Dict(cat => ii for (ii, cat) in enumerate(categories))

        all_inds   = Int[]
        all_peaks  = Float64[]
        for cat in categories
            idx = cat_to_idx[cat]
            for v in cond_maxima[cat]
                push!(all_inds, idx)
                push!(all_peaks, v)
            end
        end

        fig = Figure(size=(3.5*72, 3*72))
        ax  = CairoMakie.Axis(fig[1, 1];
            xlabel        = "",
            ylabel        = "Peak biofilm biomass (a.u.)",
            xticks        = (1:length(categories), categories),
            xticklabelrotation = 45,
        )
        
        means = [ mean(cond_maxima[cat])      for cat in categories]
        mins  = [ minimum(cond_maxima[cat])   for cat in categories]
        maxs  = [ maximum(cond_maxima[cat])   for cat in categories]
        crossbar!(ax, Int.(1:length(categories)), 
                  means, mins, maxs; 
                  color=(:white, 0), midlinecolor=:black)

        beeswarm!(
            ax,
            all_inds,
            all_peaks;
            strokecolor = :black,
            color       = (:white,0),
            strokewidth = 1,
            algorithm   = UniformJitter(),
        )

        ax.rightspinevisible = false
        ax.topspinevisible = false
        ylims!(ax, 0, maximum(all_peaks) * 1.1)
        CairoMakie.save(plot_name, fig)
    end
end

function interval_to_minutes(s)
    s = strip(s)
    m = match(r"Interval\s+(\d+):(\d+):(\d+)", s)
    d, hh, mm = parse.(Int, m.captures)
    return d*24*60 + hh*60 + mm
end

function extract_interval(dir)
	metadata = CSV.read(joinpath(dir,"metadata.csv"), DataFrame;
                        header=false)
    col3 = metadata[!, 3]
    row = findfirst(c -> (c isa AbstractString) && occursin(r"Interval\s+\d+:\d+:\d+", strip(c)),
                                    col3)
    if row === nothing
        error("No “Interval …” entry found in column 3")
    end
    minutes = interval_to_minutes(col3[row])
    return minutes
end

function plotting_main()
    config = JSON.parsefile("experiment_config.json")
    images_directories  = config["images_directory"]
    protocol_path = joinpath(images_directories[1], "protocol.csv")
    protocol = CSV.File(protocol_path, header=true)

    if isfile("plot_options.json")
        plot_options = JSON.parsefile("plot_options.json")
        if haskey(plot_options, "full_plate")
            if plot_options["plot_time_traces"] == "true" && plot_options["plot_peaks"] == "true"
                interval = extract_interval(images_directories[1])
                plot_all_line(images_directories, protocol, interval)
                plot_all_scatter(images_directories, protocol)
            elseif plot_options["plot_time_traces"] == "true"
                interval = extract_interval(images_directories[1])
                plot_all_line(images_directories, protocol, interval)
            elseif plot_options["plot_peaks"] == "true"
                plot_all_scatter(images_directories, protocol)
            end
        elseif haskey(plot_options, "plot_wells")
            plot_wells = plot_options["wells"]
            if plot_options["plot_time_traces"] == "true" && plot_options["plot_peaks"] == "true"
                interval = extract_interval(images_directories[1])
                plot_wells_line(plot_wells, images_directories, protocol, interval)
                plot_wells_scatter(plot_wells, images_directories, protocol)    
            elseif plot_options["plot_time_traces"] == "true"
                interval = extract_interval(images_directories[1])
                plot_wells_line(plot_wells, images_directories, protocol, interval)
            elseif plot_options["plot_peaks"] == "true"
                plot_wells_scatter(plot_wells, images_directories, protocol)
            end
        elseif haskey(plot_options, "screen")
            upper_FC    = plot_options["upper_FC"]
            lower_FC    = plot_options["lower_FC"]
            upper_final = plot_options["upper_final"]
            if plot_options["plot_time_traces"] == "true" && plot_options["plot_peaks"] == "true"
                interval = extract_interval(images_directories[1])
                plot_screen_line(upper_FC, lower_FC, upper_final, images_directories, protocol, interval)
                plot_screen_scatter(upper_FC, lower_FC, images_directories, protocol)
            elseif plot_options["plot_time_traces"] == "true"
                interval = extract_interval(images_directories[1])
                plot_screen_line(upper_FC, lower_FC, upper_final, images_directories, protocol, interval)
            elseif plot_options["plot_peaks"] == "true"
                plot_screen_scatter(upper_FC, lower_FC, images_directories, protocol)
            end
        elseif haskey(plot_options, "plot_conditions")
            plot_conditions = plot_options["conditions"]
            conditions_dict = config["conditions"]
            if plot_options["plot_time_traces"] == "true" && plot_options["plot_peaks"] == "true"
                interval = extract_interval(images_directories[1])
                plot_conditions_line(plot_conditions, conditions_dict, images_directories, protocol, interval)
                plot_conditions_jitter(plot_conditions, conditions_dict, images_directories, protocol)
            elseif plot_options["plot_time_traces"] == "true"
                interval = extract_interval(images_directories[1])
                plot_conditions_line(plot_conditions, conditions_dict, images_directories, protocol, interval)
            elseif plot_options["plot_peaks"] == "true"
                plot_conditions_jitter(plot_conditions, conditions_dict, images_directories, protocol)
            end
        end
    end
end
