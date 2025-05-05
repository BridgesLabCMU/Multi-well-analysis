function show_plot_options(title::String)
    c = Condition()
    result = Ref{Symbol}()
    trace_cb = GtkCheckButton("Plot time traces")
    peak_cb  = GtkCheckButton("Plot peaks")

    win = GtkWindow(title * " Options")
    v   = GtkBox(:v, 10); push!(win, v)
    push!(v, trace_cb, peak_cb)

    h = GtkBox(:h, 6)
    back = GtkButton("Back"); nxt = GtkButton("Next")
    push!(h, back, nxt); push!(v, h)

    signal_connect(back, "clicked") do _
        destroy(win)
        result[] = :back
        notify(c)
    end
    signal_connect(nxt, "clicked") do _
        if !(trace_cb.active || peak_cb.active)
            Gtk4.warn_dialog(() -> nothing,
                "Please select at least one option.", win)
            return
        end
        destroy(win)
        result[] = :next
        notify(c)
    end

    show(win)
    wait(c)
    return (result[], trace_cb.active, peak_cb.active)
end

function show_condition_selector()
    conf = JSON.parse(read("experiment_config.json", String))
    if !haskey(conf, "conditions")
        Gtk4.warn_dialog(() -> nothing,
            "No “conditions” field found in experiment_config.json.\n" *
            "Please run the conditions selector GUI first.",
            nothing)
        return nothing
    end

    names = collect(keys(conf["conditions"]))
    csel  = Condition()
    sel   = Ref{Vector{String}}()
    win   = GtkWindow("Select Conditions")
    vbox  = GtkBox(:v, 6); push!(win, vbox)

    cbs = Dict{String,GtkCheckButton}()
    for nm in names
        cb = GtkCheckButton(nm)
        push!(vbox, cb)
        cbs[nm] = cb
    end

    hbox = GtkBox(:h, 6)
    back = GtkButton("Back"); nxt = GtkButton("Next")
    push!(hbox, back, nxt); push!(vbox, hbox)

    signal_connect(back, "clicked") do _
        destroy(win)
        sel[] = []
        notify(csel)
    end

    signal_connect(nxt, "clicked") do _
        chosen = [ nm for nm in names if cbs[nm].active ]

        if length(chosen) > 8
            Gtk4.warn_dialog(() -> nothing,
                "Plotting >8 conditions on the same plot is not allowed.",
                win)
            return
        end

        if isempty(chosen)
            Gtk4.warn_dialog(() -> nothing,
                "Please select at least one condition.",
                win)
            return
        end

        destroy(win)
        sel[] = chosen
        notify(csel)
    end

    show(win)
    wait(csel)
    return sel[]
end

function show_stage2()
    c2 = Condition()
    res = Ref{Symbol}()
    win = GtkWindow("Plot all wells")
    v   = GtkBox(:v, 10); push!(win, v)

    cb = GtkCheckButton("Plot all wells")
    push!(v, cb)

    h = GtkBox(:h, 6)
    back = GtkButton("Back"); nxt = GtkButton("Next")
    push!(h, back, nxt); push!(v, h)

    signal_connect(back, "clicked") do _
        destroy(win); res[] = :back; notify(c2)
    end
    signal_connect(nxt, "clicked") do _
        destroy(win); res[] = cb.active ? :yes : :no; notify(c2)
    end

    show(win); wait(c2)
    return res[]
end

function show_stage3b()
    c4 = Condition()
    res = Ref{Symbol}()
    win = GtkWindow("Plot Specific Options")
    v   = GtkBox(:v, 10); push!(win, v)

    cond_cb  = GtkCheckButton("Plot specific conditions")
    wells_cb = GtkCheckButton("Plot specific wells")
    push!(v, cond_cb, wells_cb)

    h = GtkBox(:h, 6)
    back = GtkButton("Back"); nxt = GtkButton("Next")
    push!(h, back, nxt); push!(v, h)

    signal_connect(back, "clicked") do _
        destroy(win); res[] = :back; notify(c4)
    end
    signal_connect(nxt, "clicked") do _
        if xor(cond_cb.active, wells_cb.active)
            destroy(win)
            res[] = cond_cb.active ? :conditions : :wells
            notify(c4)
        else
            Gtk4.warn_dialog(() -> nothing,
                "Select exactly one option.", win)
        end
    end

    show(win); wait(c4)
    return res[]
end

function show_stage4()
    # Load user-selected plate geometry
    conf = JSON.parse(read("experiment_config.json", String))
    geom = get(conf, "plate_geometry", "96-well")

    # Map geometry string to rows and columns entirely within this function
    plate_map = Dict(
        "6-well"   => (2, 3),
        "12-well"  => (3, 4),
        "24-well"  => (4, 6),
        "48-well"  => (6, 8),
        "96-well"  => (8, 12),
        "384-well" => (16, 24)
    )
    rows, cols = haskey(plate_map, geom) ? plate_map[geom] : plate_map["96-well"]

    # Build the dialog
    c5 = Condition()
    win = GtkWindow("Select Specific Wells")
    v   = GtkBox(:v, 10)
    push!(win, v)

    lbl = GtkLabel("Plate 1 of 10")
    push!(v, lbl)

    grid = GtkGrid(column_homogeneous=true, row_homogeneous=true,
                          column_spacing=2, row_spacing=2)
    btns = Matrix{GtkToggleButton}(undef, rows, cols)

    for r in 1:rows, c in 1:cols
        well_id = string(Char('A' + r - 1), c)
        b = GtkToggleButton(well_id)
        b.width_request = 40
        b.height_request = 40
        btns[r, c] = b
        grid[c, r] = b
    end
    push!(v, grid)

    h = GtkBox(:h, 6)
    save_btn = GtkButton("Save Plate")
    done_btn = GtkButton("Done")
    push!(h, save_btn, done_btn)
    push!(v, h)

    cur = Ref(Vector{Vector{String}}())

    signal_connect(save_btn, "clicked") do _
        sel = [string(Char('A' + r - 1), c) for r in 1:rows, c in 1:cols if btns[r,c].active]
        if length(cur[]) < 10
            push!(cur[], sel)
            @idle_add begin
                lbl.label = "Plate $(length(cur[])+1) of 10"
                for b in btns; b.active = false; end
                false
            end
        else
            Gtk4.warn_dialog(() -> nothing,
                "Cannot save more than 10 plates.", win)
        end
    end

    signal_connect(done_btn, "clicked") do _
        destroy(win)
        notify(c5)
    end

    show(win)
    wait(c5)
    return cur[]
end

function run_plot_gui(options_path::String="plot_options.json")
    # single loop
    if !Gtk4.GLib.is_loop_running()
        @async Gtk4.GLib.glib_main()
    end

    c1 = Condition(); vals_ref = Ref{Vector{Float64}}()
    win1 = GtkWindow("Screen Configuration"); v1 = GtkBox(:v,10)
    push!(win1,v1)
    screen_cb = GtkCheckButton("Screen"); push!(v1,screen_cb)
    labels = ["Upper FC...","Lower FC...","Upper FC final..."]
    entries = GtkEntry[]; rows = GtkBox[]
    for txt in ["Upper FC for peak relative to mean peak for plate (default = 1.5)",
                "Lower FC for peak relative to mean peak for plate (default = 0.5)",
                "Upper FC of final biomass relative to mean peak of the plate (default = 0.1)"]
        row=GtkBox(:h,6); lbl=GtkLabel(txt); ent=GtkEntry(); ent.hexpand=true
        push!(row,lbl,ent); row.visible=false; push!(v1,row);
        push!(rows,row); push!(entries,ent)
    end
    signal_connect(screen_cb,"toggled") do cb for r in rows; r.visible=cb.active; end end
    btn = GtkButton("Next"); push!(v1,btn)
    signal_connect(btn,"clicked") do _
        if screen_cb.active
            vals=Float64[]; for ent in entries
                txt=strip(ent.text)
                if isempty(txt)
                    Gtk4.warn_dialog(() -> nothing,
                      "Fill all three fields first.",win1); return
                end
                try push!(vals, parse(Float64, txt))
                catch
                    Gtk4.warn_dialog(() -> nothing,
                      "\"$txt\" is not a valid number.",win1); return
                end
            end; vals_ref[] = vals
        end
        destroy(win1); notify(c1)
    end
    show(win1); wait(c1)
    if screen_cb.active
        action, tt, pk = show_plot_options("Screen")
        if action==:back; return run_plot_gui(options_path); end
        open(options_path,"w") do io
            JSON.print(io, Dict(
                "screen" => "True",
                "upper_FC" => vals_ref[][1],
                "lower_FC" => vals_ref[][2],
                "upper_final" => vals_ref[][3],
                "plot_time_traces" => string(tt),
                "plot_peaks" => string(pk)
            ),4)
        end
        return
    end

    a2 = show_stage2()
    if a2==:back; return run_plot_gui(options_path); end

    if a2==:yes
        a3, tt, pk = show_plot_options("Full Plate")
        if a3==:back; return run_plot_gui(options_path); end
        open(options_path,"w") do io
            JSON.print(io, Dict(
                "full_plate"=>"True",
                "plot_time_traces"=>string(tt),
                "plot_peaks"=>string(pk)
            ),4)
        end
        return
    else
        a4 = show_stage3b()
        if a4==:back; return run_plot_gui(options_path)
        elseif a4==:conditions
            act, tt, pk = show_plot_options("Specific Conditions")
            if act==:back; return run_plot_gui(options_path); end
            chosen = show_condition_selector()
            if chosen===nothing; return run_plot_gui(options_path); end
            open(options_path,"w") do io
                JSON.print(io, Dict(
                    "plot_conditions"=>"True",
                    "plot_time_traces"=>string(tt),
                    "plot_peaks"=>string(pk),
                    "conditions"=>chosen
                ),4)
            end
            return
        else
			# Final stage: plot specific wells
			local act, tt, pk = show_plot_options("Wells")
			if act == :back
				return run_plot_gui(options_path)
			end

			# Use show_stage4() which now reads plate_geometry from JSON
			local wells = show_stage4()
			open(options_path, "w") do io
				JSON.print(io, Dict(
					"plot_wells"       => "True",
					"plot_time_traces" => string(tt),
					"plot_peaks"       => string(pk),
					"wells"            => wells
				), 4)
			end
            return
        end
    end
end
