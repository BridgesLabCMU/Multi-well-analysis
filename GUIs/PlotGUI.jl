using Gtk4
using Gtk4.GLib: @idle_add
using JSON

# — Helper: write plot_options.json —
function save_plot_options(dict::Dict)
    open("plot_options.json", "w") do io
        JSON.print(io, dict, 4)
    end
end

# — Stage 1: Screen GUI —
function run_screen_gui()
    win  = GtkWindow("Screen Configuration")
    vbox = GtkBox(:v, 10)
    push!(win, vbox)

    screen_cb = GtkCheckButton("Screen")
    push!(vbox, screen_cb)

    labels = [
      "Upper FC for peak relative to mean peak for plate (default = 1.5)",
      "Lower FC for peak relative to mean peak for plate (default = 0.5)",
      "Upper FC of final biomass relative to mean peak of the plate (default = 0.1)"
    ]
    rows, entries = GtkBox[], GtkEntry[]
    for txt in labels
        row = GtkBox(:h, 6)
        lbl = GtkLabel(txt)
        ent = GtkEntry(); ent.hexpand = true
        push!(row, lbl, ent)
        row.visible = false
        push!(vbox, row)
        push!(rows, row)
        push!(entries, ent)
    end

    signal_connect(screen_cb, "toggled") do cb
        for row in rows
            row.visible = cb.active
        end
    end

    next_btn = GtkButton("Next")
    push!(vbox, next_btn)
    signal_connect(next_btn, "clicked") do _
        if screen_cb.active
            vals = Float64[]
            for ent in entries
                txt = strip(ent.text)
                if isempty(txt)
                    Gtk4.warn_dialog(() -> nothing,
                      "Please fill in all three numeric fields before proceeding.", win)
                    return
                end
                try
                    push!(vals, parse(Float64, txt))
                catch
                    Gtk4.warn_dialog(() -> nothing,
                      "\"$txt\" is not a valid number.", win)
                    return
                end
            end
            save_plot_options(Dict(
                "screen"      => "True",
                "upper_FC"    => vals[1],
                "lower_FC"    => vals[2],
                "upper_final" => vals[3]
            ))
            destroy(win)
        else
            destroy(win)
            run_plot_full_plate_gui()
        end
    end

    show(win)
end

# — Stage 2: Full‐Plate vs. Specific Options —
function run_plot_full_plate_gui()
    win  = GtkWindow("Plot Full Plate")
    vbox = GtkBox(:v, 10)
    push!(win, vbox)

    plot_cb = GtkCheckButton("Plot full plate")
    push!(vbox, plot_cb)

    btns   = GtkBox(:h, 6)
    back   = GtkButton("Back")
    nxt    = GtkButton("Next")
    push!(btns, back, nxt)
    push!(vbox, btns)

    signal_connect(back, "clicked") do _
        destroy(win)
        run_screen_gui()
    end

    signal_connect(nxt, "clicked") do _
        if plot_cb.active
            destroy(win)
            run_full_plate_options_gui()
        else
            destroy(win)
            run_plot_specific_options_gui()
        end
    end

    show(win)
end

# — Stage 3a: Full‐Plate Options —
function run_full_plate_options_gui()
    win  = GtkWindow("Full Plate Options")
    vbox = GtkBox(:v, 10)
    push!(win, vbox)

    tt_cb   = GtkCheckButton("Plot time traces")
    peak_cb = GtkCheckButton("Plot peaks")
    push!(vbox, tt_cb, peak_cb)

    btns   = GtkBox(:h, 6)
    back   = GtkButton("Back")
    nxt    = GtkButton("Next")
    push!(btns, back, nxt)
    push!(vbox, btns)

    signal_connect(back, "clicked") do _
        destroy(win)
        run_plot_full_plate_gui()
    end

    signal_connect(nxt, "clicked") do _
        if !(tt_cb.active || peak_cb.active)
            Gtk4.warn_dialog(() -> nothing,
              "Please select at least one option before proceeding.", win)
            return
        end
        save_plot_options(Dict(
            "full_plate"         => "True",
            "plot_time_traces"   => tt_cb.active ? "True" : "False", 
            "plot_peaks"         => peak_cb.active ? "True" : "False" 
        ))
        destroy(win)
    end

    show(win)
end

# — Stage 3b: Specific Conditions or Wells —
function run_plot_specific_options_gui()
    win  = GtkWindow("Plot Specific Options")
    vbox = GtkBox(:v, 10)
    push!(win, vbox)

    spec_cb  = GtkCheckButton("Plot specific conditions")
    wells_cb = GtkCheckButton("Plot specific wells")
    push!(vbox, spec_cb, wells_cb)

    btns = GtkBox(:h, 6)
    back = GtkButton("Back")
    nxt  = GtkButton("Next")
    push!(btns, back, nxt)
    push!(vbox, btns)

    signal_connect(back, "clicked") do _
        destroy(win)
        run_plot_full_plate_gui()
    end

    signal_connect(nxt, "clicked") do _
        if wells_cb.active && !spec_cb.active
            destroy(win)
            run_plate_selector_no_labels()
        elseif spec_cb.active && !wells_cb.active
            if !isfile("experiment_config.json")
                Gtk4.warn_dialog(() -> nothing,
                  "Cannot find experiment_config.json.", win)
                return
            end
            cfg = try JSON.parsefile("experiment_config.json") catch e
                Gtk4.warn_dialog(() -> nothing,
                  "Error parsing JSON: $e", win)
                return
            end
            if !haskey(cfg, "conditions") || !(cfg["conditions"] isa Dict)
                Gtk4.warn_dialog(() -> nothing,
                  "No \"conditions\" field found in experiment_config.json.", win)
                return
            end
            destroy(win)
            run_conditions_selector_gui()
        else
            Gtk4.warn_dialog(() -> nothing,
              "Select exactly one option.", win)
        end
    end

    show(win)
end

# — Stage 4: Conditions Selector —
function run_conditions_selector_gui()
    win     = GtkWindow("Select Conditions")
    vbox    = GtkBox(:v, 10); push!(win, vbox)

    # Scrolled area + vertical box for the checkboxes
    scroll  = GtkScrolledWindow(); scroll.vexpand = true
    cond_box = GtkBox(:v, 6)
    scroll.child = cond_box
    push!(vbox, scroll)

    # Load your conditions dict
    cfg    = JSON.parsefile("experiment_config.json")
    conds  = cfg["conditions"]  # assume Dict{String,Any}

    # Store the checkbuttons in a vector so we can read them later
    checkboxes = GtkCheckButton[]
    for key in keys(conds)
        cb = GtkCheckButton(key)
        push!(cond_box, cb)
        push!(checkboxes, cb)
    end

    # Back / Next buttons
    btns    = GtkBox(:h, 6)
    back_bt = GtkButton("Back")
    next_bt = GtkButton("Next")
    push!(btns, back_bt, next_bt)
    push!(vbox, btns)

    signal_connect(back_bt, "clicked") do _
        destroy(win)
        run_plot_specific_options_gui()
    end

    signal_connect(next_bt, "clicked") do _
        # Simply iterate our stored checkboxes
        selected = [ cb.label for cb in checkboxes if cb.active ]

        if isempty(selected)
            Gtk4.warn_dialog(() -> nothing,
                "No conditions selected. Please pick at least one.", win)
            return
        end

        save_plot_options(Dict(
            "plot_conditions" => "True",
            "conditions"      => selected
        ))
        destroy(win)
    end

    show(win)
end

# — Stage 5: Specific Wells Selector —
function run_plate_selector_no_labels()
    win  = GtkWindow("Select Specific Wells")
    vbox = GtkBox(:v, 10)
    push!(win, vbox)

    plate_lbl = GtkLabel("Plate 1 of 10")
    push!(vbox, plate_lbl)

    rows, cols = 8, 12
    grid = GtkGrid(column_homogeneous=true, row_homogeneous=true,
                   column_spacing=2,     row_spacing=2)
    buttons = Matrix{GtkToggleButton}(undef, rows, cols)
    for r in 1:rows, c in 1:cols
        lbl = string(Char('A'+r-1), c)
        btn = GtkToggleButton(lbl)
        btn.width_request  = 40
        btn.height_request = 40
        buttons[r, c] = btn
        grid[c, r]    = btn
    end
    push!(vbox, grid)

    btns = GtkBox(:h, 6)
    back = GtkButton("Back")
    save = GtkButton("Save Plate")
    done = GtkButton("Done")
    push!(btns, back, save, done)
    push!(vbox, btns)

    current = Ref(Vector{Vector{String}}())

    signal_connect(back, "clicked") do _
        destroy(win)
        run_plot_specific_options_gui()
    end

    signal_connect(save, "clicked") do _
        sel = [ string(Char('A'+r-1), c)
                for r in 1:rows, c in 1:cols if buttons[r, c].active ]
        if length(current[]) < 10
            push!(current[], sel)
            @idle_add begin
                plate_lbl.label = "Plate $(length(current[])+1) of 10"
                for b in buttons; b.active = false; end
                false
            end
        else
            Gtk4.warn_dialog(() -> nothing,
              "Cannot save more than 10 plates.", win)
        end
    end

    signal_connect(done, "clicked") do _
        save_plot_options(Dict(
            "plot_wells" => "True",
            "wells"      => current[]
        ))
        destroy(win)
    end

    show(win)
end

# — Launch the flow —
run_screen_gui()
