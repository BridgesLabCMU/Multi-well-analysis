function run_conditions_gui(config_path::AbstractString = "experiment_config.json")
    # Load existing config
    cfg = Dict{String,Any}()
    if isfile(config_path)
        try
            cfg = JSON.parsefile(config_path)
        catch e
            @warn "Failed to parse $config_path, starting fresh: $e"
        end
    end

    # Determine plate geometry
    geom = get(cfg, "plate_geometry", "96-well")
    plate_map = Dict(
        "6-well"   => (2, 3),
        "12-well"  => (3, 4),
        "24-well"  => (4, 6),
        "48-well"  => (6, 8),
        "96-well"  => (8, 12),
        "384-well" => (16, 24)
    )
    if haskey(plate_map, geom)
        rows, cols = plate_map[geom]
    else
        error("Unsupported plate geometry: $geom. Supported: $(collect(keys(plate_map)))")
    end

    # Initialize data structures
    data_dict      = Dict{String, Vector{Vector{String}}}()
    all_labels     = Ref(Vector{String}())
    current_plates = Ref(Vector{Vector{String}}())

    # Load previous conditions if present
    if haskey(cfg, "conditions") && cfg["conditions"] isa Dict
        for (k, v) in cfg["conditions"]
            data_dict[k] = deepcopy(v)
            push!(all_labels[], k)
        end
    end

    # Build UI
    win    = GtkWindow("$geom Plate Selector")
    vbox   = GtkBox(:v, 5)
    push!(win, vbox)

    entry = GtkEntry()
    entry.placeholder_text = "Enter label here"
    push!(vbox, entry)

    plate_label = GtkLabel("Plate 1 of 10")  # adjust max plates if needed
    push!(vbox, plate_label)

    grid    = GtkGrid(column_homogeneous=true, row_homogeneous=true,
                      column_spacing=2, row_spacing=2)
    buttons = Matrix{GtkToggleButton}(undef, rows, cols)
    for r in 1:rows, c in 1:cols
        lab = string(Char('A' + r - 1), c)
        btn = GtkToggleButton(lab)
        btn.width_request  = 40
        btn.height_request = 40
        buttons[r, c] = btn
        grid[c, r] = btn
    end

    # Disable already used wells
    function update_grid!()
        plate_idx = length(current_plates[]) + 1
        used = Set{String}()
        for plates in values(data_dict)
            if length(plates) >= plate_idx
                for w in plates[plate_idx]
                    push!(used, w)
                end
            end
        end
        for r in 1:rows, c in 1:cols
            let btn = buttons[r,c], well = string(Char('A'+r-1), c)
                btn.sensitive = !(well in used)
            end
        end
    end

    push!(vbox, grid)
    update_grid!()

    # Buttons
    save_plate_btn = GtkButton("Save Plate")
    add_btn        = GtkButton("Add Selection")
    done_btn       = GtkButton("Done")
    btn_box        = GtkBox(:h, 5)
    push!(btn_box, save_plate_btn, add_btn, done_btn)
    push!(vbox, btn_box)

    label_dropdown = GtkDropDown(all_labels[])
    delete_btn     = GtkButton("Delete Label")
    del_box        = GtkBox(:h, 5)
    push!(del_box, label_dropdown, delete_btn)
    push!(vbox, del_box)

    # Enable/disable controls
    function update_sensitivities!()
        has_label = !isempty(strip(entry.text))
        save_plate_btn.sensitive = has_label
        add_btn.sensitive        = length(current_plates[]) >= 1
        done_btn.sensitive       = !isempty(all_labels[])
        delete_btn.sensitive     = !isempty(all_labels[])
    end

    signal_connect(entry, "changed")   do _ update_sensitivities!() end
    for b in buttons
        signal_connect(b, "toggled")  do _ update_sensitivities!() end
    end
    update_sensitivities!()
    signal_connect(save_plate_btn, "clicked") do _
        key = strip(entry.text)
        if isempty(key)
            Gtk4.warn_dialog(() -> nothing,
                "Please enter a non-empty label before saving a plate.", win)
            return
        end
        if haskey(data_dict, key)
            Gtk4.warn_dialog(() -> nothing,
                "Label $key already exists. Choose a different name.", win)
            return
        end

        sel = String[]
        for r in 1:rows, c in 1:cols
            if buttons[r, c].active
                push!(sel, string(Char('A'+r-1), c))
            end
        end

        if length(current_plates[]) < 10
            push!(current_plates[], sel)
            @idle_add begin
                plate_label.label = "Plate $(length(current_plates[])+1) of 10"
                for btn in buttons; btn.active = false; end
				update_grid!()
                update_sensitivities!()
                false
            end
        else
            Gtk4.warn_dialog(() -> nothing,
                "Already saved 10 plates—no more can be added.", win)
        end
    end

    signal_connect(add_btn, "clicked") do _
        if any(b.active for b in buttons)
            Gtk4.warn_dialog(() -> nothing,
                "You have selected wells but have not saved this plate.", win)
            return
        end

        key = strip(entry.text)
        if isempty(key)
            Gtk4.warn_dialog(() -> nothing,
                "Please enter a non-empty label before finalizing.", win)
            return
        end

        data_dict[key] = copy(current_plates[])
        push!(all_labels[], key)

        @idle_add begin
            current_plates[]  = Vector{Vector{String}}()
            plate_label.label = "Plate 1 of 10"
            entry.text        = ""
            for btn in buttons; btn.active = false; end
            mdl = Gtk4.model(label_dropdown)
            push!(mdl, key)
            update_grid!()
            update_sensitivities!()
            false
        end
    end

    signal_connect(delete_btn, "clicked") do _
        idx = Gtk4.get_gtk_property(label_dropdown, :selected, Int)
        if idx < 0
            Gtk4.warn_dialog(() -> nothing,
                "No label selected to delete.", win)
            return
        end

        key = all_labels[][idx+1]
        pop!(data_dict, key)
        splice!(all_labels[], idx+1)

        @idle_add begin
            mdl = Gtk4.model(label_dropdown)
            deleteat!(mdl, idx+1)
            update_grid!()
            update_sensitivities!()
            false
        end
    end

    c = Condition()
    signal_connect(done_btn, "clicked") do _
        cfg2 = Dict{String,Any}()
        if isfile(config_path)
            try
                cfg2 = JSON.parsefile(config_path)
            catch
            end
        end
        cfg2["conditions"] = data_dict
        open(config_path, "w") do io
            JSON.print(io, cfg2, 4)
        end

        destroy(win)
        notify(c)
    end

    show(win)
    @async Gtk4.GLib.glib_main()
    wait(c)

    return nothing
end
