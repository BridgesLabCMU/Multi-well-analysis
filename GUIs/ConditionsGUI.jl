using Gtk4
using Gtk4.GLib: @idle_add
using JSON

const rows = 8
const cols = 12

# — Data structures —
const data_dict      = Dict{String, Vector{Vector{String}}}()  # final: label ⇒ plates
const all_labels     = Ref(Vector{String}())                   # label order
const current_plates = Ref(Vector{Vector{String}}())           # buffer up to 10 plates

# — Build GUI —
win = GtkWindow("96-Well Plate Selector")
vbox = GtkBox(:v, 5)

# 1) Label entry
entry = GtkEntry(); entry.placeholder_text = "Enter label here"
push!(vbox, entry)

# 2) Status label
plate_label = GtkLabel("Plate 1 of 10")
push!(vbox, plate_label)

# 3) 96-well grid
grid = GtkGrid(column_homogeneous=true, row_homogeneous=true,
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
push!(vbox, grid)

# 4) Action buttons
btn_box        = GtkBox(:h, 5)
save_plate_btn = GtkButton("Save Plate")
add_btn        = GtkButton("Add Selection")
done_btn       = GtkButton("Done with GUI")
push!(btn_box, save_plate_btn, add_btn, done_btn)
push!(vbox, btn_box)

# 5) Delete-label UI with GtkDropDown
label_dropdown = GtkDropDown(all_labels[])   # model from all_labels[] :contentReference[oaicite:2]{index=2}
delete_btn     = GtkButton("Delete Label")
del_box        = GtkBox(:h, 5)
push!(del_box, label_dropdown, delete_btn)
push!(vbox, del_box)

# — Helpers — 

function update_sensitivities!()
    has_label = !isempty(strip(entry.text))
    save_plate_btn.sensitive = has_label
    add_btn.sensitive        = length(current_plates[]) ≥ 1
    done_btn.sensitive       = length(data_dict) ≥ 1
    delete_btn.sensitive     = !isempty(all_labels[])
end

signal_connect(entry, "changed") do _ update_sensitivities!() end
for b in buttons
    signal_connect(b, "toggled") do _ update_sensitivities!() end
end

# — Callbacks — 

# SAVE PLATE: only disallow duplicate final labels
signal_connect(save_plate_btn, "clicked") do _
    key = strip(entry.text)
    if haskey(data_dict, key)
        Gtk4.warn_dialog(() -> nothing,
          "Label $key already exists. Choose a different name before saving.",
          win)
        return
    end

    # capture zero-or-more wells
    sel = [ string(Char('A'+r-1), c)
            for r in 1:rows, c in 1:cols if buttons[r, c].active ]

    if length(current_plates[]) < 10
        push!(current_plates[], sel)
        @idle_add begin
            plate_label.label = "Plate $(length(current_plates[])+1) of 10"
            for btn in buttons; btn.active = false; end
            update_sensitivities!()
            false
        end
    else
        Gtk4.warn_dialog(() -> nothing,
          "You have already saved 10 plates—no more can be added.",
          win)
    end
end

# ADD SELECTION: finalize and add to dropdown model
signal_connect(add_btn, "clicked") do _
    if any(b.active for b in buttons)
        Gtk4.warn_dialog(() -> nothing,
          "You have selected wells but have not saved this plate.",
          win)
        return
    end

    key = strip(entry.text)
    if isempty(key)
        Gtk4.warn_dialog(() -> nothing,
          "Please enter a non-empty label before finalizing.",
          win)
        return
    end

    # record data
    data_dict[key] = copy(current_plates[])
    push!(all_labels[], key)

    @idle_add begin
        # clear UI
        current_plates[]  = Vector{Vector{String}}()
        plate_label.label = "Plate 1 of 10"
        entry.text        = ""
        for btn in buttons; btn.active = false; end

        # update drop-down model in place
        mdl = Gtk4.model(label_dropdown)
        push!(mdl, key)

        update_sensitivities!()
        false
    end
end

# DELETE LABEL: remove from data and dropdown
signal_connect(delete_btn, "clicked") do _
    idx = Gtk4.get_gtk_property(label_dropdown, :selected, Int)
    if idx < 0
        Gtk4.warn_dialog(() -> nothing,
          "No label selected to delete.",
          win)
        return
    end

    key = all_labels[][idx+1]
    pop!(data_dict, key)
    splice!(all_labels[], idx+1)

    @idle_add begin
        # remove from dropdown model
        mdl = Gtk4.model(label_dropdown)
        deleteat!(mdl, idx+1)
        update_sensitivities!()
        false
    end
end

# DONE WITH GUI
signal_connect(done_btn, "clicked") do _
    # 1) Attempt to load the existing config, or start fresh
    config = Dict{String,Any}()
    if isfile("experiment_config.json")
        try
            config = JSON.parsefile("experiment_config.json")
        catch e
            @warn "Failed to parse experiment_config.json, starting fresh: $e"
        end
    end

    # 2) Insert (or overwrite) the data_dict under a top-level key
    config["conditions"] = data_dict

    # 3) Write it back out (overwriting)
    open("experiment_config.json", "w") do io
        JSON.print(io, config, 4)
    end

    # 4) Finally close the GUI
    destroy(win)
end

# — Initialize & show —
update_sensitivities!()
push!(win, vbox)
show(win)
