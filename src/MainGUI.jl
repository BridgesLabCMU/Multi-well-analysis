# TODO: 1. fix the dropdown to enable alternative media selection

function display_images!(stack, masks, overlay)
    normalized = similar(stack)
	fpMax = maximum(stack)
	fpMin = minimum(stack)
	fpMean = (fpMax - fpMin) / 2.0 + fpMin
	normalized = normalize_local_contrast_output(normalized, stack, copy(stack), 101, fpMean)
	normalized = Gray{N0f8}.(normalized)
    @inbounds for i in CartesianIndices(normalized)
        gray_val = RGB{N0f8}(normalized[i], normalized[i], normalized[i])
        overlay[i] = masks[i] ? RGB{N0f8}(0,1,1) : gray_val
    end
    imshow(overlay)
end

function preprocess_noreg!(img_stack, normalized_stack, blockDiameter, sig)       
    @inbounds for t in 1:size(img_stack, 3)
        img = img_stack[:,:,t]
        img_copy = img_stack[:,:,t] 
        img_normalized = normalize_local_contrast(img, img_copy, 
                                    blockDiameter)
        normalized_stack[:,:,t] = imfilter(img_normalized, Kernel.gaussian(sig))
    end
end

function timelapse_test!(test_images, fixed_thresh, blockDiameter, sig)
    trial_image = load(test_images[1]; lazyio=true)
    height, width = size(trial_image)
    ntimepoints = length(test_images)
    stack = Array{eltype(trial_image)}(undef, height, width, ntimepoints)
    read_images!(ntimepoints, stack, test_images)
    stack = Float64.(stack)
    normalized_stack = similar(stack)
    preprocess_noreg!(stack, normalized_stack, blockDiameter, sig)
    masks = zeros(Bool, size(stack))
    compute_mask!(normalized_stack, masks, fixed_thresh, ntimepoints)
    overlay = zeros(RGB{N0f8}, size(stack)...)
    display_images!(stack, masks, overlay)
    return nothing
end

function image_test!(image_path, fixed_thresh, blockDiameter, sig)
    image = load(image_path)
    image_copy = copy(image)
    image_normalized = normalize_local_contrast(image, image_copy, blockDiameter)
    image_normalized = imfilter(image_normalized, Kernel.gaussian(sig))
    mask = image_normalized .> fixed_thresh
    overlay = zeros(RGB{N0f8}, size(image)...)
    display_images!(Float64.(image), mask, overlay)
    return nothing
end

function with_help(w::GtkWidget, help_text::AbstractString; spacing=4)
    help_lbl = GtkLabel("❓")
    set_gtk_property!(help_lbl, :has_tooltip,  true)
    set_gtk_property!(help_lbl, :tooltip_text, help_text)
    set_gtk_property!(help_lbl, :margin_start, 4)
    row = GtkBox(:h, spacing)
    set_gtk_property!(row, :hexpand, true)
    set_gtk_property!(w, :hexpand, true)
    push!(row, w)
    push!(row, help_lbl)

    return row
end

function GUI_main()
    win = GtkWindow("Experiment Configuration", 600, 300)
	scrolled = GtkScrolledWindow()
	vbox      = GtkBox(:v)
	scrolled[] = vbox
	win[]      = scrolled
    
    directories = []
    test_images = []
    bulk_data = ""
    Imin = ""
    Imax = ""
    fixed_thresh = 0.03 

	notes_buffer = GtkTextBuffer()
	notes_view   = GtkTextView(buffer = notes_buffer)
	notes_view.wrap_mode = 3 
	sw_notes = GtkScrolledWindow()
	sw_notes.hexpand = true
	sw_notes.vexpand = false
	push!(vbox, GtkLabel("Experiment notes:"))
	push!(vbox, with_help(sw_notes, "Write notes to describe what you did in the experiment. For your own use."))
	sw_notes.child = notes_view

	media_options_path = "media_options.txt"
	media_lines = readlines(media_options_path)
	media_dropdown = GtkComboBoxText()
	for opt in media_lines
		push!(media_dropdown, opt)
	end
	media_dropdown.active = 0 
	push!(vbox, GtkLabel("Media:"))
	push!(vbox, with_help(media_dropdown, "Select the media you used in this experiment."))
    
	plate_geometries = ["6-well", "12-well", "24-well", "48-well", "96-well", "384-well"]
    plate_geom_box = GtkBox(:v, 4)
    push!(vbox, GtkLabel("Plate geometry (select exactly one):"))
    for geom in plate_geometries
        cb = GtkCheckButton(geom)
        push!(plate_geom_box, cb)
    end
    push!(vbox, plate_geom_box)

	select_dirs_button = GtkButton("Select experiment folders")
    push!(vbox, with_help(select_dirs_button, "Select the folders that contain images. Be sure the images are not contained in a sub-folder. Multiple choice is possible, and you can select one, click OK, select another, etc."))
	dir_dropdown   = GtkDropDown(directories)
	remove_dir_btn = GtkButton("Remove folder")
	remove_dir_btn.sensitive = false 
    clear_all_btn = GtkButton("Clear all")
    clear_all_btn.sensitive = false
	push!(vbox, GtkLabel("Selected experiment folders:"))
	hbox_dirs = GtkBox(:h, 4)
	push!(hbox_dirs, dir_dropdown, remove_dir_btn, clear_all_btn)
	push!(vbox, hbox_dirs) 

    select_bulk_data_button = GtkButton("Select bulk data file")
    push!(vbox, with_help(select_bulk_data_button, "Select the file containing bulk data and metadata for the experiment. This can be generated in Gen5 post-experiment. Ask around if you don't know how to generate it."))
    
    select_Imin_button = GtkButton("Select Imin file")
    push!(vbox, with_help(select_Imin_button, "Optional. Select .tif file containing a shutter-closed image. If selected, the analysis pipeline will convert the pixel values to 'optical density' values."))
    
    select_Imax_button = GtkButton("Select Imax file (for single image analysis)")
    push!(vbox, with_help(select_Imax_button, "Optional. Select .tif file containing a media only (blank) image. Necessary, in conjunction with Imin, if you want to convert the pixel values to 'optical density' values. For timelapses, not strictly necessary if the first timepoint is effectively a blank."))
    
	block_adj   = GtkAdjustment(101, 101, 501, 2, 20, 0) 
	block_spin  = GtkSpinButton(block_adj, 1, 0)
	push!(vbox, GtkLabel("Block diameter:"))
	push!(vbox, with_help(block_spin, "Select the block diameter,which will determine the size of the local averaging operation. Use 101 pixels if you don't know what to set. Otherwise, the maximum linear dimension of an object you want to consider 'biofilm' is a good bet."))

    adj = GtkAdjustment(0.04, 0., 1., 0.00001, 1.0, 0.0) 
    spin_button = GtkSpinButton(adj, 0.00001, 4)
	Gtk4.value(spin_button, 0.04)
    push!(vbox, GtkLabel("Threshold:"))
	push!(vbox, with_help(spin_button,
		"This controls the fixed threshold above which pixels are considered “signal”. You can fine‑tune it between 0.0 and 1.0, either by manual typing or with the +/-. Smaller values will lead to more pixels being considered biofilm. 0.04 is a reasonable starting value."))

    
    select_test_images_button = GtkButton("Select one tif or a tif series to test the threshold")
    push!(vbox, with_help(select_test_images_button, "Optional. Select a tif or selection of tifs (timelapse) or a single tif series to test the masking results for the chosen threshold. This selection is not remembered, so if you want to re-test with a different threshold, you need to select the images again."))

	signal_connect(media_dropdown, "changed") do widget, others...
	  idx = media_dropdown.active
	  str = Gtk4.active_text(media_dropdown)
	  println("Active element is \"$str\" at index $idx")
	end

    function select_directories(button)
        dlg = GtkFileChooserDialog(
            "Select Folders",
            win,
            Gtk4.FileChooserAction_SELECT_FOLDER,
                                (("_Cancel", Gtk4.ResponseType_CANCEL),
                                 ("_Open", Gtk4.ResponseType_ACCEPT)))
        
        dlgp = GtkFileChooser(dlg)
        Gtk4.G_.set_select_multiple(dlgp, true)
        signal_connect(dlg, "response") do widget, response_id
            if response_id == Gtk4.ResponseType_ACCEPT
                filename_list = Gtk4.G_.get_files(dlgp)
                sel = String[Gtk4.GLib.G_.get_path(Gtk4.GFile(f)) for f in Gtk4.GListModel(filename_list)]
				mdl = Gtk4.model(dir_dropdown)
				for path in sel
					push!(directories, path)
					push!(mdl, path)
				end
                remove_dir_btn.sensitive = true
				clear_all_btn.sensitive = true
            end
			destroy(dlg)
        end
        show(dlg)
    end

    signal_connect(select_dirs_button, "clicked") do widget
        select_directories(widget)
    end

	signal_connect(remove_dir_btn, "clicked") do _
		idx = Gtk4.get_gtk_property(dir_dropdown, :selected, Int)
		if idx < 0
			Gtk4.warn_dialog(() -> nothing,
				"No folder selected to remove.", win)
		else
			# delete from your directories array…
			deleteat!(directories, idx+1)
			# …and from the dropdown model
			mdl = Gtk4.model(dir_dropdown)
			deleteat!(mdl, idx+1)
			clear_all_btn.sensitive   = !isempty(directories)
            remove_dir_btn.sensitive = !isempty(directories)
		end
	end

	signal_connect(clear_all_btn, "clicked") do _
		empty!(directories)
		mdl = Gtk4.model(dir_dropdown)
		empty!(mdl)

		remove_dir_btn.sensitive = false
		clear_all_btn.sensitive  = false
	end
    
    function select_bulk_data(button)
        dlg = GtkFileChooserDialog(
            "Select bulk data file",
            win,
            Gtk4.FileChooserAction_OPEN,
                                (("_Cancel", Gtk4.ResponseType_CANCEL),
                                 ("_Open", Gtk4.ResponseType_ACCEPT)))
        
        dlgp = GtkFileChooser(dlg)
        Gtk4.G_.set_select_multiple(dlgp, false)
        signal_connect(dlg, "response") do widget, response_id
            if response_id == Gtk4.ResponseType_ACCEPT
                filename = Gtk4.G_.get_file(dlgp)
                sel = Gtk4.GLib.G_.get_path(Gtk4.GFile(filename))
                bulk_data = sel
            end
			destroy(dlg)
        end
        show(dlg)
    end

    signal_connect(select_bulk_data_button, "clicked") do widget
        select_bulk_data(widget)
    end

    function select_Imin(button)
        dlg = GtkFileChooserDialog(
            "Select Imin file",
            win,
            Gtk4.FileChooserAction_OPEN,
                                (("_Cancel", Gtk4.ResponseType_CANCEL),
                                 ("_Open", Gtk4.ResponseType_ACCEPT)))
        
        dlgp = GtkFileChooser(dlg)
        Gtk4.G_.set_select_multiple(dlgp, false)
        signal_connect(dlg, "response") do widget, response_id
            if response_id == Gtk4.ResponseType_ACCEPT
                filename = Gtk4.G_.get_file(dlgp)
                sel = Gtk4.GLib.G_.get_path(Gtk4.GFile(filename))
                Imin = sel
            end
			destroy(dlg)
        end
        show(dlg)
    end

    signal_connect(select_Imin_button, "clicked") do widget
        select_Imin(widget)
    end

    function select_Imax(button)
        dlg = GtkFileChooserDialog(
            "Select Imax file",
            win,
            Gtk4.FileChooserAction_OPEN,
                                (("_Cancel", Gtk4.ResponseType_CANCEL),
                                 ("_Open", Gtk4.ResponseType_ACCEPT)))
        
        dlgp = GtkFileChooser(dlg)
        Gtk4.G_.set_select_multiple(dlgp, false)
        signal_connect(dlg, "response") do widget, response_id
            if response_id == Gtk4.ResponseType_ACCEPT
                filename = Gtk4.G_.get_file(dlgp)
                sel = Gtk4.GLib.G_.get_path(Gtk4.GFile(filename))
                Imax = sel
            end
			destroy(dlg)
        end
        show(dlg)
    end

    signal_connect(select_Imax_button, "clicked") do widget
        select_Imax(widget)
    end

    function select_test_images(button)
        dlg = GtkFileChooserDialog(
            "Select test image(s)",
            win,
            Gtk4.FileChooserAction_OPEN,
                                (("_Cancel", Gtk4.ResponseType_CANCEL),
                                 ("_Open", Gtk4.ResponseType_ACCEPT)))
        
        dlgp = GtkFileChooser(dlg)
        Gtk4.G_.set_select_multiple(dlgp, true)
        signal_connect(dlg, "response") do widget, response_id
            if response_id == Gtk4.ResponseType_ACCEPT
                filename_list = Gtk4.G_.get_files(dlgp)
                sel = String[Gtk4.GLib.G_.get_path(Gtk4.GFile(f)) for f in Gtk4.GListModel(filename_list)]
                push!(test_images, sel...)
            end
			destroy(dlg)
        end
        show(dlg)
    end
    
    signal_connect(select_test_images_button, "clicked") do widget
        select_test_images(widget)
    end
    
    test_button = GtkButton("Test")
    push!(vbox, with_help(test_button, "Generate a mask of the selected test images with the chosen threshold. Might need to wait a bit for large tif series."))

    function on_test(button)
        fixed_thresh = Gtk4.value(spin_button)
        ntimepoints = length(test_images)
        if ntimepoints == 1
            image_test!(test_images[1], fixed_thresh, 101, 2)
        else
            timelapse_test!(sort(test_images, lt=natural), fixed_thresh, 101, 2)
        end
        empty!(test_images)
    end
    
    signal_connect(test_button, "clicked") do widget
        on_test(widget)
    end

    dust_correction_checkbox = GtkCheckButton("Dust correction")
    batch_processing_checkbox = GtkCheckButton("Timelapse imaging")
    push!(vbox, with_help(dust_correction_checkbox, "Check if you performed timelapse imaging and dust is an issue in the first frame. This will cause the mask to ignore areas that are very dark in the first frame for the entire timelapse. Do not check if you are analyzing a single timepoint!"))
    push!(vbox, with_help(batch_processing_checkbox, "Check if performing timelapse imaging. Leave unchecked otherwise."))

    done_button = GtkButton("Done")
    push!(vbox, with_help(done_button, "Generate .json file with selected options."))

    function on_done(button)
        # Read value from spin button
        fixed_thresh = spin_button.value
        blockDiameter = block_spin.value

        # Validate exactly one geometry
		selected = [cb for cb in plate_geom_box if cb.active]
		if length([cb for cb in plate_geom_box if cb.active]) != 1
			Gtk4.warn_dialog(() -> nothing,
				"Please select exactly one plate geometry.", win)
			return false
		end
        chosen_geom = Gtk4.get_gtk_property(selected[1], :label, String)
        
        if dust_correction_checkbox.active == true
            dust_correction = "True" 
        else
            dust_correction = "False" 
        end
        if batch_processing_checkbox.active == true
            batch_processing = "True" 
        else
            batch_processing = "False" 
        end
		selected_media = Gtk4.active_text(media_dropdown)
		notes_text = get_gtk_property(notes_buffer, 
                                      :text, String)

        config = Dict(
			"media" => selected_media,
			"notes" => notes_text,
            "plate_geometry" => chosen_geom,
            "images_directory" => directories,
            "bulk_data" => bulk_data,
            "Imin_path" => Imin,
            "Imax_path" => Imax,
            "fixed_thresh" => fixed_thresh,
            "blockDiam" => blockDiameter,
            "dust_correction" => dust_correction,
            "batch_processing" => batch_processing 
        )
        
        open("experiment_config.json", "w") do f
            JSON.print(f, config)
        end
		return true
    end
    if !isinteractive()
        # non-interactive: block until Done
        c = Condition()
        signal_connect(done_button, "clicked") do widget
            if on_done(widget)
                destroy(win)
                notify(c)
            end
        end

        show(win)
        @async Gtk4.GLib.glib_main()
        wait(c)
    else
        # interactive: just close on Done
        signal_connect(done_button, "clicked") do widget
            if on_done(widget)
                destroy(win)
            end
        end
        show(win)
    end
end
