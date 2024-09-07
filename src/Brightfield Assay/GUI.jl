using Gtk4, JSON

win = GtkWindow("Experiment Configuration", 400, 200)
vbox = GtkBox(:v)
push!(win, vbox)
directories = []
select_dirs_button = GtkButton("Select Folders")
push!(vbox, select_dirs_button)

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
            push!(directories, sel...)
            destroy(dlg)
        end
    end
    show(dlg)
end

signal_connect(select_dirs_button, "clicked") do widget
    select_directories(widget)
end

dust_correction_checkbox = GtkCheckButton("Dust Correction")
batch_processing_checkbox = GtkCheckButton("Batch Processing")
push!(vbox, dust_correction_checkbox)
push!(vbox, batch_processing_checkbox)

done_button = GtkButton("Done")
push!(vbox, done_button)

function on_done(button)
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
    config = Dict(
        "directories" => directories,
        "dust_correction" => dust_correction,
        "batch_processing" => batch_processing 
    )
    
    open("experiment_config.json", "w") do f
        write(f, JSON.json(config, 4))
    end
    
    close(win)
end

signal_connect(done_button, "clicked") do widget
    on_done(widget)
end

show(win)
