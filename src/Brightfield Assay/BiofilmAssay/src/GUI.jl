using Gtk4
using JSON
using NaturalSort: sort, natural
using TiffImages: load
using IntegralArrays: IntegralArray
using IntervalSets: width, leftendpoint, rightendpoint, Interval, ±
using Images: imfilter, mapwindow, adjust_histogram!, LinearStretching, Gray, N0f16, N0f8, Kernel, warp, axes, KernelFactors, RGB 
using StatsBase: quantile
using ImageView

function read_images!(files, arr, height, width, ntimepoints)
    @inbounds for t in 1:ntimepoints
        @views file = files[t]
        arr[:, :, t] = load(file)
    end
    return nothing 
end

function compute_mask!(stack, masks, fixed_thresh, ntimepoints)
    @inbounds for t in 1:ntimepoints
        @views masks[:,:,t] = stack[:,:,t] .> fixed_thresh
    end
end

function display_images!(stack, masks, overlay)
	flat_stack = vec(stack)
    img_min = quantile(flat_stack, 0.0035)
    img_max = quantile(flat_stack, 0.9965)
    adjust_histogram!(stack, LinearStretching(src_minval=img_min, src_maxval=img_max, dst_minval=0, dst_maxval=1))
	stack = Gray{N0f8}.(stack)
    @inbounds for i in CartesianIndices(stack)
        gray_val = RGB{N0f8}(stack[i], stack[i], stack[i])
        overlay[i] = masks[i] ? RGB{N0f8}(0,1,1) : gray_val
    end
    imshow(overlay)
end

function normalize_local_contrast(img, img_copy, blockDiameter)
	img = 1 .- img
	img_copy = 1 .- img_copy
	img_copy = imfilter(img_copy, Kernel.gaussian(blockDiameter))
	img = img - img_copy
    return img 
end

function stack_preprocess!(img_stack, normalized_stack, blockDiameter, sig)       
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
    read_images!(test_images, stack, height, width, ntimepoints)
    stack = Float64.(stack)
    normalized_stack = similar(stack)
    stack_preprocess!(stack, normalized_stack, blockDiameter, sig)
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
    display_images!(image, mask, overlay)
    return nothing
end

function main()
    win = GtkWindow("Experiment Configuration", 200, 300)
    vbox = GtkBox(:v)
    push!(win, vbox)
    
    directories = []
    test_images = []
    Imins = []
    Imaxes = []
    fixed_thresh = 0.03 

    select_dirs_button = GtkButton("Select experiment folders")
    push!(vbox, select_dirs_button)
    
    select_Imin_button = GtkButton("Select Imin file")
    push!(vbox, select_Imin_button)
    
    select_Imax_button = GtkButton("Select Imax file (for single image analysis)")
    push!(vbox, select_Imax_button)
    
    adj = GtkAdjustment(0.03, 0., 1., 0.001, 1.0, 0.0) 
    spin_button = GtkSpinButton(adj, 0.001, 3)
	Gtk4.value(spin_button, 0.03)
    push!(vbox, GtkLabel("Threshold:"))
    push!(vbox, spin_button)
    
    select_test_images_button = GtkButton("Select one tif or a tif series to test the threshold")
    push!(vbox, select_test_images_button)

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

    function select_Imin(button)
        dlg = GtkFileChooserDialog(
            "Select Imin file",
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
                push!(Imins, sel...)
                destroy(dlg)
            end
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
        Gtk4.G_.set_select_multiple(dlgp, true)
        signal_connect(dlg, "response") do widget, response_id
            if response_id == Gtk4.ResponseType_ACCEPT
                filename_list = Gtk4.G_.get_files(dlgp)
                sel = String[Gtk4.GLib.G_.get_path(Gtk4.GFile(f)) for f in Gtk4.GListModel(filename_list)]
                push!(Imaxes, sel...)
                destroy(dlg)
            end
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
                destroy(dlg)
            end
        end
        show(dlg)
    end
    
    signal_connect(select_test_images_button, "clicked") do widget
        select_test_images(widget)
    end
    
    test_button = GtkButton("Test")
    push!(vbox, test_button)

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
    batch_processing_checkbox = GtkCheckButton("Batch processing")
    push!(vbox, dust_correction_checkbox)
    push!(vbox, batch_processing_checkbox)

    done_button = GtkButton("Done")
    push!(vbox, done_button)

    function on_done(button)
        # Read value from spin button
        fixed_thresh = Gtk4.value(spin_button) 
        
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
            "images_directory" => directories,
            "Imin_path" => Imins,
            "Imax_path" => Imaxes,
            "fixed_thresh" => fixed_thresh,
            "dust_correction" => dust_correction,
            "batch_processing" => batch_processing 
        )
        
        open("experiment_config.json", "w") do f
            write(f, JSON.json(config, 4))
        end
    end

    if !isinteractive()
        c = Condition()
        signal_connect(done_button, "clicked") do widget
            on_done(widget)
            notify(c)
        end
        @async Gtk4.GLib.glib_main()
        wait(c)
    else
        signal_connect(done_button, "clicked") do widget
            on_done(widget)
            close(win)
        end
    end

    show(win)
end

main()
