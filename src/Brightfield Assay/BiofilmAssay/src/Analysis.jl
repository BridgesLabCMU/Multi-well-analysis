using Images: imfilter, mapwindow, adjust_histogram!, LinearStretching, Gray, N0f16, N0f8,
              Kernel, warp, axes, KernelFactors, RGB 
using ImageMorphology: label_components, component_lengths
using StatsBase: mean, median, quantile
using TiffImages: load, save
using NaturalSort: sort, natural
using DataFrames: DataFrame
using CSV: write
using JSON: parsefile
using IntegralArrays: IntegralArray
using CoordinateTransformations: Translation
using IntervalSets: width, leftendpoint, rightendpoint, Interval, ±
using SubpixelRegistration: phase_offset
using FLoops

round_odd(x) = div(x, 2) * 2 + 1
compmax(x) = length(x) > 1 ? maximum(x[1:end]) : 0

function read_images!(directory, height, width, 
                        ntimepoints, arr, files)
    @inbounds for t in 1:ntimepoints
        @views file = files[t]
        arr[:, :, t] = load("$directory/$file")
    end
    return nothing 
end

function dust_correct!(masks)
    rows, cols, nframes = size(masks)
    for col in 1:cols
        for row in 1:rows
            if masks[row, col, 1]  
                should_suppress = false
                for t in 2:nframes
                    if !masks[row, col, t]
                        should_suppress = true
                        break
                    end
                end
                if should_suppress
					masks[row, col, :] .= false
                end
            end
        end
    end
    return nothing
end

function normalize_local_contrast(img, img_copy, blockDiameter)
	img_copy = imfilter(img_copy, Kernel.gaussian(blockDiameter))
	img = img - img_copy
    return img 
end

function crop(img_stack)
    @views mask = .!any(isnan.(img_stack), dims=3)[:,:,1]
    @views mask_i = any(mask, dims=2)[:,1]
    @views mask_j = any(mask, dims=1)[1,:]
    i1 = findfirst(mask_i)
    i2 = findlast(mask_i)
    j1 = findfirst(mask_j)
    j2 = findlast(mask_j)
    cropped_stack = img_stack[i1:i2, j1:j2, :]
    return cropped_stack, (i1, i2, j1, j2)
end

function stack_preprocess(img_stack, normalized_stack, registered_stack, blockDiameter, nframes, mxshift, sig, Imins, Imaxes)       
    shifts = (0.0, 0.0) 
    @inbounds for t in 1:nframes
        img = img_stack[:,:,t]
        img_copy = img_stack[:,:,t] 
        img_normalized = normalize_local_contrast(img, img_copy, 
                                    blockDiameter)
        normalized_stack[:,:,t] = imfilter(img_normalized, Kernel.gaussian(sig))
        if t == 1
            registered_stack[:,:,t] = normalized_stack[:,:,t]
        else
            moving = normalized_stack[:,:,t]
            fixed = normalized_stack[:,:,t-1]
            shift, _, _ = phase_offset(fixed, moving, upsample_factor=1)
            if sqrt(shift[1]^2 + shift[2]^2) >= mxshift
                shift = Translation(shifts[1], shifts[2])
                registered_stack[:,:,t] = warp(moving, shift, axes(fixed))
                img_stack[:,:,t] = warp(img_stack[:,:,t], shift, axes(fixed))
            else
                shift = Tuple([-1*shift[1], -1*shift[2]])
                shift = shift .+ shifts
                shifts = shift
                shift = Translation(shift[1], shift[2])
                registered_stack[:,:,t] = warp(moving, shift, axes(fixed))
                img_stack[:,:,t] = warp(img_stack[:,:,t], shift, axes(fixed))
            end
        end
    end
    processed_stack, crop_indices = crop(registered_stack)
    row_min, row_max, col_min, col_max = crop_indices
    img_stack = img_stack[row_min:row_max, col_min:col_max, :]
	Imin_image = nothing
	Imax_image = nothing
    if Imins != nothing
        Imin_image = Imins[1][row_min:row_max, col_min:col_max]
    end
    if Imaxes != nothing
        Imax_image = Imaxes[1][row_min:row_max, col_min:col_max]
    end
    return img_stack, processed_stack, Imin_image, Imax_image
end

function compute_mask!(stack, masks, fixed_thresh, ntimepoints)
    @inbounds for t in 1:ntimepoints
        @views masks[:,:,t] = (1 .- stack[:,:,t]) .> fixed_thresh
    end
end

function output_images!(images, stack, masks, overlay, Imin, Imax, dir, filename)
    stack[stack .< 0] .= 0
	stack = Gray{N0f8}.(stack)
    save("$dir/Processed images/$filename.tif", stack)
    @inbounds for i in CartesianIndices(stack)
        gray_val = RGB{N0f8}(stack[i], stack[i], stack[i])
        overlay[i] = masks[i] ? RGB{N0f8}(0,1,1) : gray_val
    end
    save("$dir/Processed images/$filename"*"mask.tif", overlay)
    if length(size(stack)) == 2
        save("$dir/Processed images/$filename"*"OD_image.tif", (-1 .* log10.((images .- Imin) ./ (Imax .- Imin))) .* masks)
    else
        OD_images = similar(images)
        for t in 1:size(stack,3)
            OD_images[:,:,t] = (-1 .* log10.((images[:,:,t] .- Imin) ./ (Imax .- Imin))) .* masks[:,:,t]
        end
        save("$dir/Processed images/$filename"*"OD_image.tif", OD_images)
    end
end

function extract_base_and_ext(filename::String)
    matches = match(r"^(.*\D)\d+(\.[^\.]+)$", filename)
    if isnothing(matches)
        error("Filename $filename does not satisfy batch criteria")
    end
    base = matches.captures[1]
    ext = matches.captures[2]
    return base, ext
end

function timelapse_processing(images, blockDiameter, ntimepoints, shift_thresh, fixed_thresh, sig, dust_correction, dir, filename, Imins, Imaxes)
    images = Float32.(images)
    normalized_stack = similar(images)
    registered_stack = similar(images)
    images, output_stack, Imin, Imax = stack_preprocess(images, normalized_stack, registered_stack, blockDiameter, ntimepoints, shift_thresh, sig, Imins, Imaxes)
    masks = zeros(Bool, size(images))
    compute_mask!(output_stack, masks, fixed_thresh, ntimepoints)
    if dust_correction == "True"
        dust_correct!(masks)
    end
    overlay = zeros(RGB{N0f8}, size(output_stack)...)
    biomasses = zeros(Float32, ntimepoints)
    let images = images
        if Imin != nothing
            if Imax == nothing
				Imax = images[:,:,1]
			end
            output_images!(images, output_stack, masks, overlay, Imin, Imax, dir, filename)
            @floop for t in 1:ntimepoints
                @inbounds biomasses[t] = @views mean((-1 .* log10.((images[:,:,t] .- Imin) ./ (Imax .- Imin))) .* masks[:,:,t])
            end
        else
            @floop for t in 1:ntimepoints
                @inbounds biomasses[t] = @views mean((1 .- images[:,:,t]) .* masks[:,:,t])
            end
        end
    end
    return biomasses
end

function image_processing(image, blockDiameter, fixed_thresh, sig, dir, filename, Imins, Imaxes)
    ntimepoints = 1
    image = Float32.(image)
    img_copy = image 
    img_normalized = normalize_local_contrast(image, img_copy, blockDiameter)
    normalized_blurred = imfilter((1 .- img_normalized), Kernel.gaussian(sig))
    mask = normalized_blurred .> fixed_thresh
    overlay = zeros(RGB{N0f8}, size(image)...)
    if Imins != nothing && Imaxes != nothing
        for Imin in Imins
			if size(Imin) == size(image)
				for Imax in Imaxes
					if size(Imax) == size(image)
                        output_images!(image, img_normalized, mask, overlay, Imin, Imax,dir, filename)
						biomass = mean((-1 .* log10.((image .- Imin) ./ (Imax .- Imin))) .* mask)
						return biomass
					end
				end
			end
		end
    else
        biomass = mean((1 .- image) .* mask)
        return biomass
    end
end

function main()
    config = parsefile("experiment_config.json")
    images_directories  = config["images_directory"]
    dust_correction = config["dust_correction"]
    batch = config["batch_processing"]
    fixed_thresh = config["fixed_thresh"] 
    Imin_paths = config["Imin_path"]
    Imax_paths = config["Imax_path"]
    sig = 2
    blockDiameter = 101 
    shift_thresh = 50

    Imins = nothing
    if length(Imin_paths) > 0
        Imins = [load(Imin_path) for Imin_path in Imin_paths]
    end
    
    Imaxes = nothing
    if length(Imax_paths) > 0
        Imaxes = [load(Imax_path) for Imax_path in Imax_paths]
    end

    @inbounds for k in eachindex(images_directories)
        dir = images_directories[k]
        if isdir("$dir/Processed images")
            rm("$dir/Processed images"; recursive = true)
        end
        if isdir("$dir/Numerical data")
            rm("$dir/Numerical data"; recursive = true)
        end
        mkdir("$dir/Processed images")
        mkdir("$dir/Numerical data")
        BF_output_file = "$dir/Numerical data/biomass.csv"

        analyzed = []
        biomass_data = []
        columns = []
		files = [f for f in readdir(dir) if occursin(r"\.tif$", f)]
        for file in files
            if file ∉ analyzed
                test_image = load("$dir/$file"; lazyio=true)
                img_dims = size(test_image)
                if length(filter(x -> x != 1, img_dims)) == 3
                    height, width, ntimepoints = img_dims
                    images = load("$dir/$file")
                    target_base, target_ext = extract_base_and_ext(file)
                    push!(columns, target_base)
                    push!(biomass_data, timelapse_processing(images, 
                                                             blockDiameter,
                                                             ntimepoints,
                                                             shift_thresh,
                                                             fixed_thresh,
                                                             sig,
                                                             dust_correction,
                                                             dir, target_base, Imins, Imaxes))
                    push!(analyzed, file)
                elseif length(filter(x -> x != 1, img_dims)) == 2
                    target_base, target_ext = extract_base_and_ext(file)
                    matching_files = filter(file_name -> begin
                        base, ext = extract_base_and_ext(file_name)
                        base == target_base && ext == target_ext
                    end, files)
                    push!(columns, target_base)
                    if length(matching_files) > 1 && batch == "True"
                        timelapse_files = sort(matching_files, lt=natural)
                        ntimepoints = length(timelapse_files)
                        height, width = img_dims
                        images = Array{eltype(test_image), 3}(undef, height, width, ntimepoints)
                        read_images!(dir, height, width, ntimepoints, images, timelapse_files)
                        push!(biomass_data, timelapse_processing(images, 
                                                                 blockDiameter,
                                                                 ntimepoints,
                                                                 shift_thresh,
                                                                 fixed_thresh,
                                                                 sig,
                                                                 dust_correction,
                                                                 dir, target_base, Imins, Imaxes))
                        for f in matching_files
                            push!(analyzed, f)
                        end
                    else
                        height, width = img_dims
                        image = load("$dir/$file")
                        push!(biomass_data, image_processing(image, 
                                                             blockDiameter,
                                                             fixed_thresh, sig,
                                                             dir, target_base, Imins, Imaxes))
                        push!(analyzed, file)
                    end
                else
                    error("Number of image dimensions must be either 2 or 3")
                end
            end
        end
		max_length = maximum(length, biomass_data)
		padded = [vcat(l, fill(Missing, max_length - length(l))) for l in biomass_data]
		df = DataFrame(padded, Symbol.(columns))
        write(BF_output_file, df)
    end # loop over directories
end
main()
