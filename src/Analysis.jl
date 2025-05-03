round_odd(x) = isodd(x) ? x : x + 1 
compmax(x) = length(x) > 1 ? maximum(x[1:end]) : 0

function phase_offset(source::AbstractArray, target::AbstractArray; kwargs...)
    plan = plan_fft(source)
    return phase_offset(plan, plan * source, plan * target; kwargs...)
end

function phase_offset(
    plan,
    source_freq::AbstractMatrix{<:Complex{T}},
    target_freq;
    upsample_factor = 1,
    normalize = false,
) where {T}
    image_product = @. source_freq * conj(target_freq)
    if normalize
        @. image_product /= max(abs(image_product), eps(T))
    end
    if isone(upsample_factor)
        cross_correlation = ifft!(image_product)
    else
        cross_correlation = plan \ image_product
    end
    maxima, maxidx = @compat findmax(abs, cross_correlation)
    shape = size(source_freq)
    midpoints = map(ax -> (first(ax) + last(ax)) / T(2), axes(source_freq))
    idxoffset = map(first, axes(cross_correlation))
    shift = @. T(ifelse(maxidx.I > midpoints, maxidx.I - shape, maxidx.I) - idxoffset)

    isone(upsample_factor) &&
        return (; shift, calculate_stats(maxima, source_freq, target_freq)...)

    shift = @. round(shift * upsample_factor) / T(upsample_factor)
    upsample_region_size = ceil(upsample_factor * T(1.5))
    dftshift = div(upsample_region_size, 2)
    sample_region_offset = @. dftshift - shift * upsample_factor
    cross_correlation = upsampled_dft(
        image_product,
        upsample_region_size,
        upsample_factor,
        sample_region_offset,
    )
    maxima, maxidx = @compat findmax(abs, cross_correlation)
    shift = @. shift + (maxidx.I - dftshift - idxoffset) / T(upsample_factor)

    stats = calculate_stats(maxima, source_freq, target_freq)
    return (; shift, stats...)
end

function upsampled_dft(
    data::AbstractMatrix{T},
    region_size,
    upsample_factor,
    offsets,
) where {T<:Complex}
    shiftrange = 1:region_size
    idxoffset = map(first, axes(data))
    sample_rate = inv(T(upsample_factor))
    freqs = fftfreq(size(data, 2), sample_rate)
    kernel = @. cis(-T(2π) * (shiftrange - offsets[2] - idxoffset[2]) * freqs')

    _data = kernel * data'

    freqs = fftfreq(size(data, 1), sample_rate)
    kernel = @. cis(T(2π) * (shiftrange - offsets[1] - idxoffset[1]) * freqs')
    _data = kernel * _data'
    return _data
end

function calculate_stats(crosscor_maxima, source_freq, target_freq)
    source_amp = mean(abs2, source_freq)
    target_amp = mean(abs2, target_freq)
    error = 1 - abs2(crosscor_maxima) / (source_amp * target_amp)
    phasediff = atan(imag(crosscor_maxima), real(crosscor_maxima))
    return (; error, phasediff)
end

function read_images!(ntimepoints, arr, files)
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

function mean_filter!(X, length_scale)
    iX = IntegralArray(X)
    @inbounds for i in CartesianIndex(1,1):CartesianIndex(size(X))
        x, y = i.I
        x_int = x±length_scale
        y_int = y±length_scale
        x_int = Interval(max(leftendpoint(x_int), 1),
                         min(rightendpoint(x_int), size(X)[1]))
        y_int = Interval(max(leftendpoint(y_int), 1),
                         min(rightendpoint(y_int), size(X)[2]))
        X[i] = iX[x_int, y_int]/(IntervalSets.width(x_int)*IntervalSets.width(y_int))
    end
    return nothing
end

function normalize_local_contrast_output(normalized, images, images_copy, blockDiameter, fpMean)
	length_scale = Int((blockDiameter-1)/2)
    if length(size(images)) == 3
        for t in 1:size(images, 3)
            img = images[:,:,t]
            img_copy = images_copy[:,:,t]
            mean_filter!(img_copy, length_scale)
            img = img - img_copy
            img .+= fpMean
            @. img[img < 0.0] = 0.0
            @. img[img > 1.0] = 1.0
            normalized[:,:,t] = img  
        end
    elseif length(size(images)) == 2
        img_copy = images_copy
        mean_filter!(img_copy, length_scale)
        images = images - img_copy
        images .+= fpMean
        @. images[images < 0.0] = 0.0
        @. images[images > 1.0] = 1.0
        normalized .= images
    else
        error("Dimension error in normalize_local_contrast_output")
    end
    return normalized 
end

function normalize_local_contrast(img, img_copy, blockDiameter)
	img = 1 .- img
	img_copy = 1 .- img_copy
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

function stack_preprocess(img_stack, normalized_stack, registered_stack, blockDiameter, nframes, mxshift, sig, Imin, Imax)       
    shifts_array = []
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
                push!(shifts_array, shift)
                registered_stack[:,:,t] = warp(moving, shift, axes(fixed))
                img_stack[:,:,t] = warp(img_stack[:,:,t], shift, axes(fixed))
            else
                shift = Tuple([-1*shift[1], -1*shift[2]])
                shift = shift .+ shifts
                shifts = shift
                shift = Translation(shift[1], shift[2])
                push!(shifts_array, shift)
                registered_stack[:,:,t] = warp(moving, shift, axes(fixed))
                img_stack[:,:,t] = warp(img_stack[:,:,t], shift, axes(fixed))
            end
        end
    end
    processed_stack, crop_indices = crop(registered_stack)
    row_min, row_max, col_min, col_max = crop_indices
    img_stack = img_stack[row_min:row_max, col_min:col_max, :]
    if Imin != nothing
        Imin = Imin[row_min:row_max, col_min:col_max]
    end
    if Imax != nothing
        Imax = Imax[row_min:row_max, col_min:col_max]
    end
    return img_stack, processed_stack, Imin, Imax, crop_indices, shifts_array
end

function write_OD_images!(OD_images, dir, filename)
    if length(size(OD_images)) == 3 
        for t in 1:size(OD_images, 3)
            save("$dir/Processed images/$filename"*"_OD"*"_t"*string(t)*".tif", 
                        OD_images[:,:,t])
        end
    elseif length(size(OD_images)) == 2
        save("$dir/Processed images/$filename"*"_OD.tif", 
                    OD_images)
    end
end

function output_images!(stack, masks, overlay, OD_images, dir, filename, blockDiameter)
    filename = split(filename, "/")[end]
    normalized = similar(stack)
	fpMax = maximum(stack)
	fpMin = minimum(stack)
	fpMean = (fpMax - fpMin) / 2.0 + fpMin
	normalized = normalize_local_contrast_output(normalized, stack, copy(stack), blockDiameter, fpMean)
	normalized = Gray{N0f8}.(normalized)
    save("$dir/Processed images/$filename.tif", normalized)
    @inbounds for i in CartesianIndices(normalized)
        gray_val = RGB{N0f8}(normalized[i], normalized[i], normalized[i])
        overlay[i] = masks[i] ? RGB{N0f8}(0,1,1) : gray_val
    end
    save("$dir/Processed images/$filename"*"mask.tif", overlay)
    if OD_images != nothing
        write_OD_images!(OD_images, dir, filename)
    end
end

function extract_base_and_ext(filename::String, batch)
    # The batch criterion is that there is a number at the end of the filename
    if batch == "True"
        matches = match(r"^(.*\D)\d+(\.[^\.]+)$", filename)
        if isnothing(matches)
            base, ext = splitext(filename)
        else
            base = matches.captures[1]
            ext = matches.captures[2]
        end
    else
        base, ext = splitext(filename)
    end
    return base, ext
end

function timelapse_processing(images, blockDiameter, ntimepoints, shift_thresh, fixed_thresh, sig, dust_correction, dir, filename, Imin, Imax)
    images = Float64.(images)
    normalized_stack = similar(images)
    registered_stack = similar(images)
    images, output_stack, Imin, Imax, crop_indices, shifts = stack_preprocess(images, normalized_stack, registered_stack, blockDiameter, ntimepoints, shift_thresh, sig, Imin, Imax)
    masks = zeros(Bool, size(images))
    compute_mask!(output_stack, masks, fixed_thresh, ntimepoints)
    if dust_correction == "True"
        dust_correct!(masks)
    end
    overlay = zeros(RGB{N0f8}, size(output_stack)...)
    biomasses = zeros(Float64, ntimepoints)
    OD_images = nothing
    if Imin != nothing
        OD_images = Array{Gray{Float32}, 3}(undef, size(images))
        for t in 1:ntimepoints
            if Imax != nothing
                OD_images[:,:,t] = @views (-1 .* log10.((images[:,:,t] .- Imin) ./ (Imax .- Imin)))
            else
                OD_images[:,:,t] = @views (-1 .* log10.((images[:,:,t] .- Imin) ./ (images[:,:,1] .- Imin)))
            end
            @inbounds biomasses[t] = @views Float64(mean(OD_images[:,:,t] .* masks[:,:,t]))
        end
    else
        for t in 1:ntimepoints
            @inbounds biomasses[t] = @views Float64(mean((1 .- images[:,:,t]) .* masks[:,:,t]))
        end
    end
    output_images!(images, masks, overlay, OD_images, dir, filename, blockDiameter)
    return shifts, crop_indices, masks, biomasses
end

function image_processing(image, blockDiameter, fixed_thresh, sig, dir, filename, Imin, Imax)
    ntimepoints = 1
    image = Float64.(image)
    img_copy = image 
    img_normalized = normalize_local_contrast(image, img_copy, blockDiameter)
    normalized_blurred = imfilter(img_normalized, Kernel.gaussian(sig))
    mask = normalized_blurred .> fixed_thresh
    overlay = zeros(RGB{N0f8}, size(image)...)
    OD_image = nothing
    biomass = nothing
    if Imin != nothing && Imax != nothing
        OD_image = Array{Gray{Float32}, 2}(undef, size(image))
        OD_image .= (-1 .* log10.((image .- Imin) ./ (Imax .- Imin)))
        biomass = Float64(mean(OD_image .* mask))
    else
        biomass = Float64(mean((1 .- image) .* mask))
    end
    output_images!(image, mask, overlay, OD_image, dir, filename, blockDiameter)
    return mask, biomass
end

function filter_same_well(target, paths)
    fn = basename(target)
    well = first(split(fn, '_'))
    return filter(p -> first(split(basename(p), '_')) == well, paths)
end

function register_and_crop(raw_stack, shifts_array, crop_indices) 
    n_rows, n_cols, n_frames = size(raw_stack)
    registered = similar(raw_stack)
    base_axes = axes(raw_stack[:,:,1])
    registered[:,:,1] = raw_stack[:,:,1]
    for t in 2:n_frames
        registered[:,:,t] = warp(raw_stack[:,:,t], shifts_array[t-1], base_axes)
    end
    row_min, row_max, col_min, col_max = crop_indices
    cropped_stack = registered[row_min:row_max, col_min:col_max, :]
    return cropped_stack
end

function analysis_main()
    config = JSON.parsefile("experiment_config.json")
    images_directories  = config["images_directory"]
    dust_correction = config["dust_correction"]
    batch = config["batch_processing"]
    fixed_thresh = config["fixed_thresh"] 
    Imin_path = config["Imin_path"]
    Imax_path = config["Imax_path"]
    blockDiameter = round_odd(config["blockDiam"]) 
    shift_thresh = 50
    sig = 2

    # Read in protocol
    if !isfile(joinpath(images_directories[1], "protocol.csv"))
        error("Protocol file not found when running Analysis.jl")
    end
    protocol = CSV.File(joinpath(images_directories[1], "protocol.csv"), header=true)

    Imin = nothing
    if Imin_path != ""
        Imin = load(Imin_path)
    end
    
    Imax = nothing
    if Imax_path != ""
        Imax = load(Imax_path)
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

		tif_paths = [f for f in readdir(dir, join=true) if occursin(".tif", f)]
		bf_reads = filter(r -> 
			r.action    == "Imaging Read" &&
			r.channel   == "Bright Field", protocol)

		organized_files = []
        mags = []
		for bf in eachrow(bf_reads)
            bf = bf[1]
			mag = bf.magnification
			same_mag_steps = 
			  protocol.step[
				(protocol.action .== "Imaging Read") .& 
				(protocol.magnification .== mag)
			  ]
            same_mag_channel = protocol.channel[
                (protocol.action .== "Imaging Read") .& 
                (protocol.magnification .== mag)
            ]

			mag_files = [] 
            for (l,s) in enumerate(same_mag_steps)
				pat = "_0$(s)_"
                append!(mag_files, filter(p -> (occursin(pat, p) && occursin(same_mag_channel[l], p)), tif_paths))
			end
            push!(mags, mag)   
			push!(organized_files, mag_files)
		end

        for i in eachindex(organized_files) 
            @views mag = mags[i]
            @views bf_files = [f for f in organized_files[i] if occursin("Bright Field", f)]
            BF_output_file = "$dir/Numerical data/$(mag)_BF_biomass.csv"
            analyzed = []
            biomass_data = []
            columns = []

            @views non_bf_files = [f for f in organized_files[i] if !occursin("Bright Field", f)]
			CH_RE = r"_0(\d+)_(\d+)_([^_]+)_(\d+)_"

			groups = Tuple{Int,Int,String,Int}[]
			for f in non_bf_files
				m = match(CH_RE, f)
				if m === nothing
					@warn "skipping non-matching file" f
					continue
				end
				i1 = parse(Int, m.captures[1])
				i2 = parse(Int, m.captures[2])
				ch =       m.captures[3]
				i3 = parse(Int, m.captures[4])
				key = (i1, i2, ch, i3)
				if key ∉ groups
					push!(groups, key)
				end
			end

			channels = [ ch for (_, _, ch, _) in groups ]
			files_by_channel = [
				filter(f -> begin
					m = match(CH_RE, f)
					( parse(Int,m.captures[1]),
					  parse(Int,m.captures[2]),
					  m.captures[3],
					  parse(Int,m.captures[4]) ) == key
				end, non_bf_files)
				for key in groups
			]

			other_output_files = [
				joinpath(dir, "Numerical data",
						 string(mag, "_", i1, "_", i2, "_", ch, "_", i3, "_biomass.csv"))
				for (i1, i2, ch, i3) in groups
			]
			
			non_bf_data     = [ [] for _ in groups ]
			channel_columns = [ [] for _ in groups ]

            for file in bf_files
                if file ∉ analyzed
                    test_image = load(file; lazyio=true)
                    img_dims = size(test_image)
                    if length(filter(x -> x != 1, img_dims)) == 3
                        height, width, ntimepoints = img_dims
                        images = load(file)
                        target_base, target_ext = extract_base_and_ext(file, batch)
                        push!(columns, target_base)
                        shifts, crop_indices, mask, biomass = timelapse_processing(images, 
                                                                 blockDiameter,
                                                                 ntimepoints,
                                                                 shift_thresh,
                                                                 fixed_thresh,
                                                                 sig,
                                                                 dust_correction,
                                                                 dir, target_base, Imin, Imax)
                        push!(biomass_data, biomass)
                        push!(analyzed, file)
                        # Check for corresponding non-bf files
                        # if they exist, warp and crop, apply mask, calculate average signal, output
                        for j in eachindex(files_by_channel)
                            channel = channels[j]
                            # Check if there are channel images corresponding to the same well
                            channel_files = filter_same_well(file, files_by_channel[j])
                            if length(channel_files) > 0
                                channel_file = channel_files[1]
                                images = load(channel_file)
                                target_base, target_ext = extract_base_and_ext(channel_file, batch)
                                push!(channel_columns[j], target_base)
                                # Register and crop images
                                channel_registered = register_and_crop(images, shifts, crop_indices)
                                channel_float = Float64.(channel_registered)
                                # Apply mask, perform background subtraction, calculate average, output registered images
                                n_frames = size(channel_registered, 3)
                                channel_biomass = zeros(Float64, n_frames) 
                                for t in 1:n_frames
                                    @views background = mean(channel_float[:,:,t][.!mask[:,:,t]])
                                    signal = @views !any(mask[:,:,t]) ? 0 : mean((channel_float[:,:,t] .- background).*mask[:,:,t])
                                    channel_biomass[t] = signal
                                end
                                push!(non_bf_data[j], channel_biomass)
                                channel_filename = split(target_base, "/")[end]
                                save("$dir/Processed images/$channel_filename.tif", channel_registered)
                            end
                        end
                    elseif length(filter(x -> x != 1, img_dims)) == 2
                        target_base, target_ext = extract_base_and_ext(file, batch)
                        matching_files = filter(file_name -> begin
                            base, ext = extract_base_and_ext(file_name, batch)
                            base == target_base && ext == target_ext
                        end, bf_files)
                        push!(columns, target_base)
                        if length(matching_files) > 1 && batch == "True"
                            timelapse_files = sort(matching_files, lt=natural)
                            ntimepoints = length(timelapse_files)
                            height, width = img_dims
                            images = Array{eltype(test_image), 3}(undef, height, width, ntimepoints)
                            read_images!(ntimepoints, images, timelapse_files)
                            shifts, crop_indices, mask, biomass = timelapse_processing(images, 
                                                                     blockDiameter,
                                                                     ntimepoints,
                                                                     shift_thresh,
                                                                     fixed_thresh,
                                                                     sig,
                                                                     dust_correction,
                                                                     dir, target_base, Imin, Imax)
                            push!(biomass_data, biomass)
                            for f in matching_files
                                push!(analyzed, f)
                            end
                            # Check for corresponding non-bf files
                            # if they exist, warp and crop, apply mask, calculate average signal, output
                            for j in eachindex(files_by_channel)
                                channel = channels[j]
                                # Check if there are channel images corresponding to the same well
                                channel_files = filter_same_well(file, files_by_channel[j])
                                if length(channel_files) > 0
                                    channel_file = channel_files[1]
                                    channel_test = load(channel_file; lazyio=true)
                                    channel_dims = size(channel_test)
                                    ntimepoints = length(channel_files)
                                    images = Array{eltype(channel_test), 3}(undef, height, width, ntimepoints)
                                    read_images!(ntimepoints, images, channel_files)
                                    target_base, target_ext = extract_base_and_ext(channel_file, batch)
                                    push!(channel_columns[j], target_base)
                                    # Register and crop images
                                    channel_registered = register_and_crop(images, shifts, crop_indices)
                                    channel_float = Float64.(channel_registered)
                                    # Apply mask, perform background subtraction, calculate average, output registered images
                                    n_frames = size(channel_registered, 3)
                                    channel_biomass = zeros(Float64, n_frames) 
                                    for t in 1:n_frames
                                        @views background = mean(channel_float[:,:,t][.!mask[:,:,t]])
                                        signal = @views !any(mask[:,:,t]) ? 0 : mean((channel_float[:,:,t] .- background).*mask[:,:,t])
                                        channel_biomass[t] = signal
                                    end
                                    push!(non_bf_data[j], channel_biomass)
                                    channel_filename = split(target_base, "/")[end]
                                    save("$dir/Processed images/$channel_filename.tif", channel_registered)
                                end
                            end
                        else
                            height, width = img_dims
                            image = load(file)
                            mask, biomass = image_processing(image, 
                                                         blockDiameter,
                                                         fixed_thresh, sig,
                                                         dir, target_base, Imin, Imax)
                            push!(biomass_data, biomass)
                            push!(analyzed, file)
                            for j in eachindex(files_by_channel)
                                channel = channels[j]
                                channel_files = filter_same_well(file, files_by_channel[j])
                                if length(channel_files) > 0
                                    channel_file = channel_files[1]
                                    image = load(channel_file)
                                    image_float = Float64.(image)
                                    target_base, target_ext = extract_base_and_ext(channel_file, batch)
                                    push!(channel_columns[j], target_base)
                                    @views background = mean(image_float[.!mask])
                                    channel_biomass = !any(mask) ? 0 : mean((image_float .- background).*mask)
                                    push!(non_bf_data[j], channel_biomass)
                                    channel_filename = split(target_base, "/")[end]
                                    save("$dir/Processed images/$channel_filename.tif", image)
                                end
                            end
                        end
                    else
                        error("Number of image dimensions must be either 2 or 3")
                    end
                end
            end
            if length(biomass_data) > 0
                max_length = maximum(length, biomass_data)
                padded = [vcat(l, fill(Missing, max_length - length(l))) for l in biomass_data]
                df = DataFrame(padded, Symbol.(columns))
                write(BF_output_file, df)
            end
            for i in eachindex(non_bf_data)
                channel_data = non_bf_data[i]
                if length(channel_data) > 0
                    max_length = maximum(length, channel_data)
                    padded = [vcat(l, fill(Missing, max_length - length(l))) for l in channel_data]
                    df = DataFrame(padded, Symbol.(channel_columns[i]))
                    write(other_output_files[i], df)
                end
            end
        end # loop over mags
    end # loop over directories
end
