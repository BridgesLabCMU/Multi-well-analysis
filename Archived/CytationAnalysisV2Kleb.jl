using Images: imfilter, mapwindow, adjust_histogram!, LinearStretching, Gray, N0f16, N0f8,
              Kernel, warp, axes, KernelFactors, RGB 
using ImageMorphology: label_components, component_lengths
using StatsBase: mean, median, quantile
using TiffImages: load, save
using NaturalSort: sort, natural
using DataFrames
using CSV
using JSON: parsefile
using IntegralArrays: IntegralArray
using CoordinateTransformations: Translation
using IntervalSets: width, leftendpoint, rightendpoint, Interval, ±
using SubpixelRegistration: phase_offset
using FLoops

round_odd(x) = div(x, 2) * 2 + 1
compmax(x) = length(x) > 1 ? maximum(x[1:end]) : 0

function read_images!(well, directory, height, width, 
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
        X[i] = iX[x_int, y_int]/(width(x_int)*width(y_int))
    end
    return nothing
end

function normalize_local_contrast_output(normalized, images, images_copy, blockDiameter, fpMean)
	length_scale = Int((blockDiameter-1)/2)
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

function crop_square(images, y1, y2, x1, x2)
    height, width, _ = size(images)
    y2 = min(y2, height)
    x2 = min(x2, width)
    crop_height = y2 - y1 + 1
    crop_width = x2 - x1 + 1
    square_size = min(crop_height, crop_width)
    y2 = y1 + square_size - 1
    x2 = x1 + square_size - 1
    return images[y1:y2, x1:x2, :]
end

function stack_preprocess(img_stack, normalized_stack, registered_stack, blockDiameter, nframes, mxshift, sig, CFP_stack, YFP_stack, texas_red_stack, CY5_stack)       
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
				if CFP_stack != nothing
					CFP_stack[:,:,t] = warp(CFP_stack[:,:,t], shift, axes(fixed))
				end
				if YFP_stack != nothing
					YFP_stack[:,:,t] = warp(YFP_stack[:,:,t], shift, axes(fixed))
				end
				if texas_red_stack != nothing
					texas_red_stack[:,:,t] = warp(texas_red_stack[:,:,t], shift, axes(fixed))
				end
				if CY5_stack != nothing
					CY5_stack[:,:,t] = warp(CY5_stack[:,:,t], shift, axes(fixed))
				end
            else
                shift = Tuple([-1*shift[1], -1*shift[2]])
                shift = shift .+ shifts
                shifts = shift
                shift = Translation(shift[1], shift[2])
                registered_stack[:,:,t] = warp(moving, shift, axes(fixed))
                img_stack[:,:,t] = warp(img_stack[:,:,t], shift, axes(fixed))
				if CFP_stack != nothing
					CFP_stack[:,:,t] = warp(CFP_stack[:,:,t], shift, axes(fixed))
				end
				if YFP_stack != nothing
					YFP_stack[:,:,t] = warp(YFP_stack[:,:,t], shift, axes(fixed))
				end
				if texas_red_stack != nothing
					texas_red_stack[:,:,t] = warp(texas_red_stack[:,:,t], shift, axes(fixed))
				end
				if CY5_stack != nothing
					CY5_stack[:,:,t] = warp(CY5_stack[:,:,t], shift, axes(fixed))
				end
            end
        end
    end
    crop1 = round(Int,0.075*size(img_stack, 2))
    crop2 = size(img_stack, 2) - crop1
    processed_stack, crop_indices = crop(registered_stack)
    row_min, row_max, col_min, col_max = crop_indices
    img_stack = img_stack[row_min:row_max, col_min:col_max, :]
    if CFP_stack != nothing
        CFP_stack = CFP_stack[row_min:row_max, col_min:col_max, :]
        CFP_stack = crop_square(CFP_stack, crop1, crop2, crop1, crop2)
    end
    if YFP_stack != nothing
        YFP_stack = YFP_stack[row_min:row_max, col_min:col_max, :]
        YFP_stack = crop_square(YFP_stack, crop1, crop2, crop1, crop2)
    end
    if texas_red_stack != nothing
        texas_red_stack = texas_red_stack[row_min:row_max, col_min:col_max, :]
        texas_red_stack = crop_square(texas_red_stack, crop1, crop2, crop1, crop2)
    end
    if CY5_stack != nothing
        CY5_stack = CY5_stack[row_min:row_max, col_min:col_max, :]
        CY5_stack = crop_square(CY5_stack, crop1, crop2, crop1, crop2)
    end
    img_stack = crop_square(img_stack, crop1, crop2, crop1, crop2)
    processed_stack = crop_square(processed_stack, crop1, crop2, crop1, crop2)
    return CFP_stack, YFP_stack, texas_red_stack, CY5_stack, img_stack, processed_stack
end

function compute_mask!(stack, masks, fixed_thresh, ntimepoints)
    @inbounds for t in 1:ntimepoints
		@views masks[:,:,t] = stack[:,:,t] .> fixed_thresh
    end
end

function output_images!(stack, masks, overlay, dir, well)
    normalized = similar(stack)
	fpMax = maximum(stack)
	fpMin = minimum(stack)
	fpMean = (fpMax - fpMin) / 2.0 + fpMin
	normalized = normalize_local_contrast_output(normalized, stack, copy(stack), 1001, fpMean) # number is block size (odd number)
	normalized = Gray{N0f8}.(normalized)
    save("$dir/results_images/$well.tif", normalized)
    @inbounds for i in CartesianIndices(normalized)
        gray_val = RGB{N0f8}(normalized[i], normalized[i], normalized[i])
        overlay[i] = masks[i] ? RGB{N0f8}(0,1,1) : gray_val
    end
    save("$dir/results_images/$well"*"mask.tif", overlay)
end

function multiple_mags(dir, bulk_file)
    df = DataFrame(CSV.File(bulk_file, header=false))
    objective_rows = passmissing(occursin).([r"Objective"], df.Column2)  
	objective_rows = [ismissing(x) ? false : x for x in objective_rows]
    objective_rows = findall(objective_rows)
    # Find all rows with "Read" and count how many before the first imaging step
    read_rows = passmissing(==).(["Read"], df.Column1)
    read_rows = [ismissing(x) ? false : x for x in read_rows]
    read_rows = findall(read_rows)
    function generate_flags(n_objectives, count)
        start = 2 + count
        return ["_$(i<10 ? "0$i" : string(i))_1_" for i in start:(start+n_objectives-1)]
    end
    count = 0
	for row in read_rows
		if df[row, 2] == "Image Single Image" || occursin("x", df[row, 2])
			break 
		end
		count += 1
	end
    if length(objective_rows) > 1
        flags = generate_flags(length(objective_rows), count)
        BF_output_files = []
        for i in 1:length(objective_rows)
            objective = match(r"(\d{1,2}x)", df[objective_rows[i], 2]).match
            push!(BF_output_files, "$dir/results_data/BF_imaging_"*objective*".csv")
        end
    else
        if count == 0
            flags = ["_02_1_"]
        elseif count == 1
            flags = ["_03_1_"]
        elseif count == 2
            flags = ["_04_1_"]
        elseif count == 3
            flags = ["_05_1_"]
        end
        BF_output_files = ["$dir/results_data/BF_imaging.csv"]
    end
    return flags, BF_output_files
end

function main()
    config = parsefile("experiment_config.json")
    images_directories  = config["images_directory"]
    all_conditions = config["conditions"]
    sig = config["sig"]
    fixed_thresh = float(config["fixed_thresh"])
    blockDiameter = config["blockDiameter"] 
    shift_thresh = config["shift_thresh"]
    dust_correction = config["dust_correction"]

    @inbounds for k in eachindex(images_directories)
        bulk_file = images_directories[k]*"/metadata.csv" 
        dir = images_directories[k]
        conditions = all_conditions[k]
        if isdir("$dir/results_images")
            rm("$dir/results_images"; recursive = true)
        end
        if isdir("$dir/results_data")
            rm("$dir/results_data"; recursive = true)
        end
        mkdir("$dir/results_images")
        mkdir("$dir/results_data")
        flags, BF_output_files = multiple_mags(dir, bulk_file)
        @inbounds for (l, flag) in enumerate(flags)
            BF_output_file = BF_output_files[l] 
            num_wells = sum(length(v) for v in values(conditions))
            all_wells = vcat(values(conditions)...)
            files = [f for f in readdir(dir) if occursin(r"\.tif$", f) && occursin("Bright Field", f) && any(occursin(well*"_", f) && occursin(flag, f) for well in all_wells)]
            ntimepoints = div(length(files), num_wells)
            file1 = files[1]
            test_image = load("$dir/$file1")
            height, width = size(test_image)
            BF_data_matrix = Array{Float64, 2}(undef, ntimepoints, num_wells)
            CFP_data_matrix = Array{Float64, 2}(undef, ntimepoints, num_wells)
            YFP_data_matrix = Array{Float64, 2}(undef, ntimepoints, num_wells)
            texas_red_data_matrix = Array{Float64, 2}(undef, ntimepoints, num_wells)
            CY5_data_matrix = Array{Float64, 2}(undef, ntimepoints, num_wells)
            conditions = values(conditions)
            CFP_flag = false
            YFP_flag = false
            texas_red_flag = false
            CY5_flag = false
            idx = 0
            @inbounds for (i, wells) in enumerate(conditions)
                @inbounds for j in eachindex(wells)
                    idx += 1
                    CFP_images = nothing
                    YFP_images = nothing
                    texas_red_images = nothing
                    CY5_images = nothing
                    images = Array{Gray{N0f16}, 3}(undef, height, width, ntimepoints)
                    well = wells[j]
                    BF_well_files = sort([f for f in readdir(dir) if occursin(well*"_", f) && occursin("Bright Field", f) && occursin(flag, f)], 
                                 lt=natural)
                    CFP_well_files = sort([f for f in readdir(dir) if occursin(well*"_", f) && occursin("CFP", f) && occursin(flag, f)], 
                                 lt=natural)
                    YFP_well_files = sort([f for f in readdir(dir) if occursin(well*"_", f) && occursin("YFP", f) && occursin(flag, f)], 
                                 lt=natural)
                    texas_red_well_files = sort([f for f in readdir(dir) if occursin(well*"_", f) && occursin("Texas Red", f) && occursin(flag, f)], 
                                 lt=natural)
                    CY5_well_files = sort([f for f in readdir(dir) if occursin(well*"_", f) && occursin("CY5", f) && occursin(flag, f)], 
                                 lt=natural)
                    read_images!(well, dir, height, width, ntimepoints, images, BF_well_files)
                    if length(CFP_well_files) > 0
                        if !CFP_flag
                            CFP_flag = true
                        end
                        CFP_images = Array{Gray{N0f16}, 3}(undef, height, width, ntimepoints)
                        read_images!(well, dir, height, width, ntimepoints, CFP_images, CFP_well_files)
                    end
                    if length(YFP_well_files) > 0
                        if !YFP_flag
                            YFP_flag = true
                        end
                        YFP_images = Array{Gray{N0f16}, 3}(undef, height, width, ntimepoints)
                        read_images!(well, dir, height, width, ntimepoints, YFP_images, YFP_well_files)
                    end
                    if length(texas_red_well_files) > 0
                        if !texas_red_flag
                            texas_red_flag = true
                        end
                        texas_red_images = Array{Gray{N0f16}, 3}(undef, height, width, ntimepoints)
                        read_images!(well, dir, height, width, ntimepoints, texas_red_images, texas_red_well_files)
                    end
                    if length(CY5_well_files) > 0
                        if !CY5_flag
                            CY5_flag = true
                        end
                        CY5_images = Array{Gray{N0f16}, 3}(undef, height, width, ntimepoints)
                        read_images!(well, dir, height, width, ntimepoints, CY5_images, CY5_well_files)
                    end
                    images = Float64.(images)
                    normalized_stack = similar(images)
                    registered_stack = similar(images)
                    ##### CHANGE THRESHOLD #####
                    CFP_images, YFP_images, texas_red_images, CY5_images, images, output_stack = stack_preprocess(images, normalized_stack, 
                                                                                                                  registered_stack, blockDiameter[2], 
                                                                                                                  ntimepoints, shift_thresh, 
                                                                                                                  sig, CFP_images, YFP_images, 
                                                                                                                  texas_red_images, CY5_images)
                    masks = zeros(Bool, size(images))
                    compute_mask!(output_stack, masks, fixed_thresh, ntimepoints)
                    if dust_correction == "True"
                        dust_correct!(masks)
                    end
                    overlay = zeros(RGB{N0f8}, size(output_stack)...)
                    if length(flags) > 1
                        str = BF_output_files[l]
                        well = well*"_"*str[findlast(isequal('_'), str)+1:length(str)-4]
                    end
                    output_images!(images, masks, overlay, dir, well)
                    let images = images
                        @floop for t in 1:ntimepoints
                            @inbounds signal = @views mean((1 .- images[:,:,t]) .* masks[:,:,t])
                            @inbounds BF_data_matrix[t, idx] = signal 
                        end
                    end
                    if CFP_flag
                        let CFP_images = CFP_images
                            @floop for t in 1:ntimepoints
                                @inbounds signal = @views !any(masks[:,:,t]) ? 0 : mean(CFP_images[:,:,t]) - mean(CFP_images[:,:,t][.!masks[:,:,t]])
                                @inbounds CFP_data_matrix[t, idx] = signal 
                            end
                        end
                    end
                    if YFP_flag
                        let YFP_images = YFP_images
                            @floop for t in 1:ntimepoints
                                @inbounds signal = @views !any(masks[:,:,t]) ? 0 : mean(YFP_images[:,:,t]) - mean(YFP_images[:,:,t][.!masks[:,:,t]])
                                @inbounds YFP_data_matrix[t, idx] = signal 
                            end
                        end
                    end
                    if texas_red_flag
                        let texas_red_images = texas_red_images
                            @floop for t in 1:ntimepoints
                                @inbounds signal = @views !any(masks[:,:,t]) ? 0 : mean(texas_red_images[:,:,t]) - mean(texas_red_images[:,:,t][.!masks[:,:,t]])
                                @inbounds texas_red_data_matrix[t, idx] = signal 
                            end
                        end
                    end
                    if CY5_flag
                        let CY5_images = CY5_images
                            @floop for t in 1:ntimepoints
                                @inbounds signal = @views !any(masks[:,:,t]) ? 0 : mean(CY5_images[:,:,t]) - mean(CY5_images[:,:,t][.!masks[:,:,t]])
                                @inbounds CY5_data_matrix[t, idx] = signal 
                            end
                        end
                    end
                end # loop over wells for a condition 
            end # loop over conditions 
            df = DataFrame(BF_data_matrix, Symbol.(all_wells))
            df .= ifelse.(isnan.(df), 0, df)
            CSV.write(BF_output_file, df)
            if CFP_flag
                CFP_output_file = "$dir/results_data/CFP_imaging.csv"
                df = DataFrame(CFP_data_matrix, Symbol.(all_wells))
                df .= ifelse.(isnan.(df), 0, df)
                CSV.write(CFP_output_file, df)
            end
            if YFP_flag
                YFP_output_file = "$dir/results_data/YFP_imaging.csv"
                df = DataFrame(YFP_data_matrix, Symbol.(all_wells))
                df .= ifelse.(isnan.(df), 0, df)
                CSV.write(YFP_output_file, df)
            end
            if texas_red_flag
                texas_red_output_file = "$dir/results_data/texas_red_imaging.csv"
                df = DataFrame(texas_red_data_matrix, Symbol.(all_wells))
                df .= ifelse.(isnan.(df), 0, df)
                CSV.write(texas_red_output_file, df)
            end
            if CY5_flag
                CY5_output_file = "$dir/results_data/CY5_imaging.csv"
                df = DataFrame(CY5_data_matrix, Symbol.(all_wells))
                df .= ifelse.(isnan.(df), 0, df)
                CSV.write(CY5_output_file, df)
            end
        end # loop over flags
    end # loop over directories
end
main()
