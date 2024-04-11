using Images: imfilter, mapwindow, Gray, N0f16, N0f8,
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

function thresh(μ, t₀)
    if μ < t₀/1.5 
        return (t₀+0.03) + ((t₀-0.03)-(t₀+0.03))*μ
    else
        return (t₀+0.06) + ((t₀-0.06)-(t₀+0.06))*μ
    end
end

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

function normalize_local_contrast(img, img_copy, blockDiameter, fpMean)
    length_scale = Int((blockDiameter-1)/2)
    mean_filter!(img_copy, length_scale)
    img -= img_copy 
    img .+= fpMean 
    @. img[img < 0.0] = 0.0
    @. img[img > 1.0] = 1.0
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

function stack_preprocess(img_stack, normalized_stack, registered_stack,
                        blockDiameter, fpMean, nframes, mxshift, sig,
                        CFP_stack, YFP_stack, texas_red_stack, CY5_stack)       
    shifts = (0.0, 0.0) 
    @inbounds for t in 1:nframes
        img = img_stack[:,:,t]
        img_copy = img_stack[:,:,t] 
        img_normalized = normalize_local_contrast(img, img_copy, 
                                    blockDiameter, fpMean)
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
    if CFP_stack != nothing
        CFP_stack = CFP_stack[row_min:row_max, col_min:col_max, :]
    end
    if YFP_stack != nothing
        YFP_stack = YFP_stack[row_min:row_max, col_min:col_max, :]
    end
    if texas_red_stack != nothing
        texas_red_stack = texas_red_stack[row_min:row_max, col_min:col_max, :]
    end
    if CY5_stack != nothing
        CY5_stack = CY5_stack[row_min:row_max, col_min:col_max, :]
    end
    return CFP_stack, YFP_stack, texas_red_stack, CY5_stack, img_stack, processed_stack
end

function recompute_mask(img, img_copy, raw_img, sig, blockDiameter, fpMean, fixed_thresh)
    img_normalized = normalize_local_contrast(img, img_copy, 
                                              blockDiameter, fpMean)
    img = imfilter(img_normalized, Kernel.gaussian(sig))
    plank_mask = img .> fixed_thresh 
    plank = plank_mask .* raw_img 
    flattened_plank = vec(plank)
    plank_pixels = filter(x -> x != 0, flattened_plank)
    plank_avg = mean(plank_pixels)
    threshold = thresh(plank_avg, fixed_thresh)
    mask = img .<= threshold
    return mask 
end

function compute_mask!(stack, raw_stack, masks, 
                        sig, fixed_thresh, ntimepoints, blockDiameters, fpMean)
    @inbounds for t in 1:ntimepoints
        @views img = stack[:,:,t]
        plank_mask = img .> fixed_thresh 
        @views plank = plank_mask .* raw_stack[:,:,t] 
        flattened_plank = vec(plank)
        plank_pixels = filter(x -> x != 0, flattened_plank)
        plank_avg = mean(plank_pixels)
        threshold = thresh(plank_avg, fixed_thresh)
        mask = img .<= threshold
        if plank_avg > 0.45 
            masks[:,:,t] = mask
        else
            @views new_mask = (mask .| masks[:,:,t-1]) .- masks[:,:,t-1]
            ccs = label_components(new_mask)
            lengths = component_lengths(ccs)
            largest_cc = compmax(lengths)
            @views prev_ccs = label_components(masks[:,:,t-1])
            prev_lengths = component_lengths(prev_ccs)
            largest_cluster = compmax(prev_lengths)
            if largest_cc > largest_cluster && largest_cluster < 5000 
                img = raw_stack[:,:,t]
                img_copy = raw_stack[:,:,t] 
                mask = recompute_mask(img, img_copy, img, sig, blockDiameters[2], 
                                      fpMean, fixed_thresh)
                masks[:,:,t] = mask 
                @views new_mask = (mask .| masks[:,:,t-1]) .- masks[:,:,t-1]
                ccs = label_components(new_mask)
                lengths = component_lengths(ccs)
                largest_cc = compmax(lengths)
            else
                masks[:,:,t] = mask 
            end
            if largest_cc < 1000 && largest_cluster < 100
                masks[:,:,t] = zeros(Bool, size(img))
            end
        end
    end
end

function output_images!(stack, masks, overlay, dir, well)
    @inbounds for i in CartesianIndices(stack)
        gray_val = RGB{N0f8}(stack[i], stack[i], stack[i])
        overlay[i] = masks[i] ? RGB{N0f8}(1, 0, 0) : gray_val
    end
    save("$dir/results_images/$well"*"mask.tif", overlay)
end

function main()
    config = parsefile("experiment_config.json")
    images_directories  = config["images_directory"]
    all_conditions = config["conditions"]
    sig = config["sig"]
    blockDiameter = config["blockDiameter"] 
    shift_thresh = config["shift_thresh"]
    dust_correction = config["dust_correction"]

    @inbounds for k in eachindex(images_directories)
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
        BF_output_file = "$dir/results_data/BF_imaging.csv"
        num_wells = sum(length(v) for v in values(conditions))
        files = [f for f in readdir(dir) if occursin(r"\.tif$", f) && occursin("Bright Field", f)]
        ntimepoints = div(length(files), num_wells)
        file1 = files[1]
        test_image = load("$dir/$file1"; lazyio=true)
        height, width = size(test_image)
        BF_data_matrix = Array{Float64, 2}(undef, ntimepoints, num_wells)
        CFP_data_matrix = Array{Float64, 2}(undef, ntimepoints, num_wells)
        YFP_data_matrix = Array{Float64, 2}(undef, ntimepoints, num_wells)
        texas_red_data_matrix = Array{Float64, 2}(undef, ntimepoints, num_wells)
        CY5_data_matrix = Array{Float64, 2}(undef, ntimepoints, num_wells)
        all_wells = vcat(values(conditions)...)
        conditions = values(conditions)
        CFP_flag = false
        YFP_flag = false
        texas_red_flag = false
        CY5_flag = false
        @inbounds for (i, wells) in enumerate(conditions)
            @inbounds for j in eachindex(wells)
                CFP_images = nothing
                YFP_images = nothing
                texas_red_images = nothing
                CY5_images = nothing
                images = Array{Gray{N0f16}, 3}(undef, height, width, ntimepoints)
                well = wells[j]
                BF_well_files = sort([f for f in readdir(dir) if occursin(well*"_", f) && occursin("Bright Field", f)], 
                             lt=natural)
                CFP_well_files = sort([f for f in readdir(dir) if occursin(well*"_", f) && occursin("CFP", f)], 
                             lt=natural)
                YFP_well_files = sort([f for f in readdir(dir) if occursin(well*"_", f) && occursin("YFP", f)], 
                             lt=natural)
                texas_red_well_files = sort([f for f in readdir(dir) if occursin(well*"_", f) && occursin("Texas Red", f)], 
                             lt=natural)
                CY5_well_files = sort([f for f in readdir(dir) if occursin(well*"_", f) && occursin("CY5", f)], 
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
                fpMax = maximum(images) 
                fpMin = minimum(images) 
                fpMean = (fpMax - fpMin) / 2.0 + fpMin
                fixed_thresh = fpMean - 0.04
                CFP_images, YFP_images, texas_red_images, CY5_images, images, output_stack = stack_preprocess(images, normalized_stack, 
                                                                                                              registered_stack, blockDiameter[1], 
                                                                                                              fpMean, ntimepoints, shift_thresh, 
                                                                                                              sig, CFP_images, YFP_images, 
                                                                                                              texas_red_images, CY5_images)
                masks = zeros(Bool, size(images))
                compute_mask!(output_stack, images, masks, sig, 
                              fixed_thresh, ntimepoints, blockDiameter, fpMean)
                if dust_correction == "True"
                    dust_correct!(masks)
                end
                output_stack = Gray{N0f8}.(output_stack)
                overlay = zeros(RGB{N0f8}, size(output_stack)...)
                output_images!(output_stack, masks, overlay, dir, well)
                let images = images
                    @floop for t in 1:ntimepoints
                        @inbounds signal = @views mean((1 .- images[:,:,t]) .* masks[:,:,t])
                        @inbounds BF_data_matrix[t, i*j] = signal 
                    end
                end
                if CFP_flag
                    let CFP_images = CFP_images
                        @floop for t in 1:ntimepoints
                            @inbounds signal = @views mean(CFP_images[:,:,t] .* masks[:,:,t]) - mean(CFP_images[:,:,t] .* .!masks[:,:,t])
                            @inbounds CFP_data_matrix[t, i*j] = signal 
                        end
                    end
                end
                if YFP_flag
                    let YFP_images = YFP_images
                        @floop for t in 1:ntimepoints
                            @inbounds signal = @views mean(YFP_images[:,:,t] .* masks[:,:,t]) - mean(YFP_images[:,:,t] .* .!masks[:,:,t])
                            @inbounds YFP_data_matrix[t, i*j] = signal 
                        end
                    end
                end
                if texas_red_flag
                    let texas_red_images = texas_red_images
                        @floop for t in 1:ntimepoints
                            @inbounds signal = @views mean(texas_red_images[:,:,t] .* masks[:,:,t]) - mean(texas_red_images[:,:,t] .* .!masks[:,:,t])
                            @inbounds texas_red_data_matrix[t, i*j] = signal 
                        end
                    end
                end
                if CY5_flag
                    let CY5_images = CY5_images
                        @floop for t in 1:ntimepoints
                            @inbounds signal = @views mean(CY5_images[:,:,t] .* masks[:,:,t]) - mean(CY5_images[:,:,t] .* .!masks[:,:,t])
                            @inbounds CY5_data_matrix[t, i*j] = signal 
                        end
                    end
                end
            end # loop over wells for a condition 
        end # loop over conditions 
        df = DataFrame(BF_data_matrix, Symbol.(all_wells))
        df .= ifelse.(isnan.(df), 0, df)
        write(BF_output_file, df)
        if CFP_flag
            CFP_output_file = "$dir/results_data/CFP_imaging.csv"
            df = DataFrame(CFP_data_matrix, Symbol.(all_wells))
            df .= ifelse.(isnan.(df), 0, df)
            write(CFP_output_file, df)
        end
        if YFP_flag
            YFP_output_file = "$dir/results_data/YFP_imaging.csv"
            df = DataFrame(YFP_data_matrix, Symbol.(all_wells))
            df .= ifelse.(isnan.(df), 0, df)
            write(YFP_output_file, df)
        end
        if texas_red_flag
            texas_red_output_file = "$dir/results_data/texas_red_imaging.csv"
            df = DataFrame(texas_red_data_matrix, Symbol.(all_wells))
            df .= ifelse.(isnan.(df), 0, df)
            write(texas_red_output_file, df)
        end
        if CY5_flag
            CY5_output_file = "$dir/results_data/CY5_imaging.csv"
            df = DataFrame(CY5_data_matrix, Symbol.(all_wells))
            df .= ifelse.(isnan.(df), 0, df)
            write(CY5_output_file, df)
        end
    end # loop over directories
end
main()
