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

using LaTeXStrings
using PythonCall
using PythonPlot
using StatsBase 

round_odd(x) = div(x, 2) * 2 + 1
compmax(x) = length(x) > 1 ? maximum(x[1:end]) : 0

function read_images(directory, files)
    file = files[1]
    img = load("$directory/$file")
    return img
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

function stack_preprocess(img_stack, normalized_stack, registered_stack, blockDiameter, nframes, mxshift, sig)       
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
    height,width, _ = size(img_stack)
    
    if (height > width) 
        img_stack = permutedims(img_stack, [2, 1, 3])
        processed_stack = permutedims(processed_stack, [2, 1, 3])
        height,width, _ = size(img_stack)
        row_min, row_max, col_min, col_max = col_min, col_max, row_min, row_max
    end

    img_stack = img_stack[row_min:row_max, col_min:col_max, :]
    img_stack = img_stack[Int(round(.15 * height)):Int(round(.85 * height)), Int(round(.15 * width)):Int(round(.85 * width)), :]
    processed_stack = processed_stack[Int(round(.15 * height)):Int(round(.85 * height)), Int(round(.15 * width)):Int(round(.85 * width)), :]
    return img_stack, processed_stack
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
	normalized = normalize_local_contrast_output(normalized, stack, copy(stack), 101, fpMean)
	normalized = Gray{N0f8}.(normalized)
    save("$dir/results_images/$well.tif", normalized)
    @inbounds for i in CartesianIndices(normalized)
        gray_val = RGB{N0f8}(normalized[i], normalized[i], normalized[i])
        overlay[i] = masks[i] ? RGB{N0f8}(0,1,1) : gray_val
    end
    save("$dir/results_images/$well"*"mask.tif", overlay)
end

function plottingfunc(df, csv_output_file, output_file, acquisition_frequency)
    t = range(0, stop=nrow(df) - 1, length=nrow(df)) .÷ acquisition_frequency
    PythonPlot.matplotlib.rcParams["figure.figsize"] = [4, 3]
    PythonPlot.matplotlib.rcParams["lines.linewidth"] = 2.5
    PythonPlot.matplotlib.rcParams["axes.titlesize"] = 20
    PythonPlot.matplotlib.rcParams["axes.labelsize"] = 20
    PythonPlot.matplotlib.rcParams["xtick.labelsize"] = 14
    PythonPlot.matplotlib.rcParams["ytick.labelsize"] = 14
    PythonPlot.matplotlib.rcParams["legend.fontsize"] = 12

    final_values = [df[end, i] for i in 1:ncol(df)]
    peak_values = [maximum(df[:, i]) for i in 1:ncol(df)]
    avg_final_value = mean(final_values)
    std_final_value = std(final_values)
    avg_peak_value = mean(peak_values)
    std_peak_value = std(peak_values)

    fig, ax = pyplot.subplots()
    
    sig_wells = []

    for i in 1:ncol(df)
        final_val = final_values[i]
        peak_val = peak_values[i]
        col_name = names(df)[i]
        if final_val >= avg_peak_value * 0.1 && peak_val >= avg_peak_value * 1.5
            color = "green"
        elseif final_val >= avg_peak_value * 0.1
            color = "red"
        elseif peak_val >= avg_peak_value * 1.5
            color = "blue"
        elseif peak_val <= avg_peak_value * 0.5
            color = "purple"
        else
            color = "black"
        end
        ax.plot(t, df[:, i], color=color, alpha=0.2)
        
        if color != "black"
            normalized_peak_val = peak_val / avg_peak_value
            normalized_final_val = final_val / avg_peak_value
            push!(sig_wells, (col_name, normalized_peak_val, normalized_final_val))
        end
    end

    avg_line = reduce(+, eachcol(df)) ./ ncol(df)
    ax.plot(t, avg_line, color="black", linewidth=2, label="Average")
    ax.set_ylabel("Biofilm biomass (a.u.)")
    ax.set_xlabel("Time (h)")
    pyplot.locator_params(axis='x', min_n_ticks=5)
    pyplot.locator_params(axis='y', min_n_ticks=5)
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    pyplot.tight_layout()
    savefig(output_file)

    df_sig = DataFrame(
        ColumnName = [row[1] for row in sig_wells],
        NormalizedPeakValue = [row[2] for row in sig_wells],
        NormalizedFinalValue = [row[3] for row in sig_wells]
    )
	csv_output = replace(csv_output_file, "BF_imaging" => "significant_wells")
	CSV.write(csv_output, df_sig)
end

function multiple_mags(dir, bulk_file)
    df = DataFrame(CSV.File(bulk_file, header=false))
    objective_rows = passmissing(occursin).([r"Objective"], df.Column2)  
	objective_rows = [ismissing(x) ? false : x for x in objective_rows]
    objective_rows = findall(objective_rows)
    if length(objective_rows) > 1
        flags = ["_03_1_", "_04_1_"]
        plot_output_files = []
        BF_output_files = []
        for i in 1:length(objective_rows)
            objective = match(r"(\d{1,2}x)", df[objective_rows[i], 2]).match
            push!(plot_output_files, "$dir/results_data/plot_"*objective*".svg")
            push!(BF_output_files, "$dir/results_data/BF_imaging_"*objective*".csv")
        end
    else
        flags = ["_03_1_"]
        plot_output_files = ["$dir/results_data/plot.svg"]
        BF_output_files = ["$dir/results_data/BF_imaging.csv"]
    end
    return flags, BF_output_files, plot_output_files
end

function main()
    config = parsefile("experiment_config.json")
    images_directories  = config["images_directory"]
    acquisition_frequency = config["acquisition_frequency"]
    bulk_file = "BF_10x_Time_Biospa_MJ_Drawer2 10-Oct-2024 15-53-34BF_10x_Time_Biospa_MJ_241011_.csv"
    sig = 2
    blockDiameter = [501, 101] 
    shift_thresh = 50
    dust_correction = "True" 

    @inbounds for k in eachindex(images_directories)
        dir = images_directories[k]

        if isdir("$dir/results_images")
            rm("$dir/results_images"; recursive = true)
        end
        if isdir("$dir/results_data")
            rm("$dir/results_data"; recursive = true)
        end
        mkdir("$dir/results_images")
        mkdir("$dir/results_data")
        
        flags, BF_output_files, plot_output_files = multiple_mags(dir, bulk_file)
        for (i, flag) in enumerate(flags)
            num_wells = 96 
            letters = ["A", "B", "C", "D", "E", "F", "G", "H"]
            numbers = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
            all_wells = [x*y for x in letters for y in numbers]
            files = [f for f in readdir(dir) if occursin(r"^[A-H]\d{1,2}__raw\.tif$", f)]
            file1 = files[1]
            test_image = load("$dir/$file1")
            height, width, ntimepoints = size(test_image)
            if (height > width) 
                test_image = permutedims(test_image, [2, 1, 3])
                height, width, ntimepoints = size(test_image)
            end

            test_image = test_image[Int(round(.15 * height)):Int(round(.85 * height)), Int(round(.15 * width)):Int(round(.85 * width)), :]
            BF_data_matrix = Array{Float64, 2}(undef, ntimepoints, num_wells)
            idx = 0
            nlayers = 96 * ntimepoints
            @inbounds for j in eachindex(all_wells)
                idx += 1
                well = all_wells[j]
                BF_well_files = [f for f in readdir(dir) if occursin(well, f)]
                images = read_images(dir, BF_well_files)
                images = Float64.(images)
                normalized_stack = similar(images)
                registered_stack = similar(images)
                fixed_thresh = 0.02
                images, output_stack = stack_preprocess(images, normalized_stack, registered_stack, blockDiameter[2], ntimepoints, shift_thresh, 
                                                                                                              sig)
                masks = zeros(Bool, size(images))
                compute_mask!(output_stack, masks, fixed_thresh, ntimepoints)
                if dust_correction == "True"
                    dust_correct!(masks)
                end
                overlay = zeros(RGB{N0f8}, size(output_stack)...)
                if length(flags) > 1
					str = BF_output_files[i]
                    well = well*"_"*str[findlast(isequal('_'), str)+1:length(str)-4]
                end
                output_images!(images, masks, overlay, dir, well)
                let images = images
                    @floop for t in 1:ntimepoints
                        @inbounds signal = @views mean((1 .- images[:,:,t]) .* masks[:,:,t])
                        @inbounds BF_data_matrix[t, idx] = signal 
                    end
                end
            end # loop over wells for a condition 
            df = DataFrame(BF_data_matrix, Symbol.(all_wells))
            df .= ifelse.(isnan.(df), 0, df)
            plottingfunc(df, BF_output_files[i], plot_output_files[i], acquisition_frequency)
            CSV.write(BF_output_files[i], df)
        end # loop over magnifications
    end # loop over directories
end
main()