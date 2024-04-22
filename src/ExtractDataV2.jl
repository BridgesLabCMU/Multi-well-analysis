using JSON: parsefile
using DataFrames
using CSV

function extract_matrix(df, i, row_indices, valid_columns)
    row_index = row_indices[i]
    if i == length(row_indices)
        results_rows = passmissing(occursin).([r"Results"], df.Column1)
        results_rows = [ismissing(x) ? false : x for x in results_rows]
        results_row_indices = findall(results_rows)
        if length(results_row_indices) == 0
            df_subset = df[row_index+2:end, :]
        else
            df_subset = df[row_index+2:results_row_indices[1]-1, :]
        end
    else
        next_row_index = row_indices[i+1]
        df_subset = df[row_index+2:next_row_index-1, :]
    end
    delete!(rename!(df_subset, Symbol.(Vector(df_subset[1,:]))), 1)
    df_subset = df_subset[:, intersect(valid_columns, names(df_subset))]
    for (name, col) in pairs(eachcol(df_subset))
        if ismissing(df_subset[1, name])
            select!(df_subset, Not(name))
        end
    end
    dropmissing!(df_subset)
    return df_subset
end

function determine_filename(pattern)
    if occursin(r"600", pattern)
        return "OD600.csv"
    elseif occursin(r"Lum", pattern)
        return "lum.csv"
    elseif occursin(r"485/20,528/20", pattern)
        return "YFP.csv"
    elseif occursin(r"620/40,680/30", pattern)
        return "CY5.csv"
    else
        return nothing 
    end
end

function main()
    config = parsefile("experiment_config.json")
    bulk_data = config["bulk_data"]
    experiment_dir = config["experiment_directory"]
    images_dirs = config["images_directory"]
    
    for j in eachindex(bulk_data)
        bulk_data_dir = bulk_data[j][1:end-4] 
        if j < length(images_dirs) && isdir(images_dirs[j])
            images_dir = images_dirs[j]
            output_dir = images_dir 
        else
            if isdir(bulk_data_dir)
                rm(bulk_data_dir; recursive = true)
            end
            mkdir(bulk_data_dir)
            output_dir = bulk_data_dir
        end
        df = DataFrame(CSV.File(bulk_data, header=false))
        metadata_rows = passmissing(occursin).([r"Delay after plate movement"], df.Column2)  
        metadata_rows = [ismissing(x) ? false : x for x in metadata_rows]
        metadata_rows2 = passmissing(occursin).([r"Read Height"], df.Column2)
        metadata_rows2 = [ismissing(x) ? false : x for x in metadata_rows2]
        metadata_row_indices = findall(metadata_rows)
        metadata_row_indices2 = findall(metadata_rows2)
        if length(metadata_row_indices) == 0
            metadata_row_indices = metadata_row_indices2
        end
        metadata_df = df[1:metadata_row_indices[1], :]
        CSV.write("$output_dir/metadata.csv", metadata_df, header=false)

        well_names = [string(c, r) for c in 'A':'H', r in 1:12]
        rows = passmissing(occursin).([r"Read \d+"], df.Column1)
        rows = [ismissing(x) ? false : x for x in rows]
        row_indices = findall(rows)
        
        reads = [df[index, :Column1] for index in row_indices]
        for i in eachindex(row_indices)
            data_matrix = extract_matrix(df, i, row_indices, well_names)
            filename = determine_filename(reads[i])
            if filename != nothing
                CSV.write("$output_dir/$filename", data_matrix)
            end
        end

        if length(row_indices) == 0 && any(!all(ismissing, row) for row in eachrow(df[metadata_row_indices[1]+1:end, :]))  
            row_indices = metadata_row_indices .+ 2
            data_matrix = extract_matrix(df, 1, row_indices, well_names)
            filename = determine_filename(df[row_indices[1], :Column1])
            if filename != nothing
                CSV.write("$output_dir/$filename", data_matrix)
            end
        end
    end
end

main()
