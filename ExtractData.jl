using JSON: parsefile
using DataFrames
using CSV

function extract_matrix(df, i, row_indices, valid_columns)
    row_index = row_indices[i] 
    if i == length(row_indices)
        df_subset = df[row_index+2:end, :]
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
    if occursin(r"Read \d+:600", pattern)
        return "OD600.csv"
    elseif occursin(r"Read \d+:Lum", pattern)
        return "lum.csv"
    elseif occursin(r"Read \d+:485/20,528/20", pattern)
        return "YFP.csv"
    elseif occursin(r"Read \d+:620/40,680/30", pattern)
        return "CY5.csv"
    else
        return nothing 
    end
end

function main()
    config = parsefile("experiment_config.json")
    bulk_data  = config["bulk_data"]
    images_dir  = config["directory"]
    
    bulk_data_dir = dirname(bulk_data)
    if isdir(images_dir)
        if isdir("$images_dir/Bulk data")
            rm("$images_dir/Bulk data"; recursive = true)
        end
        mkdir("$images_dir/Bulk data")
        output_dir = "$images_dir/Bulk data" 
    else
        if isdir("$bulk_data_dir/Bulk data")
            rm("$bulk_data_dir/Bulk data"; recursive = true)
        end
        mkdir("$bulk_data_dir/Bulk data")
        output_dir = "$bulk_data_dir/Bulk data"
    end

    df = DataFrame(CSV.File(bulk_data, header=false))
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
end

main()
