function ExtractMeta_main()
    config = JSON.parsefile("experiment_config.json")
    bulk_data = config["bulk_data"]
    images_dirs = config["images_directory"]
    output_dir = images_dirs[1] 
    df = CSV.read(bulk_data, DataFrame, header=false; delim=',', ntasks=1)
    metadata_rows = passmissing(occursin).([r"Delay after plate movement"], df.Column2)  
    metadata_rows = [ismissing(x) ? false : x for x in metadata_rows]
    metadata_rows2 = passmissing(occursin).([r"Read Height"], df.Column2)
    metadata_rows2 = [ismissing(x) ? false : x for x in metadata_rows2]
    metadata_row_indices = findall(metadata_rows)
    metadata_row_indices2 = findall(metadata_rows2)
    if length(metadata_row_indices) == 0
        metadata_row_indices = metadata_row_indices2
    end
    metadata_df = df[1:metadata_row_indices[end], :]
    CSV.write("$output_dir/metadata.csv", metadata_df, header=false)

    # Extract protocol 
	steps = DataFrame(
		step          = Int[],
		action        = String[],
		magnification = Union{String,Missing}[],
		channel       = Union{String,Missing}[],
	)

	step_counter = 0

	for i in 1:nrow(df)
		c1 = df[i, 1]
		c2 = df[i, 2]

		if !ismissing(c1) && c1 == "Set Temperature"
			step_counter += 1
			push!(steps, (step_counter, "Set Temperature", missing, missing))

		elseif !ismissing(c1) && c1 == "Read"
			step_counter += 1

			if !ismissing(c2) && occursin(r"^\d+x$", c2)
				push!(steps, (step_counter, "Imaging Read", c2, missing))

			elseif !ismissing(c2) && c2 == "Image Single Image"
				objective_val = missing
				if i + 2 <= nrow(df)
					obj_cell = df[i+2, 2]
					if !ismissing(obj_cell)
						m = match(r"^Objective:\s*(.*)", obj_cell)
						objective_val = m === nothing ? missing : m.captures[1]
					end
				end
				push!(steps, (step_counter, "Imaging Read", objective_val, missing))

			else
				push!(steps, (step_counter, "Read (non-imaging)", missing, missing))
			end

		elseif !ismissing(c2) && occursin(r"^Channel \d+:\s*", c2)
			channel_name = replace(c2, r"^Channel \d+:\s*" => "")
			channel_name = replace(channel_name, r"\s+\d+$" => "")

			last_idx = findlast(==( "Imaging Read"), steps.action)
			if last_idx === nothing
				@warn "Channel found without preceding Imaging Read: $channel_name"
				continue
			end

			if ismissing(steps[last_idx, :channel])
				steps[last_idx, :channel] = channel_name
			else
				push!(steps, (
					steps[last_idx, :step],
					"Imaging Read",
					steps[last_idx, :magnification],
					channel_name
				))
			end
		end
	end
	CSV.write(joinpath(output_dir, "protocol.csv"), steps; writeheader=true)
end
