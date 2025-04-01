using ImageTracking
using Images
using CoordinateTransformations, StaticArrays
using JSON
using CSV, DataFrames
using LinearAlgebra
using TiffImages

function main()
    config = JSON.parsefile("experiment_config.json")
    images_directories  = config["images_directory"]
    for directory in images_directories
        data_path = joinpath(directory, "results_data")
        image_paths = [f for f in readdir(joinpath(directory, 
                                       "results_images"),
					  join=true) if !occursin("mask", f) && occursin("tif", f)]
        image_displacements = []
        for image_path in image_paths
            source_name = basename(image_path)
            images = Gray{Float32}.(TiffImages.load(image_path))
            height, width, frames = size(images)

            algorithm = LucasKanade(10, window_size = 20,
                                        pyramid_levels = 4,
                                        eigenvalue_threshold = 5e-5)
            avg_displacements = Array{Float64, 1}(undef, frames-1)
            area = height*width
            for t in 1:frames-1
                @views corners = imcorner(images[:,:,t],
                                          method=shi_tomasi)
                I = findall(!iszero, corners)
                r, c = (getindex.(I, 1), getindex.(I, 2))
                points = map((ri, ci) -> SVector{2}(Float64(ri),
                                                    Float64(ci)),
                             r, c)
                @views flow, indicator = optical_flow(images[:,:,t],
                                                  images[:,:,t+1], 
                                                      points,
                                                      algorithm)

                valid_flow = flow[indicator]
                if !isempty(valid_flow)
                    avg_displacements[t] = sum(sqrt.(getindex.(valid_flow, 2).^2 + getindex.(valid_flow, 1).^2)) / area
                else
                    avg_displacements[t] = 0 
                end
            end
            push!(image_displacements, avg_displacements)
        end
        df = DataFrame(Dict(Symbol(name) => data for (name, data) in zip(image_paths, image_displacements)))
        CSV.write(joinpath(data_path, "average_displacement.csv"),
                  df)
    end
end

main()
