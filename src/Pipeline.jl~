module Pipeline

export pipeline

using Gtk4, JSON
using NaturalSort
using CairoMakie, Makie, SwarmMakie
using TiffImages: load, save
using IntegralArrays
using IntervalSets
using Images
using ImageView
using StatsBase
using DataFrames
using CSV: write
using CoordinateTransformations
using AbstractFFTs
using Compat
using FFTW
using CSV
using PlotlyJS
using ImageTracking
using StaticArrays
using LinearAlgebra
using ProgressMeter

include("MainGUI.jl")  
include("ConditionsGUI.jl")
include("PlotGUI.jl")     
include("Analysis.jl")   
include("Plotting.jl")  
include("OpticalFlow.jl")  
include("ExtractMetadata.jl")

function pipeline()
    if !isinteractive()
        @async Gtk4.GLib.start_main_loop()
    end
    config_file = "experiment_config.json"

    # 1. Experiment configuration
    if isfile(config_file)
        if ask_dialog("Configuration file exists: overwrite?", nothing)
            GUI_main()
        else
            @info "Skipping GUI_main; using existing configuration."
        end
    else
        GUI_main()
    end

    if isfile(config_file)
        cfg = JSON.parsefile(config_file; use_mmap=false)
        if haskey(cfg, "bulk_data")
            ExtractMeta_main()
        else
            error("Cannot extract metadata: bulk file path not supplied to the configuration GUI")
        end
    else
        error("experiment_config.json not found: cannot proceed without running the configuration GUI")
    end

    # 2. Conditions GUI
    if ask_dialog("Do you want to run the Conditions GUI now?", nothing)
        if isfile(config_file)
            cfg = JSON.parsefile(config_file; use_mmap=false)
            if haskey(cfg, "conditions")
                if ask_dialog("Conditions already defined; overwrite?", nothing)
                    run_conditions_gui(config_file)
                else
                    @info "Skipping conditions GUI; keeping existing conditions."
                end
            else
                run_conditions_gui(config_file)
            end
        else
            error("Cannot run Conditions GUI: $config_file not found.")
        end
    end

    # 3. Plot Settings GUI
    plot_file = "plot_options.json"
    if ask_dialog("Do you want to run the Plot Settings GUI now?", nothing)
        if isfile(plot_file)
            if ask_dialog("Plot options file exists: overwrite?", nothing)
                run_plot_gui(plot_file)
            else
                @info "Skipping Plot GUI; using existing options."
            end
        else
            run_plot_gui(plot_file)
        end
    end

    # 4a. Brightfield biomass analysis 
	data_dir = ""
    if isfile(config_file)
		cfg = JSON.parsefile(config_file; use_mmap=false)
        first_dir = cfg["images_directory"][1]
        data_dir = joinpath(first_dir, "Numerical data")
    else
        error("Cannot run analysis: $config_file not found.")
    end
    if isdir(data_dir)
        if ask_dialog("Analysis results found; re-run analysis?", nothing)
            analysis_main()
        else
            @info "Skipping analysis; using previous results."
        end
    else
        analysis_main()
    end

    # 4b. Optical Flow
	if ask_dialog("Do you want to perform optical flow on your data?", nothing)
        optflow_main()
	end

    # 5. Plotting
    plotting_main()

    # 6. Move files
	if ask_dialog("Do you want to move your files out of the Multi-well-analysis folder? Only select no if you will return to these files quickly", nothing)
        if isfile(config_file)
            cfg = JSON.parsefile(config_file, use_mmap=false)
            first_dir = cfg["images_directory"][1]
            plot_dir = joinpath(first_dir, "Plots")
            if !isdir(plot_dir)
                mkpath(plot_dir)
            end
            for file in readdir(pwd())
                if occursin(".html", file) || occursin(".pdf", file)
                    mv(joinpath(pwd(), file), joinpath(plot_dir, file))
                elseif occursin(".csv", file) || occursin(".json", file)
                    mv(joinpath(pwd(), file), joinpath(first_dir, file))
                end
            end
        else
            error("Cannot move files: $config_file not found.")
        end
    end
end

end # module
