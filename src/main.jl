"""
    compute_all_culverts(path::String)

Parameters:

**`path`**: The file path to an INI file containing run parameters
"""
function compute_all_culverts(path::String)
    ## Parse INI Arguments
    cfg = Omniscape.parse_cfg(path)

    culvert_id_field = Symbol(cfg["culvert_id_field"])
    culvert_shapefile_path = cfg["culvert_shapefile_path"]
    culvert_ids = readdlm(cfg["culvert_ids_path"], ',', String)
    landcover_path = cfg["landcover_path"]
    os_radius_pixels = parse(Int, cfg["os_radius_pixels"])
    summary_radius_pixels = parse(Int, cfg["summary_radius_pixels"])
    os_source_threshold = parse(Float64, cfg["os_source_threshold"])
    reclass_table_path = cfg["reclass_table_path"]
    source_strength_path = cfg["source_strength_path"]
    os_block_size = parse(Int, cfg["os_block_size"])
    mask_values = float.(eval(Meta.parse(cfg["mask_values"])))
    output_folder = cfg["output_folder"]

    ## Side length derivation - assumes 0.5m resolution (hard coded)
    input_resolution = 0.5
    side_length_meters = (os_radius_pixels + summary_radius_pixels) * input_resolution * 2

    ## Setup output file structure and write empty file to store metadata
    # Create directory
    dir_suffix = 2
    while isdir(string(output_folder, "_$(dir_suffix)"))
        dir_suffix+=1
    end
    if isdir(output_folder)
        @warn string("A folder already exists at $(output_folder).",
                     " Creating a new folder",
                     ", $(string(output_folder, "_$(dir_suffix)")), ",
                     "and saving outputs there.")
        output_folder = string(output_folder, "_$(dir_suffix)")
    end
    mkpath(output_folder)

    ## Check that all ids in culvert_ids are in the shapefile
    culvert_shapefile_ids = string.(getproperty(Shapefile.Table(culvert_shapefile_path), culvert_id_field))
    bad_ids = culvert_ids[(!).(map(x -> in(x, culvert_shapefile_ids), culvert_ids))]
    if !isempty(bad_ids)
        @warn("The following user-provided culvert IDs are missing from the culverts shapefile and will be skipped. $(bad_ids)")

        missing_str_file = string(output_folder, "/missing_structures.txt")
        touch(missing_str_file)

        open(missing_str_file, "a") do f
          writedlm(f, bad_ids, ",")
        end

        culvert_ids = culvert_ids[map(x -> in(x, culvert_shapefile_ids), culvert_ids)]
    end

    # Create metadata CSV
    results_file = string(output_folder, "/results_summary.csv")
    touch(results_file)

    # write column headers to CSV
    summary_radius_str = string(Int(round(summary_radius_pixels * input_resolution)), "m")

    column_names = hcat([String(culvert_id_field)],
                        ["x_utm"],
                        ["y_utm"],
                        [string("mean_current_within_", summary_radius_str)],
                        [string("stdev_current_within_", summary_radius_str)],
                        [string("min_current_within_", summary_radius_str)],
                        [string("max_current_within_", summary_radius_str)])

    open(results_file, "w") do io
        writedlm(io, column_names, ",")
    end

    ## Build Dict for Omnsicape ini params
    os_args = Dict{String, String}()
    os_args["radius"] = string(os_radius_pixels)
    os_args["block_size"] = string(os_block_size)
    os_args["source_threshold"] = string(os_source_threshold)
    os_args["source_from_resistance"] = "false"
    os_args["reclassify_resistance"] = "true"
    os_args["solver"] = cfg["solver"]
    reclass_table = Omniscape.read_reclass_table(reclass_table_path, Float64)

    total_solves = length(culvert_ids)
    for i in 1:total_solves
        @info "Solving culvert $(culvert_ids[i]), $(i) of $(total_solves)"
        try
            os_args["project_name"] = string(output_folder, "/", culvert_ids[i])
            local output = compute_connectivity(culvert_shapefile_path,
                                                culvert_id_field,
                                                String(culvert_ids[i]),
                                                landcover_path,
                                                side_length_meters,
                                                source_strength_path,
                                                mask_values,
                                                os_args,
                                                reclass_table)

            cum_current, geotransform, x_coord, y_coord = output
            # Clip current map for summary, should always be an odd number of rows and
            # columns based on internal code
            culvert_row_y = Int(ceil((geotransform[4] - y_coord) / input_resolution)) + 1
            culvert_col_x = Int(floor((x_coord - geotransform[1]) / input_resolution)) + 1
            current_clipped = Omniscape.clip(cum_current, x = culvert_col_x,
                                             y = culvert_row_y,
                                             distance = summary_radius_pixels)

            current_mean = mean(skipmissing(current_clipped))
            current_sd = sqrt(var(skipmissing(current_clipped)))
            current_min = minimum(skipmissing(current_clipped))
            current_max = maximum(skipmissing(current_clipped))


            row_values = hcat([String(culvert_ids[i])],
                              [x_coord],
                              [y_coord],
                              [current_mean],
                              [current_sd],
                              [current_min],
                              [current_max])

            open(results_file, "a") do io
              writedlm(io, row_values, ",")
            end
        catch e
            error_message = sprint(showerror, e, catch_backtrace())
            @warn "Solve failed for culvert with ID $(culvert_ids[i]). Skipping."

            # Keep track of failed structures and their error messages
            if !isfile(string(output_folder, "/failed_structures.csv"))
                touch(string(output_folder, "/failed_structures.csv"))
                open(string(output_folder, "/failed_structures.csv"), "a") do f
                  writedlm(f, ["Culvert ID" "Error Message"], ",")
                end
            end

            open(string(output_folder, "/failed_structures.csv"), "a") do f
              writedlm(f, [culvert_ids[i] String(error_message)], ",")
            end
        end
        println(" ")
    end

    println("Done! Find all the output in $(string(pwd(), "/", output_folder))")

    nothing
end
