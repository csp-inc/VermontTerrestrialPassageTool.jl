"""
    get_culvert_coords(culvert_shapefile_path,
                       culvert_id_field),
                       culvert_id)

Get the X and Y coordinates for a culvert by ID

Parameters:
- `culvert_shapefile_path`: The file path to the shapefile containing culverts
- `culvert_id_field`: A Symbol for the field in the culvert shapefile to use for
uniquely identifying culverts. e.g. `Symbol("UID_FIELD")`
- `culvert_id`: The unique ID (of the field specified by `culvert_id_field`)
of the desired culvert

Returns a Tuple in the form of (X\\_COORDINATE, Y\\_COORDINATE)
"""
function get_culvert_coords(culvert_shapefile_path::String,
                            culvert_id_field::Symbol,
                            culvert_id::String)
    # Load the culverts shapefile and grab Table representation
    culverts_table = Shapefile.Table(culvert_shapefile_path)
    idxs = getproperty(culverts_table, culvert_id_field) .== culvert_id

    n_records = sum(idxs)
    if n_records == 0
        @error("No culvert exists with $(:culvert_id_field) = $(culvert_id)")
    elseif n_records > 1
        @warn("More than one culvert exists with $(:culvert_id_field) = $(culvert_id). Using the first matching record.")
    end

    culvert = Shapefile.shapes(culverts_table)[findfirst(idxs)]

    coord_x_column = culvert.x
    coord_y_row = culvert.y

    return coord_x_column, coord_y_row
end

"""
    extract_landcover(raster_path, side_length, x_coordinate, y_coordinate)

Extract a square from a landcover raster of side length `side_length`, centered on
the point at [y_coordinate, x_coordinate]. `side_length` will be automatically
adjusted (by no more that the width of single pixel) if needed to ensure that
the result is an odd number of pixels in width and height. This makes it play
nicer with `Omniscape.clip`.

Parameters:
- `raster_path`: The file path to the GeoTIFF raster containing landcover info.
- `side_length`: The side length of the square to extract in units of the projection of the raster.
- `x_coordinate`: X coordinate of the center of the square in units of the projection of the raster.
- `y_coordinate`: Y coordinate of the center of the square in units of the projection of the raster.

Returns a `GeoData.GeoArray`
"""
function extract_landcover(raster_path::String,
                           side_length::Number,
                           x_coordinate::Number,
                           y_coordinate::Number)
    # Load the GDALarray (GeoData.jl)
    landcover = GDALarray(raster_path)
    # Get dims to use for subsetting
    lc_dims = dims(landcover)
    lc_size = size(landcover)

    # Get dim values
    x_dim = lc_dims[1].val
    y_dim = lc_dims[2].val

    if iseven(Int(round(side_length / step(x_dim))))
        side_length = side_length + step(x_dim)
    end

    # Extract full extent of raster
    # Need to add step to get extent because entries in
    # x_dim and y_dim correspond to left and top edges of pixel, so need to
    # incorporate width of last row and column of pixels to get full extent
    west_bound = minimum(x_dim)
    east_bound = maximum(x_dim) + step(x_dim)
    south_bound = minimum(y_dim) + step(y_dim)
    north_bound = maximum(y_dim)

    # snap center point of square (that will be extracted from `landcover`)
    # to the nearest pixel center
    pixel_center_coord_x_column = maximum(x_dim[x_dim .< x_coordinate]) + step(x_dim)/2
    pixel_center_coord_y_row = maximum(y_dim[y_dim .< y_coordinate]) + step(y_dim)/2

    # Get extent boundary for subsetting the landcover raster
    # Use min and max to ensure extent isn't beyond the full raster's extent
    landcover_west_bound = max(pixel_center_coord_x_column - side_length / 2, west_bound)
    landcover_east_bound = min(pixel_center_coord_x_column +  side_length / 2, east_bound)
    landcover_south_bound = max(pixel_center_coord_y_row -  side_length / 2, south_bound)
    landcover_north_bound = min(pixel_center_coord_y_row +  side_length / 2, north_bound)

    # Add tiny padding to extent bounds and subset to get a GeoData.GeoArray
    # The padding ensures all pixels within the bounds are extracted
    landcover_ga = landcover[Lat(Between(landcover_south_bound - EXTENT_NUDGE, landcover_north_bound + EXTENT_NUDGE)),
                             Lon(Between(landcover_west_bound - EXTENT_NUDGE, landcover_east_bound + EXTENT_NUDGE))]

    # landcover_array = Array{Float64, 2}(landcover_ga.data[:,:,1])

    return permutedims(landcover_ga, (Lat, Lon, Band)), abs(step(y_dim)), pixel_center_coord_x_column, pixel_center_coord_y_row
end


"""
    prep_source_strength(source_strength_path, landcover_ga)

Parameters:
- `source_strength_path`: The file path to the layer representing source
strength.
- `landcover_ga`: The landcover GeoData.GeoArray to use as a snap raster for
resampling and masking the source strength raster.
"""
function prep_source_strength(source_strength_path::String,
                              landcover_ga::GeoArray)
    wkt = crs(landcover_ga).val
    lc_dims = dims(landcover_ga)

    # Get dim values
    # expects a north-side-up GeoArray, so 1st (rows) dim should
    # correspond to Latitude, y
    y_dim = lc_dims[1].val
    x_dim = lc_dims[2].val

    # Extract bounds
    # Need to add step to last entries to get extent because entries in
    # x_dim and y_dim correspond to left and top edge of pixel, so need to
    # incorporate width of last row and column of pixels to get full extent
    west_bound = minimum(x_dim)
    east_bound = maximum(x_dim) + step(x_dim)
    south_bound = minimum(y_dim) + step(y_dim)
    north_bound = maximum(y_dim)

    ArchGDAL.read(source_strength_path) do source
        ArchGDAL.gdalwarp([source], ["-t_srs", "$(wkt)",
                                     "-te", "$(west_bound + EXTENT_NUDGE)",
                                     "$(south_bound + EXTENT_NUDGE)",
                                     "$(east_bound + EXTENT_NUDGE)",
                                     "$(north_bound + EXTENT_NUDGE)",
                                     "-te_srs", "$(wkt)",
                                     "-tr", "$(step(x_dim))", "$(step(y_dim))",
                                     "-r", "cubic"]) do warped
            band = ArchGDAL.getband(warped, 1)

            geotransform = ArchGDAL.getgeotransform(warped)

            # Extract the array
            array_t = ArchGDAL.read(band)

            # This handles UInt tiff rasters that can still have negative NoData values
            # Need to convert the NoData value to Int64 in these cases
            if eltype(array_t) <: Integer
                ras_type = Int64
            else
                ras_type = eltype(array_t)
            end

            # Extract no data value, first converting it to the proper type (based on
            # the raster). Then, need to convert to Float64. Weird, yes,
            # but it's the only way I could get it to work for all raster types... -VL
            nodata_val = convert(Float64, convert(ras_type, ArchGDAL.getnodatavalue(band)))

            # Transpose the array -- ArchGDAL returns a x by y array, need y by x
            # Line to handle NaNs in datasets read from tifs
            array_t[isnan.(array_t)] .= nodata_val

            array = convert(Array{Union{Float64, Missing}, 2}, permutedims(array_t, [2, 1]))

            array[array .== nodata_val] .= missing

            array, geotransform
        end
    end
end


"""
Prep inputs and run Omniscape for a single culvert
"""
function compute_connectivity(culvert_shapefile_path::String,
                              culvert_id_field::Symbol,
                              culvert_id::String,
                              landcover_path::String,
                              side_length_meters::Number,
                              source_strength_path::String,
                              mask_values::Array{Float64, 1},
                              os_args::Dict{String, String},
                              reclass_table::Array{Union{Float64, Missing}, 2}
                              )
    # Get culvert coordinates
    x_coordinate, y_coordinate = get_culvert_coords(
        culvert_shapefile_path,
        culvert_id_field,
        culvert_id
    )

    # Prep landcover
    lc_output = extract_landcover(landcover_path, side_length_meters,
                                  x_coordinate, y_coordinate)

    landcover_ga, resolution, pixel_x_column, pixel_y_row =  lc_output

    ### read in source and warp
    source_strength, geotransform = prep_source_strength(source_strength_path, landcover_ga)

    # Convert landcover to Array{Union{Float64, Missing}, 2} for Omniscape
    landcover = Array{Union{Missing, Float64}, 2}(landcover_ga.data[:,:,1])
    lc_missingval = Float64(landcover_ga.missingval)
    landcover[landcover .== lc_missingval] .= missing

    # Mask source strength with landcover, and give source strength a missing val where
    # landcover is missing
    source_strength[(!).(coalesce.(map(x -> in(x, mask_values), landcover), true))] .= missing

    # Grab wkt
    wkt = crs(landcover_ga).val

    # Objects to pass to Omniscape.jl, landcover (resistance), source_strength,
    # geotransform, wkt

    # Clip inputs
    center_pixel = Int(size(landcover)[1] / 2 + 0.5)
    landcover_clipped = Omniscape.clip(landcover, x = center_pixel, y = center_pixel, distance = center_pixel - 1)
    sources_clipped = Omniscape.clip(source_strength, x = center_pixel, y = center_pixel, distance = center_pixel - 1)

    if isdir("$(os_args["project_name"])")
        rm("$(os_args["project_name"])", force = true, recursive = true)
    end

    cum_current = run_omniscape(os_args, landcover_clipped,
                                reclass_table = reclass_table,
                                source_strength = sources_clipped,
                                wkt = wkt,
                                geotransform = geotransform,
                                write_outputs = true)

    cum_current, geotransform, pixel_x_column, pixel_y_row
end

