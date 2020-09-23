# VermontTerrestrialPassageTool.jl

Compute and summarize connectivity around road culverts in Vermont. This package was built under an assumption that spatial data will have projections with units of meters. **All input data should be in the same projection to guarantee functionality**. You will need Julia 1.5 or higher to use this package.

## Installation
```julia
julia> using Pkg; Pkg.add(PackageSpec(url = "https://github.com/csp-inc/VermontTerrestrialPassageTool.jl"))
```

## Quickstart Guide

This package provides one user-facing function, `compute_all_culverts()`. It has one parameter (argument) that specifies the file path for an INI file that contains run options.

### The INI file
The INI file should have an entry for each of the following:

- `culvert_shapefile_path` - The file path for the shapefile containing culverts.
- `culvert_id_field` - The field the culverts shapefile that contains unique IDs for each structure
- `culvert_ids_path` - The file path to a one-column CSV file (with no column names/headers) with entries specifying the IDs of the culverts for which the analysis should be run (e.g. unique IDs corresponding to the values in `culvert_id_field` in the culverts shapefile.
- `landcover_path` - The file path to the landcover raster to be used for generating resistance surfaces and masking source strength surfaces
- `reclass_table_path` - The file path to a reclass table to be used in Omniscape.jl for converting landcover into resistance. See `reclass_table` in the Omniscape.jl docs under [Arguments](https://docs.circuitscape.org/Omniscape.jl/stable/usage/#Arguments) for more info.
- `source_strength_path` - The file path to the source strength raster to be used for generating a source strength surface. This layer will be resampled and if needed, reprojected, to match the projection and resolution of the landcover dataset.
- `mask_values` - A list of values corresponding to the landcover data set to use for masking source strength. e.g., `[1 2]` specifies that any source strength pixel where the corresponding land cover pixel is not equal to 1 _or_ 2 will be set to 0.
-`os_radius_pixels` - This argument is passed as the `radius` argument to Omniscape.jl. See `radius` in the Omniscape.jl docs under [Arguments](https://docs.circuitscape.org/Omniscape.jl/stable/usage/#Arguments) for more info.
-`os_source_threshold` - This argument is passed as the `source_threshold` argument to Omniscape.jl. See `source_threshold` in the Omniscape.jl docs under [Arguments](https://docs.circuitscape.org/Omniscape.jl/stable/usage/#Arguments) for more info.
-`os_block_size` - This argument is passed as the `block_size` argument to Omniscape.jl. See `block_size` in the Omniscape.jl docs under [Arguments](https://docs.circuitscape.org/Omniscape.jl/stable/usage/#Arguments) for more info.
- `summary_radius_pixels` - The radius **in pixels** to use for summarizing connectivity around culverts. For example, `summary_radius_pixels` of 100 will result in summaries being reported for current within 100 pixels of the structure.
- `output_folder` - The file path specifying the path for a new folder where outputs should be written. Inside this folder will be a file called results.summary.csv containing summaries for current flow around each culvert. `compute_all_culverts()` will also create a new folder inside of `output_folder` for each culvert and save current maps there.

Note that all file paths should be relative to the working directory of your Julia session.

### Example
#### Example INI file:
```ini
[Culvert Options]
culvert_id_field = TNC_Join_I
culvert_shapefile_path = ../data/source/culverts/vt_culverts_bridges_statefed_3ft_AI091720_vt_utm.shp
culvert_ids_path = test_culvert_ids.csv

[Raster inputs]
landcover_path = ../data/source/LandLandcov_BaseLC2016/LandLandcov_BaseLC2016.tif
source_strength_path = ../data/source/marten/MartenOccVT30.tif

[Omniscape options]
os_radius_pixels = 200
os_source_threshold = 0.85
os_block_size = 31
reclass_table_path = reclass_table.txt

[Run options]
mask_values = [1]
summary_radius_pixels = 100

[Output options]
output_folder = output/testrun

```
#### Using the package
Once you have installed this package and have an INI file created, to use this package simply run:
```julia
julia> using VermontTerrestrialPassageTool
julia> compute_all_culverts("<file path to your INI file>")
```
Replace `<file path to your INI file>` with the file path to your INI file.
