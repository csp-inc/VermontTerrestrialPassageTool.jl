module VermontTerrestrialPassageTool

using ArchGDAL
using DelimitedFiles
using GeoData
using Omniscape
using Shapefile
using Statistics

include("utils.jl")
include("main.jl")
include("consts.jl")

export compute_all_culverts

end
