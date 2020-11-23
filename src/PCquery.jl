module PCquery

import HTTP.IOExtras
using HTTP
using Mustache
using TikzGraphs,TikzPictures
using JLD2
using LightXML,DataFrames,Query
using LightGraphs, MetaGraphs
using Base.Threads


# import/export
export exportLRT

# graphs
export initGraph,filterVertices!,filterEdges!,filterVertices,filterEdges,initRow

# pathway related
export getTransTargs,getPathways,annotatePathway!

# search
export delimitValues

# NextProt_module
export getNextProt, annotateGraphFcn!,annotateGraphP!, annotateGraphG!

# orthodb module
export addExpression

# LRpairs_module
export getCxLR

# plotting
export plotDagLRT, plotDag

include("util.jl") # saving and loading data
include("plot.jl") # plotting
include("http.jl") # override delimiters
include("query.jl") # sparql tools for various endpoints
include("graph.jl") # find paths
include("pathway.jl") # helper functions for pathway-related sparql

include("PathwayCommons/query.jl") # will be separate modules

include("NextProt/query.jl") # will be separate modules
include("NextProt/annotate.jl") # will be separate modules

include("OrthoDB/query.jl") # will be separate modules
include("OrthoDB/annotate.jl") # will be separate modules

include("UniProt/query.jl") # queries to uniprot public endpoint

include("OMA/query.jl") # queries to OMA orthology public endpoint
end
