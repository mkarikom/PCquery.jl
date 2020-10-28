module PCquery

import HTTP.IOExtras
using HTTP
using Mustache
using LightXML,DataFrames,Query
using LightGraphs, MetaGraphs, GraphPlot
using Base.Threads

# graphs
export initGraph,filterVertices!,filterEdges!,filterVertices,filterEdges

# pathway related
export getTransTargs,getPathways

# search
export delimitValues

# NextProt_module
export getNextProt, annotateGraph!

# LRpairs_module
export getCxLR

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
