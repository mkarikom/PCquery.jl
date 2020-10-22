module PCquery

import HTTP.IOExtras
using HTTP
using LightXML,DataFrames,Query
using LightGraphs, MetaGraphs, GraphPlot
using Mustache

# graphs
export initBpGraph,initBpGraphCx,filterVertices!,filterEdges!,filterVertices,filterEdges

# pathway related
export getTransTargs,getPathways

# search
export delimitValues

# NextProt_module
export getNextProt, annotateGraph!

# LRpairs_module
export getCxLR

include("Http.jl") # override delimiters
include("Search.jl") # run http request against pc
include("Sparql.jl") # sparql tools for various endpoints
include("Graph.jl") # find paths
include("Pathways.jl") # helper functions for pathway-related sparql
include("Nested.jl") # testing functions
include("NextProt_module.jl") # will be separate modules
include("LRpairs_module.jl") # will be separate modules
include("OrthoDB_module.jl") # will be separate modules

end
