module PCquery

import HTTP.IOExtras
using HTTP
using LightXML,DataFrames,Query
using LightGraphs, MetaGraphs, GraphPlot
using Mustache

# graphs
export getPathways,initBpGraph,filterVertices!,filterEdges!

# pathway related
export getTransTargs

# search
export getNextProt, delimitValues

include("Http.jl") # override delimiters
include("Search.jl") # run http request against pc
include("Sparql.jl") # sparql tools for various endpoints
include("Graph.jl") # find paths
include("Pathways.jl") # helper functions for pathway-related sparql
#include("nestedOld.jl") # testing functions
include("Nested.jl") # testing functions
end
