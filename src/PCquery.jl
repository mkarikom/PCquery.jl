module PCquery

import HTTP.IOExtras
using LightXML,DataFrames


# http
export issafe,escapeuri,delimitedAS
# search
export resultData

include("Http.jl") # override delimiters
include("Search.jl") # run http request against pc
include("Sparql.jl") # sparql tools for various endpoints
end
