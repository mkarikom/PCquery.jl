function getNested(dbParams::Dict,entity::String)
    val = delimitValues([entity],"","<>")
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/rq")
    str = open(f->read(f, String), "$rqDir/getNestedFull.rq");
    turtle = Mustache.render(str,
            Dict{Any,Any}("ent"=>val))

    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
              "Accept" => "application/sparql-results+xml"]
    fmt = "application/sparql-results+xml"
    resp = PCquery.requestTTL(dbParams[:port],"http",fmt,
                        dbParams[:host],dbParams[:path],:POST,
                        header,turtle)
    members = parseSparqlResponse(resp)
    reflist = members
    @assert size(reflist,1) > 0 "require at least one simple ref"
    reflist
end

function hasCol(df::DataFrame,testKey::Symbol)
    cx = true
    if isnothing(findfirst(n->n==string(testKey),names(df)))
        cx = false
    end
    cx
end

function getXref(dbParams::Dict,xrefs::Vector)
    colnames = dbParams[:fields]
    coltypes = fill(Union{Missing,String},length(colnames))
    df = DataFrame(coltypes,colnames)
    vals = delimitValues(xrefs,"","<>")
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/rq")
    str = open(f->read(f, String), "$rqDir/getXref.rq");
    turtle = Mustache.render(str,
            Dict{Any,Any}("xref"=>vals))

    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
              "Accept" => "application/sparql-results+xml"]
    fmt = "application/sparql-results+xml"
    resp = PCquery.requestTTL(dbParams[:port],"http",fmt,
                        dbParams[:host],dbParams[:path],:POST,
                        header,turtle)
    parseSparqlResponse!(resp,df)
    df
end

function sortUnique(cols...)
    allvert = []
    for c in cols
        cverts = collect(skipmissing(c))
        push!(allvert,cverts...)
    end
    sort(unique(allvert))
end

function deConcat!(df::DataFrame,expandCols...)
    for c in expandCols
        df[!,c] = split.(df[!,c],",")
    end
end

# find a tree of complexes, proteins, dna, and small molecules
function recurseNestedComplexes(dbParams::Dict,
                       lRef::String,rRef::String,ctrlEntityRef::String,
                       intRef::String)

end
