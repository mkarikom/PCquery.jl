function getNested(df::DataFrame,dbParams::Dict,ind::Int,
                    simpleF::Function,nestedF::Function,
                    refTypes::Symbol,refs::Symbol)
    nestedRefInd = findall(nestedF,df[!,refTypes][ind])
    simpleRefInd = findall(simpleF,df[!,refTypes][ind])
    refList = []
    if length(nestedRefInd) == length(df[!,refTypes][ind])
        ref = df[!,refs][ind][nestedRefInd]
        vals = delimitValues(ref,"","<>")
        srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
        rqDir = string(srcDir,"/rq")
        str = open(f->read(f, String), "$rqDir/getNested.rq");
        turtle = Mustache.render(str,
                Dict{Any,Any}("entity"=>vals))

        # compose the header and execute query
        header = ["Content-Type" => "application/x-www-form-urlencoded",
                  "Accept" => "application/sparql-results+xml"]
        fmt = "application/sparql-results+xml"
        resp = PCquery.requestTTL(7200,"http",fmt,
                            dbParams[:host],dbParams[:path],:POST,
                            header,turtle)
        members = parseSparqlResponse(resp)
        push!(refList,collect(skipmissing(members[!,dbParams[:member]]))...)
        @assert all(map(simpleF,refList))
    else
        @assert length(simpleRefInd)==1 "only 1 simple ref is allowed"
        push!(refList,df[!,refs][ind][simpleRefInd]...)
    end
    refList
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
