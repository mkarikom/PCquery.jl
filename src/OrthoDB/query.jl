function getOrthoDB(dbParams::Dict)
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/OrthoDB/rq")
    str = open(f->read(f, String), string(rqDir,"/","getOrthoDB.rq"));
	turtle = str
    turtle = Mustache.render(str,
            Dict{Any,Any}("entid"=>dbParams[:entid],
						  "entpfx"=>dbParams[:entpfx],
			              "gname"=>dbParams[:gname],
						  "spname"=>dbParams[:spname]))
    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
              "Accept" => "application/sparql-results+xml"]
    fmt = "application/sparql-results+xml"
    resp = PCquery.requestTTL(dbParams[:port],dbParams[:protocol],fmt,
                        dbParams[:host],dbParams[:path],dbParams[:method],
                        header,turtle)

    ann = parseSparqlResponse(resp)
end

function preGetOrthoDB!(q,origids,dbParams,df)
	firstind = findfirst(x->x==q[1],origids)
	lastind = findfirst(x->x==q[end],origids)
	entFilter = delimitValues(q,dbParams[:entpfxs],"")
	dbParams[:entid] = entFilter
	println("retrieving ids $firstind : $lastind")
	append!(df,getOrthoDB(dbParams))
end

# serially process a list of uniprot ids
function selectOrthologs!(endids,dbParams::Dict,df::DataFrame)
	origids = copy(endids)
	q = []
	while length(endids) > 0
		if length(q) < dbParams[:maxq]
			id = popat!(endids,1,missing)
			if !ismissing(id)
				push!(q,id)
			else
				preGetOrthoDB!(q,origids,dbParams,df)
				q = []
			end
		else
			preGetOrthoDB!(q,origids,dbParams,df)
			q = []
		end
	end
	if length(q) > 0 # clear the queue
		preGetOrthoDB!(q,origids,dbParams,df)
	end
	df
end
