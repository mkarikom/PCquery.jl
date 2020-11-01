function getOMA(dbParams::Dict)
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/OMA/rq")
    str = open(f->read(f, String), string(rqDir,"/","getOMA.rq"));
	# turtle = str
    turtle = Mustache.render(str,
            Dict{Any,Any}("upid"=>dbParams[:upid]))
    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
              "Accept" => "application/sparql-results+xml"]
    fmt = "application/sparql-results+xml"
    resp = PCquery.requestTTL(dbParams[:port],dbParams[:protocol],fmt,
                        dbParams[:host],dbParams[:path],dbParams[:method],
                        header,turtle)

    ann = parseSparqlResponse(resp)
end

function selectOMA(oma_ids,dbParams::Dict)
	entFilter = delimitValues(oma_ids,"",["('","')"])
	dbParams[:upid] = entFilter
	unipr = getOMA(dbParams)
end
