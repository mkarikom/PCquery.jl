function getUniProt(dbParams::Dict)
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/rq")
    str = open(f->read(f, String), "$rqDir/getUniProt.rq");
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

function selectUniProt(uniprot_ids,upParams::Dict)
	entFilter = delimitValues(uniprot_ids,"",["('","')"])
	upParams[:upid] = entFilter
	unipr = PCquery.getUniProt(upParams)
end
