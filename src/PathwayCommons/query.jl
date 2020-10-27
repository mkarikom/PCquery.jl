function getPathways(dbParams::Dict)
    val = delimitValues(dbParams[:refs],"")
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/PathwayCommons/rq")
    str = open(f->read(f, String), string(rqDir,"/","getPC_gdb.rq"));
    turtle = Mustache.render(str,Dict{Any,Any}("pathVals"=>val))

    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
    		  "Accept" => "application/sparql-results+xml"]
    fmt = "application/sparql-results+xml"
    resp = PCquery.requestTTL(dbParams[:port],"http",fmt,
                        dbParams[:host],dbParams[:path],:POST,
                        header,turtle)
    df_pairs = PCquery.parseSparqlResponse(resp)
end

function getParticipant(dbParams::Dict,entity::String)
    val = delimitValues([entity],"","<>")
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/PathwayCommons/rq")
    str = open(f->read(f, String), string(rqDir,"/",dbParams[:rqFile]));
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
end

function getComplex(dbParams::Dict,cxref::String)
    val = delimitValues([cxref],"","<>")
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/PathwayCommons/rq")
	str = open(f->read(f, String), string(rqDir,"/","getComplex_gdb.rq"));
    turtle = Mustache.render(str,Dict{Any,Any}("cx"=>val))

    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
    		  "Accept" => "application/sparql-results+xml"]
    fmt = "application/sparql-results+xml"
    resp = PCquery.requestTTL(dbParams[:port],dbParams[:protocol],fmt,
                        dbParams[:host],dbParams[:path],dbParams[:method],
                        header,turtle)
    cxEnts = PCquery.parseSparqlResponse(resp)
end
