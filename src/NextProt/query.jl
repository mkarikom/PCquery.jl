function getNextProt(fcnParams::Dict)
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/NextProt/rq")
    str = open(f->read(f, String), string(rqDir,"/","getNextProt_gdb.rq"));
    turtle = Mustache.render(str,
            Dict{Any,Any}("goFilter"=>fcnParams[:goFilter],
						  "entityFilter"=>fcnParams[:entFilter]))
    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
              "Accept" => "application/sparql-results+xml"]
    fmt = "application/sparql-results+xml"
    resp = PCquery.requestTTL(fcnParams[:port],fcnParams[:protocol],fmt,
                        fcnParams[:host],fcnParams[:path],fcnParams[:method],
                        header,turtle)

    ann = parseSparqlResponse(resp)
end

function getNextProtGene(fcnParams::Dict)
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/NextProt/rq")
    str = open(f->read(f, String), string(rqDir,"/","getNextProtGene_gdb.rq"));
    turtle = Mustache.render(str,
            Dict{Any,Any}("unipage"=>fcnParams[:unipage]))
    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
              "Accept" => "application/sparql-results+xml"]
    fmt = "application/sparql-results+xml"
    resp = PCquery.requestTTL(fcnParams[:port],fcnParams[:protocol],fmt,
                        fcnParams[:host],fcnParams[:path],fcnParams[:method],
                        header,turtle)

    ann = parseSparqlResponse(resp)
end
