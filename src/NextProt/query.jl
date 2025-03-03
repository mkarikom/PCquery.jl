function getNextProtFcn(fcnParams::Dict)
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/NextProt/rq")
    str = open(f->read(f, String), string(rqDir,"/","getNextProtFcn_gdb.rq"));
    turtle = Mustache.render(str,
            Dict{Any,Any}("goFilter"=>fcnParams[:goFilter],
						  "entityFilter"=>fcnParams[:entFilter],
   						  "fromgraph"=>fcnParams[:fromgraph]))
    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
              "Accept" => "application/sparql-results+xml"]
    fmt = "application/sparql-results+xml"
    resp = PCquery.requestTTL(fcnParams[:port],fcnParams[:protocol],fmt,
                        fcnParams[:host],fcnParams[:path],fcnParams[:method],
                        header,turtle)

    ann = parseSparqlResponse(resp)
end

function getNextProtP(fcnParams::Dict)
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/NextProt/rq")
    str = open(f->read(f, String), string(rqDir,"/","getNextProtP_gdb.rq"));
    turtle = Mustache.render(str,
            Dict{Any,Any}("unipage"=>fcnParams[:unipage],
						  "fromgraph"=>fcnParams[:fromgraph]))
    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
              "Accept" => "application/sparql-results+xml"]
    fmt = "application/sparql-results+xml"
    resp = PCquery.requestTTL(fcnParams[:port],fcnParams[:protocol],fmt,
                        fcnParams[:host],fcnParams[:path],fcnParams[:method],
                        header,turtle)

    ann = parseSparqlResponse(resp)
end

function getNextProtG(fcnParams::Dict)
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/NextProt/rq")
    str = open(f->read(f, String), string(rqDir,"/","getNextProtG_gdb.rq"));
    turtle = Mustache.render(str,
            Dict{Any,Any}("gene"=>fcnParams[:gene],
						  "fromgraph"=>fcnParams[:fromgraph]))
    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
              "Accept" => "application/sparql-results+xml"]
    fmt = "application/sparql-results+xml"
    resp = PCquery.requestTTL(fcnParams[:port],fcnParams[:protocol],fmt,
                        fcnParams[:host],fcnParams[:path],fcnParams[:method],
                        header,turtle)

    ann = parseSparqlResponse(resp)
end
