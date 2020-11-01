## Module: test UniProt SPARQL query
# upParams = Dict(
# 	:host=>"sparql.uniprot.org",
# 	:port=>443,
# 	:protocol=>"https",
# 	:method=>:POST,
# 	:path=>"/")
#
# simpleVerts = filterVertices(g,:entIdDb,p->occursin("uniprot",p))
# uniprot_ids = map(v->props(g,v)[:entId],simpleVerts)
# uniprot = PCquery.exploreUniProt(uniprot_ids[1],upParams)


function getUniProt(dbParams::Dict)
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/UniProt/rq")
    str = open(f->read(f, String), string(rqDir,"/","getUniProt.rq"));
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

function exploreUniProt(uniprot_ids,upParams::Dict)
	entFilter = delimitValues(uniprot_ids,"",["('","')"])
	upParams[:upid] = entFilter
	unipr = getUniProt(upParams)
end
