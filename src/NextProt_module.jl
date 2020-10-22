function getNextProt(fcnParams::Dict)
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/rq")
    str = open(f->read(f, String), "$rqDir/getNextProt_gdb.rq");
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

function annotateGraph!(fcnParams::Dict,g::AbstractMetaGraph)
	simpleVerts = filterVertices(g,:entIdDb,p->occursin("uniprot",p))
	uniprot_ids = map(v->props(g,v)[:entId],simpleVerts)
	entFilter = delimitValues(map(r->string(split(r,"/")[end]),unique(uniprot_ids)),"unipage:")
	fcnParams[:entFilter] = entFilter
	ann = getNextProt(fcnParams)

	# get the shorthand and add to the result
	terms = fcnParams[:annotationMapAncTerm](ann)
	insertcols!(ann, 3, :goAncShort=>terms)

	ctrlFcnVert = map(uniprot->collect(filter_vertices(g,:ctrlEntityEntId,uniprot)),map(iri->split(iri,"/")[end],ann.uniprot))
	fcnVert = map(uniprot->collect(filter_vertices(g,:entId,uniprot)),map(iri->split(iri,"/")[end],ann.uniprot))
	insertcols!(ann, 3, :ctrlFcnVert=>ctrlFcnVert)
	insertcols!(ann, 3, :fcnVert=>fcnVert)

	a = @from i in ann begin
		@group i by {i.goAncShort} into g
		@select {goAncShort=reduce(vcat,unique(g.goAncShort)),verts=vcat(g.fcnVert...)}
		@collect DataFrame
	end

	for r in 1:length(a[!,:goAncShort])
		for v in a[r,:verts]
			set_prop!(g, v, fcnParams[:annotationKey], a[r,:goAncShort])
		end
	end
	a
end

function getNextProtGene(fcnParams::Dict)
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/rq")
    str = open(f->read(f, String), "$rqDir/getNextProtGene_gdb.rq");
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

# add gene ids to all the proteins
function annotateGraphGene!(fcnParams::Dict,g::AbstractMetaGraph)
	prot_v = filterVertices(g,:participantType,i->i=="http://www.biopax.org/release/biopax-level3.owl#Protein",a->get(a,:entId,missing))
	fcnParams[:unipage]=delimitValues(unique(prot_v.val),"unipage:")
	ann = getNextProtGene(fcnParams::Dict)

	# get the shorthand and add to the result
	upids = map(x->split(x,"/")[end],ann[!,:uniprot])
	insertcols!(ann, 1, :upid=>upids)

	for v in 1:length(prot_v.ind)
		ind = findfirst(x->x==prot_v.val[v],ann[!,:upid])
		ens = split(ann[ind,:gene],"/")[end]
		set_props!(g, prot_v.ind[v], Dict(:ensId=>ens))
		setval = props(g,prot_v.ind[v])[:ensId]
		println("success: $setval")
	end
	ann
end
