# identify biochemical reactions where the input is Dna and output is protein
function findTranscriptionBC(g::AbstractMetaGraph)
    ctrl = ["http://www.biopax.org/release/biopax-level3.owl#Catalysis",
 			"http://www.biopax.org/release/biopax-level3.owl#Control"]

    colnames = [:ctrlInd,:geneInd,:protInd,
                :ctrlRef,:geneRef,:protRef]
    coltypes = [Union{Missing,Int64},Union{Missing,Int64},Union{Missing,Int64},
                Union{Missing,String},Union{Missing,String},Union{Missing,String}]
    edgeTable = DataFrame(coltypes,colnames)

	for v in vertices(g)
		# identify biochemical rxn
		if props(g,v)[:entityType] == "http://www.biopax.org/release/biopax-level3.owl#BiochemicalReaction"
			ns_o = outneighbors(g,v)
	        ns_i = inneighbors(g,v)
			for n_o in ns_o
				for n_i in ns_i
					if props(g,n_i)[:entityType] == "http://www.biopax.org/release/biopax-level3.owl#Dna"
						if props(g,n_o)[:entityType] == "http://www.biopax.org/release/biopax-level3.owl#Protein"
							rxref = props(g,v)[:unification]
							inref = props(g,n_i)[:unification]
							outref = props(g,n_o)[:unification]
							data = tuple(v,n_i,n_o,inref,rxref,outref)
							push!(edgeTable,initRow(colnames,colnames,data))
						end
					end
				end
			end
		end
	end
	edgeTable
end

function getPathways(dbParams::Dict)
    val = delimitValues(dbParams[:refs],"")
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/rq")
    str = open(f->read(f, String), "$rqDir/pc_f0_path_gdb.rq");
    turtle = Mustache.render(str,Dict{Any,Any}("pathVals"=>val))

    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
    		  "Accept" => "application/sparql-results+xml"]
    fmt = "application/sparql-results+xml"
    host = "192.168.1.3"
    resp = PCquery.requestTTL(dbParams[:port],"http",fmt,
                        dbParams[:host],dbParams[:path],:POST,
                        header,turtle)
    df_pairs = PCquery.parseSparqlResponse(resp)
end

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

# query nextprot to get GO data on proteins
function getNP(fcnParams::Dict,g::AbstractMetaGraph,df::DataFrame)
	uniprot_ids = fcnParams[:dbFilter](df)
	entFilter = delimitValues(map(r->string(split(r,"/")[end]),uniprot_ids),"unipage:")
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/rq")
    str = open(f->read(f, String), "$rqDir/getNextProt_gdb.rq");
    turtle = Mustache.render(str,
            Dict{Any,Any}("goFilter"=>fcnParams[:fcnRefs],
						  "entityFilter"=>entFilter))
    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
              "Accept" => "application/sparql-results+xml"]
    fmt = "application/sparql-results+xml"
    resp = PCquery.requestTTL(fcnParams[:port],fcnParams[:protocol],fmt,
                        fcnParams[:host],fcnParams[:path],fcnParams[:method],
                        header,turtle)

    ann = parseSparqlResponse(resp)
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
end

function hasCol(df::DataFrame,testKey::Symbol)
    cx = true
    if isnothing(findfirst(n->n==string(testKey),names(df)))
        cx = false
    end
    cx
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

function traverseGraph(g::AbstractMetaGraph,sources::Vector,destinations::Vector)
	paths = []
	for s in sources
		t = dijkstra_shortest_paths(g, s);
		p = enumerate_paths(t,destinations)
		for pp in p
			if length(pp) > 0
				println("found path")
				push!(paths,pp)
			end
		end
	end
	paths
end

function printProps(g::AbstractMetaGraph,paths::Vector,entSymbol::Symbol,strMap::Function)
	for p in paths
		println("\n\n")
		pr = props(g,p[1])
		pstring = strMap(pr[entSymbol])
		for v in p[2:end]
			pr = props(g,v)
			pstring = string(pstring,"▶\n",strMap(pr[entSymbol]))
		end
		println(pstring)
	end
end
