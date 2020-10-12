# identify biochemical reactions where the input is Dna and output is protein
function findTranscriptionBC(g::AbstractMetaGraph)
    ctrl = ["http://www.biopax.org/release/biopax-level3.owl#Catalysis",
 			"http://www.biopax.org/release/biopax-level3.owl#Control"]

    colnames = [:geneInd,:ctrlInd,:protInd,
                :geneRef,:ctrlRef,:protRef]
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
function getNP(fcnParams::Dict,df::DataFrame)
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

    annotation = parseSparqlResponse(resp)
	# PCquery.deConcat!(annotations,:goTerms,:entries)
	# map the annotations to terms and get the nodes

	ann = @from i in annotations begin
		@group i by i.uniprot into g
		@select {uniprot=key(g),goAncestor=unique(g.goAncestor),
			entry=unique(g.entry),goTerm=unique(g.goTerm),isoform=unique(g.isoform)}
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
