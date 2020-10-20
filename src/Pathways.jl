# step 1. an array of interaction dictionaries
# each interaction dictionary has:
# 	1) metadata for the edge
#	2) two member dictionaries
#		2.1) source member dictionary
#		2.2) destination member dictionary
# step 2. replace interaction dictionaries with nested members by a collection of interaction dictionaries over simple members
#
function expandInteractions(df::DataFrame,nestedParams)
	refs = []
	rxns = []
	# criteria for participant edge direction
	symOut = "http://www.biopax.org/release/biopax-level3.owl#right"
	symIn = "http://www.biopax.org/release/biopax-level3.owl#left"
	for i in 1:size(df,1)
		# println("getting rxn $i")
		rxn = Dict(nestedParams[:interaction].=>collect(df[i,nestedParams[:interaction]]))
		biochemRxn = Dict(nestedParams[:biochemInteraction].=>collect(df[i,nestedParams[:biochemInteraction]]))
		# push the participant rxn dictionary to the ref list
		push!(refs,biochemRxn)
		ctrlRxn = Dict(nestedParams[:ctrlInteraction].=>collect(df[i,nestedParams[:ctrlInteraction]]))
		entRef = df[i,nestedParams[:physicalEntity][4]] # the ref for the participant entity
		nonNested = getNested(nestedParams,entRef,
						nestedParams[:physicalEntity],
						nestedParams[:simpleEntity])
		for ii in 1:length(nonNested)
			# push the participant simple entity dictionary to the ref list 
			push!(refs,nonNested[ii])
			if length(collect(skipmissing(values(ctrlRxn)))) > 0
				# push the ctrl rxn to the ref list
				push!(refs,ctrlRxn)
				ctrlEntRef = df[i,nestedParams[:ctrlPhysicalEntity][4]] # the ref for the ctrl entity
				ctrlNonNested = getNested(nestedParams,ctrlEntRef,
									nestedParams[:physicalEntity],
									nestedParams[:simpleEntity])
				for iii in 1:length(ctrlNonNested)
					# push the ctrl simple entity dictionary to the ref list
					push!(refs,ctrlNonNested[iii])
					# push the controlled rxn dictionary to the rxn list
					rxnDict = Dict(:rxn=>rxn,
								   :biochemRxn=>biochemRxn,
								   :participant=>nonNested[ii],
								   :ctrlRxn=>ctrlRxn,
								   :ctrlEntity=>ctrlNonNested[iii])
					push!(rxns,rxnDict)
				end
			else
				rxnDict = Dict(:rxn=>rxn,
							   :biochemRxn=>biochemRxn,
							   :participant=>nonNested[ii])
				# push the non-controlled rxn dictionary to the rxn list
				push!(rxns,rxnDict)
			end
		end
	end
	(unique(refs),rxns)
end

# initialize a bipartite graph where {left,right,control} ʌ {interactions} = ∅
# for edge [1,2], direction is 1->2.  define: [left,interaction], [interaction,right], [ctrl,interaction]
# because the graph is bipartite, all edges will be labeled with bp:BiochemicalReaction
function initBpGraph(df::DataFrame,nestedParams)
	# get interactions
	refs,rxns = expandInteractions(df::DataFrame,nestedParams)

	# criteria for participant edge direction
	symOut = "http://www.biopax.org/release/biopax-level3.owl#right"
	symIn = "http://www.biopax.org/release/biopax-level3.owl#left"

	# initialize the graph vertices
	g = MetaDiGraph(length(refs))
	for i in 1:length(refs)
		added = set_props!(g,i,refs[i])
		# println("updated props for vertex $i: $added")
	end

	# initialize edges (rxns)
	for i in 1:length(rxns)
		# add the primary rxn
		base_dict = rxns[i][:rxn]
		rx_dict = rxns[i][:biochemRxn]
		p_dict = rxns[i][:participant]
		rx_ind = findVertex(g,rx_dict)
		p_ind = findVertex(g,p_dict)
		if rx_dict[:intType] == "http://www.biopax.org/release/biopax-level3.owl#BiochemicalReaction"
			if base_dict[:partPred] == "http://www.biopax.org/release/biopax-level3.owl#left"
				addEdge!(g,p_ind,rx_ind,Dict(:eType=>"biochemInput"),false)
			elseif base_dict[:partPred] == "http://www.biopax.org/release/biopax-level3.owl#right"
				addEdge!(g,rx_ind,p_ind,Dict(:eType=>"biochemOutput"),false)
			end
		else
			rtype = rx_dict[:intType]
			println("rxn $i unsupported rx type $rtype")
		end

		# add the control rxn if any
		if haskey(rxns[i],:ctrlRxn)
			ctrl_rx_dict = rxns[i][:ctrlRxn]
			ctrl_p_dict = rxns[i][:ctrlEntity]
			ctrl_rx_ind = findVertex(g,ctrl_rx_dict)
			ctrl_p_ind = findVertex(g,ctrl_p_dict)
			addEdge!(g,ctrl_p_ind,ctrl_rx_ind,Dict(:eType=>"controller"),false)
			addEdge!(g,ctrl_rx_ind,rx_ind,Dict(:eType=>"catalysis"),false)
		end
	end

	Dict(:graph=>g,:vertices=>refs,:simple=>rxns)
end

# identify biochemical reactions where the input is Dna and output is protein
function getTransTargs(g::AbstractMetaGraph)
    ctrl = ["http://www.biopax.org/release/biopax-level3.owl#Catalysis",
 			"http://www.biopax.org/release/biopax-level3.owl#Control"]

    colnames = [:ctrlInd,:geneInd,:protInd,
                :ctrlRef,:geneRef,:protRef]
    coltypes = [Union{Missing,Int64},Union{Missing,Int64},Union{Missing,Int64},
                Union{Missing,String},Union{Missing,String},Union{Missing,String}]
    edgeTable = DataFrame(coltypes,colnames)

	for v in vertices(g)
		# identify biochemical rxn
		if haskey(props(g,v),:intType)
			if props(g,v)[:intType] == "http://www.biopax.org/release/biopax-level3.owl#BiochemicalReaction"
				ns_o = outneighbors(g,v)
		        ns_i = inneighbors(g,v)
				println("out neigh = ",length(ns_o),", in neigh = ",length(ns_i))
				for n_o in ns_o
					for n_i in ns_i
						if haskey(props(g,n_i),:participantType)
							if props(g,n_i)[:participantType] == "http://www.biopax.org/release/biopax-level3.owl#Dna"
								if haskey(props(g,n_o),:participantType)
									if props(g,n_o)[:participantType] == "http://www.biopax.org/release/biopax-level3.owl#Protein"
										rxref = props(g,v)[:interaction]
										inref = props(g,n_i)[:participant]
										outref = props(g,n_o)[:participant]
										data = tuple(v,n_i,n_o,inref,rxref,outref)
										push!(edgeTable,initRow(colnames,colnames,data))
									end
								else
									println("no output participant $n_o")
								end
							end
						else
							println("no input participant $n_i")
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
    str = open(f->read(f, String), "$rqDir/getPC_gdb.rq");
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
    rqDir = string(srcDir,"/rq")
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

function getNextProt(fcnParams::Dict,g::AbstractMetaGraph,df::Vector)
	simpleVerts = filterVertices!(g,:entIdDb,p->occursin("uniprot",p))
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
