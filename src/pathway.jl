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

# calls initGraph(df::DataFrame,dbParams::Dict)
# after getPathways(pathwayParams)
# given:
#  1) a set of pathway refs
#  2) dbParams
function initGraph(refs::Vector,dbParams::Dict)
	pathwayParams = copy(dbParams)
	pathwayParams[:refs] = refs
	df_pairs = getPathways(pathwayParams)
	gd = initGraph(df_pairs,dbParams)
end

# initialize the graph
function initGraph(df::DataFrame,dbParams::Dict)
	# compose the db params
	nestedParams = Dict{Symbol,Any}(
		:interaction=>[:partPred],
		:biochemInteraction=>[:interaction,:intType,:displayNameIntxn],
		:ctrlInteraction=>[:ctrlRxn,:ctrlRxnType,:ctrlRxnDir],
		:physicalEntity=>[:participantType,:participantRef,:participantLocRef,:participant,:displayName],
		:simpleEntity=>[:participantEntRef,:participantEntRefType,:entId,:entIdDb],
		:ctrlPhysicalEntity=>[:ctrlEntityType,:ctrlEntityRef,:ctrlEntityLocRef,:ctrlEntity],
		:ctrlSimpleEntity=>[:ctrlEntityEntRef,:ctrlEntityEntRefType,:ctrlEntityEntId,:ctrlEntityEntIdDb])
	# form a add pc entity-type-wise vertex property keys to the provided connection params (dbParams)
	[nestedParams[k] = v for (k,v) in dbParams]

	# get interactions
	refs,rxns = expandInteractions(df::DataFrame,nestedParams)

	# criteria for participant edge direction
	symOut = "http://www.biopax.org/release/biopax-level3.owl#right"
	symIn = "http://www.biopax.org/release/biopax-level3.owl#left"

	# initialize the graph vertices
	g = MetaDiGraph(length(refs))
	for i in 1:length(refs)
		added = set_props!(g,i,refs[i])
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

	# initialize the graph vertices
	for i in 1:length(refs)
		if haskey(refs[i],:participantType)
			if refs[i][:participantType] == "http://www.biopax.org/release/biopax-level3.owl#Complex"
				masterlist = Vector{Vector{Int64}}()
				appendTree!(g,i,masterlist,Dict(:eType=>"cxTree"))
				set_props!(g,i,Dict(:vertShell=>masterlist))
				println(string("updated props for cx $i: ",reduce(vcat,masterlist)))
			end
		end
	end

	Dict(:graph=>g,:vertices=>refs,:simple=>rxns)
end

# forms a recursive loop with decomposeComplex
# get a vector of simple physical entity dictionaries
function getNested(dbParams::Dict,entity::String,
					entNames::Vector,entNamesSimple::Vector)
	reflist = []
	dbParams[:rqFile] = "getParticipantNested_gdb.rq"
	resp = getParticipant(dbParams::Dict,entity::String)
	###
	# Stopping criteria for recursion: have found physical entity
	###
	if size(resp,1) == 0 # simple ref was provided
		dbParams[:rqFile] = "getParticipantSimple_gdb.rq"
		resp = getParticipant(dbParams::Dict,entity::String)
	else # simple refs were retrieved from nested
		indMissing = findall(n->ismissing(n),collect(resp[1,:]))
		@assert length(indMissing) == 0 "all nested participants must be resolved"
	end
	# loop over refs and recurse any complexes
	for i in 1:size(resp,1)
		simpleNames = [dbParams[:physicalEntity]...,dbParams[:simpleEntity]...]
		inds = indexin(simpleNames,Symbol.(names(resp)))
		ind_found = findall(r->!isnothing(r),inds)
		ind_notfound = findall(r->isnothing(r),inds)
		# determine if resp[i,:] is a complex or a non-complex
		if length(ind_found) < length(simpleNames)
			println("adding complex row $i")
			cxref =  resp[i,dbParams[:physicalEntity][4]]
			members = decomposeComplex(dbParams,cxref)
            cxdict = Dict([entNames...,:members] .=> [collect(resp[i,dbParams[:physicalEntity]])...,members])
			push!(reflist,cxdict)
		else
			println("adding participant row $i")
			push!(reflist,Dict([entNames...,entNamesSimple...] .=>
				collect(resp[i,[dbParams[:physicalEntity]...,dbParams[:simpleEntity]...]])))
		end
	end
	# @assert size(reflist,1) > 0 "require at least one simple ref"
    reflist
end

# forms a recursive loop with getNested
# decompose all complexes in the interaction list
function decomposeComplex(dbParams::Dict,cxref::String)
	cxRefs = getComplex(dbParams,cxref)
	members = []
	for m in 1:size(cxRefs,1)
        # need to include getNested here, temporarily omit
        nonNested = getNested(dbParams,cxRefs[m,:comp],
                                dbParams[:physicalEntity],
                                dbParams[:simpleEntity])
        for nn in 1:length(nonNested)
            push!(members,nonNested[nn])
        end
	end
	members
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

# mutate (extend) g with a directed tree rooted at the participant complex
function appendTree!(g::AbstractMetaGraph,parentInd::Int,edgeProps)
	children = props(g,parentInd)[:members]
	for c in children
		cInd = addVertex!(g,c)
		addEdge!(g,cInd,parentInd,edgeProps,false)
		if haskey(c,:members)
			appendTree!(g,cInd)
			rem_prop!(g,cInd,:members)
		end
	end
end

# mutate (extend) g with a directed tree rooted at the participant complex
function appendTree!(g::AbstractMetaGraph,parentInd::Int,masterlist,edgeProps)
	nlist = []
	traverseNodeDict!(g,parentInd,nlist,edgeProps)
	unpackNodes!(nlist,masterlist)
	push!(masterlist,[parentInd])
	reverse!(masterlist)
end


# mutate (extend) g with a directed tree rooted at the participant complex
# return a nested vector which will later be unpacked
function traverseNodeDict!(g::AbstractMetaGraph,parentInd::Int,nlist,edgeProps)
	children = props(g,parentInd)[:members]
	c_nlist=[]
	for c in children
		cInd = addVertex!(g,c)
		addEdge!(g,cInd,parentInd,edgeProps,false)
		if haskey(c,:members)
			push!(nlist,traverseNodeDict!(g,cInd,c_nlist,edgeProps))
			# rem_prop!(g,cInd,:members)
		end
		push!(nlist,cInd)
	end
	nlist
end

function unpackNodes!(nlist,masterlist)
	clist = Vector{Int64}()
	for i in nlist
		if length(i) == 1
			push!(clist,i)
		else
			unpackNodes!(i,masterlist)
		end
	end
	push!(masterlist,clist)
end
