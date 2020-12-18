# forms a recursive loop with decomposeComplex
# get a vector of simple physical entity dictionaries
function getNested(dbParams::Dict,entity::String,
					entNames::Vector,entNamesSimple::Vector;verbose=false)
	reflist = []
	isnested = true
	dbParams[:rqFile] = dbParams[:nestrqfile]
	resp = getParticipant(dbParams,entity)
	###
	# Stopping criteria for recursion: have found physical entity
	###
	if size(resp,1) == 0 # simple ref was provided
		verbose ? println("simple entity detected") : nothing
		isnested = false
		dbParams[:rqFile] = dbParams[:simprqfile]
		resp = getParticipant(dbParams,entity)
	end

	# loop over refs and recurse any complexes
	for i in 1:size(resp,1)
		cxtype = "http://www.biopax.org/release/biopax-level3.owl#Complex"
		unknowntype = "http://www.biopax.org/release/biopax-level3.owl#PhysicalEntity"
		# determine if resp[i,:] is a complex or a non-complex
		if resp.participantType[i] == cxtype
			verbose ? println("adding complex row $i") : nothing
			cxref =  resp[i,:participant]
			members = decomposeComplex(dbParams,cxref)
            cxdict = filterKeys(resp[i,:],entNames,dbParams[:physicalEntity])
			cxdict[:members] = members
			push!(reflist,cxdict)
		elseif resp.participantType[i] == unknowntype
			verbose ? println("adding generic to $i") : nothing
			cxdict = filterKeys(resp[i,:],
								[entNames...,entNamesSimple...],
								[dbParams[:physicalEntity]...,dbParams[:simpleEntity]...])
			push!(reflist,cxdict)
		elseif ismissing(resp.entId[i])
			verbose ? println("expanding family row $i") : nothing
			famref = resp[i,:participant]
			members = getNested(dbParams,famref,
								entNames,entNamesSimple)
			cxdict = filterKeys(resp[i,:],entNames,dbParams[:physicalEntity])
			cxdict[:members] = members
			push!(reflist,cxdict)
		else
			verbose ? println("adding participant row $i") : nothing
			cxdict = filterKeys(resp[i,:],
								[entNames...,entNamesSimple...],
								[dbParams[:physicalEntity]...,dbParams[:simpleEntity]...])
			push!(reflist,cxdict)
		end
	end
	# @assert size(reflist,1) > 0 "require at least one simple ref"
    reflist
end

# given a record (usually dataframe row) and a set of keys, return a dict with the found keys
function filterKeys(rec,keys,filterkeys)
	data = collect(rec[filterkeys])
	inds = findall((!ismissing).(data))
	d = Dict{Any,Any}(keys[inds].=>data[inds])
end

# step 1. an array of interaction dictionaries
# each interaction dictionary has:
# 	1) metadata for the edge
#	2) two member dictionaries
#		2.1) source member dictionary
#		2.2) destination member dictionary
# step 2. replace interaction dictionaries with nested members by a collection of interaction dictionaries over simple members
#
function expandInteractions(df::DataFrame,dbParams;verbose=false)
	# compose the db params
	ctrlExclude = ["http://www.biopax.org/release/biopax-level3.owl#Pathway","http://www.biopax.org/release/biopax-level3.owl#PhysicalEntity"]
	nestedParams = Dict{Symbol,Any}(
		:interaction=>[:partPred,:intDisplayName],
		:biochemInteraction=>[:interaction,:intType],
		:ctrlInteraction=>[:ctrlRxn,:ctrlRxnType,:ctrlRxnDir,:ctrlRxnDisplayName],
		:physicalEntity=>[:participantType,:participantLocRef,:participant,:displayName],
		:simpleEntity=>[:participantEntRef,:participantEntRefType,:entId,:entIdDb,:standardName])
	# form a add pc entity-type-wise vertex property keys to the provided connection params (dbParams)
	[nestedParams[k] = v for (k,v) in dbParams]

	refs = []
	rxns = []
	for i in 1:size(df,1)
		verbose ? println("getting rxn $i of ", size(df,1)) : nothing
		rxn = filterKeys(df[i,:],
				   nestedParams[:interaction],
				   nestedParams[:interaction])
		biochemRxn = filterKeys(df[i,:],
				   		nestedParams[:biochemInteraction],
				   		nestedParams[:biochemInteraction])

		# push the participant rxn dictionary to the ref list
		push!(refs,biochemRxn)
		ctrlRxn = filterKeys(df[i,:],
				   nestedParams[:ctrlInteraction],
				   nestedParams[:ctrlInteraction])

		entRef = df[i,:participant] # the ref for the participant entity

		nonNested = getNested(nestedParams,entRef,
						nestedParams[:physicalEntity],
						nestedParams[:simpleEntity])
		for ii in 1:length(nonNested)
			# push the participant simple entity dictionary to the ref list
			push!(refs,nonNested[ii])
			if length(collect(skipmissing(values(ctrlRxn)))) > 0 && !any(map(x->x==df[i,:ctrlEntityType],ctrlExclude))

				# push the ctrl rxn to the ref list
				push!(refs,ctrlRxn)
				ctrlEntRef = df[i,:ctrlEntity] # the ref for the ctrl entity
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
	(refs=refs,rxns=rxns)
	# (unique(refs),rxns)
end

# calls initGraph(df::DataFrame,dbParams::Dict)
# after getPathways(pathwayParams)
# given:
#  1) a set of pathway refs
#  2) dbParams
function initGraph(refs::Vector,dbParams::Dict)
	dbParams[:refs] = refs
	dbParams[:rqFile] = dbParams[:pwrqfile]
	df_pairs = getPathways(dbParams)
	refs,rxns = expandInteractions(df_pairs,dbParams)
	gd = initGraph(refs,rxns)
end

# initialize the graph
function initGraph(refs,rxns;verbose=false)
	# what reactions to consider
	rxtypes = ["http://www.biopax.org/release/biopax-level3.owl#ComplexAssembly",
			   "http://www.biopax.org/release/biopax-level3.owl#TemplateReaction",
			   "http://www.biopax.org/release/biopax-level3.owl#Transport",
			   "http://www.biopax.org/release/biopax-level3.owl#BiochemicalReaction",
			   "http://www.biopax.org/release/biopax-level3.owl#TransportWithBiochemicalReaction"]

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
		if any(occursin.(rx_dict[:intType],rxtypes))
			if base_dict[:partPred] == "http://www.biopax.org/release/biopax-level3.owl#left"
				addEdge!(g,p_ind,rx_ind,Dict(:eType=>"biochemInput"))
			elseif base_dict[:partPred] == "http://www.biopax.org/release/biopax-level3.owl#right"
				addEdge!(g,rx_ind,p_ind,Dict(:eType=>"biochemOutput"))
			elseif base_dict[:partPred] == "http://www.biopax.org/release/biopax-level3.owl#product"
				addEdge!(g,rx_ind,p_ind,Dict(:eType=>"templateProduct"))
			end
		else
			rtype = rx_dict[:intType]
			verbose ? println("rxn $i unsupported rx type $rtype") : nothing
		end

		# add the control rxn if any
		if haskey(rxns[i],:ctrlRxn)
			ctrl_rx_dict = rxns[i][:ctrlRxn]
			ctrl_p_dict = rxns[i][:ctrlEntity]
			ctrl_rx_ind = findVertex(g,ctrl_rx_dict)
			ctrl_p_ind = findVertex(g,ctrl_p_dict)
			addEdge!(g,ctrl_p_ind,ctrl_rx_ind,Dict(:eType=>"controller"))
			addEdge!(g,ctrl_rx_ind,rx_ind,Dict(:eType=>"catalysis"))
		end
	end

	# initialize the graph vertices
	for i in 1:length(refs)
		if haskey(refs[i],:participantType)
			if refs[i][:participantType] == "http://www.biopax.org/release/biopax-level3.owl#Complex"
				masterlist = Vector{Vector{Int64}}()
				appendTree!(g,i,masterlist,Dict(:eType=>"cxTree"))
				set_props!(g,i,Dict(:vertShell=>masterlist))
				verbose ? println(string("updated props for cx $i: ",reduce(vcat,masterlist))) : nothing
			end
		end
	end

	Dict(:graph=>g,:vertices=>refs,:simple=>rxns)
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

# mutate (extend) g with a directed tree rooted at the participant complex
function appendTree!(g::AbstractMetaGraph,parentInd::Int,edgeProps)
	children = props(g,parentInd)[:members]
	for c in children
		cInd = addVertex!(g,c)
		addEdge!(g,cInd,parentInd,edgeProps)
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
		addEdge!(g,cInd,parentInd,edgeProps)
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
			push!(clist,i...)
		else
			unpackNodes!(i,masterlist)
		end
	end
	push!(masterlist,clist)
end
