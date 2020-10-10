function getEdges(edgelist,edgelabels)
    list = Int.(zeros(length(edgelist)))
    for n in 1:length(edgelabels)
        list[findall(x->x==edgelabels[n],edgelist)] .= n
    end
    list
end

function loadGraph(left,right,allnodes)
    g = SimpleDiGraph(length(allnodes))
    for i in 1:length(left)
        add_edge!(g, left[i], right[i]);
    end
    g
end

function mapNodeNames(nodes::Vector,mapping::Dict)
	out = []
	# for n in nodes
	#
	# end
end
# get left,right using the value of surveilance rows to identify:
# each row will always have a "left" and a "right" (for the head and tail of a pathway these carry special importance)
# the head of the pathway as a "left" where the surveilance column is "lig_rec_complex"
# the tail of the pathway as a "right" where the surveilance column is "transcription complex"
function nodeLabels(df::DataFrame,
					left::Vector{Symbol},right::Vector{Symbol},
					leftloc::Vector,rightloc::Vector,
					keyvals...)
	len = size(df,1)
	allpairs = []
	lnodes = []
	rnodes = []
	for r in 1:len
		lnode = vec(convert(Array,df[!,left][r,:]))
		lnode = string(lnode[(!ismissing).(lnode)][1])
		lnode = string(lnode,"-",leftloc[r])

		rnode = vec(convert(Array,df[!,right][r,:]))
		rnode = string(rnode[(!ismissing).(rnode)][1])
		rnode = string(rnode,"-",rightloc[r])

		push!(allpairs,[lnode,rnode])
		push!(lnodes,lnode)
		push!(rnodes,rnode)
	end
	allnodes = sort(unique(vcat(lnodes,rnodes)))

	linds = []
	rinds = []
	for l in lnodes
		push!(linds,findfirst(x->x==l,allnodes))
	end
	for r in rnodes
		push!(rinds,findfirst(x->x==r,allnodes))
	end

	g = SimpleDiGraph(length(allnodes))

	for n in 1:length(linds)
		add_edge!(g,linds[n],rinds[n])
	end

	println(allnodes[1])

	return Dict(:g=>g,:lnodes=>lnodes,:rnodes=>rnodes,:allnodes=>allnodes)
end

# initialize a set of directed edges between pairs of vertices (bp:left/bp:right)
# add edge properties from edgeProps::Vector{Symbol}
# add vertex properties from vertexProps, where each element is a 2-element vector of left,right symbol
function initGraph(df::DataFrame,
					left::Symbol,right::Symbol,
					edgeProps::Vector,
					vertexProps::Vector)
	len = size(df,1)
	allpairs = []
	lnodes = []
	rnodes = []
	for r in 1:len
		lnode = df[!,left][r]
		rnode = df[!,right][r]
		push!(lnodes,lnode)
		push!(rnodes,rnode)
	end
	allnodes = sort(unique(vcat(lnodes,rnodes)))

	linds = []
	rinds = []
	for l in lnodes
		push!(linds,findfirst(x->x==l,allnodes))
	end
	for r in rnodes
		push!(rinds,findfirst(x->x==r,allnodes))
	end

	# define the digraph
	g = SimpleDiGraph(length(allnodes))
	for n in 1:length(linds)
		add_edge!(g,linds[n],rinds[n])
	end
	mg = MetaDiGraph(g)

	# add edge props
	for p in edgeProps
		for e in 1:len
			if !ismissing(df[!,p[1]][e])
				l = df[!,left][e]
				l_ind = findfirst(n->n==l,allnodes)
				r = df[!,right][e]
				r_ind = findfirst(n->n==r,allnodes)
				set_prop!(mg, Edge(l_ind, r_ind), p[2], df[!,p[1]][e])
			end
		end
	end

	# add vertex props
	for p in vertexProps
		for v in 1:len
			l = df[!,left][v]
			l_ind = findfirst(n->n==l,allnodes)
			r = df[!,right][v]
			r_ind = findfirst(n->n==r,allnodes)
			if !ismissing(df[!,collect(keys(p))[1]][v])
				set_prop!(mg, l_ind, collect(values(p))[1], df[!,collect(keys(p))[1]][v])
			end
			if !ismissing(df[!,collect(keys(p))[2]][v])
				set_prop!(mg, r_ind, collect(values(p))[2], df[!,collect(keys(p))[2]][v])
			end
		end
	end

	return Dict(:g=>mg,:lnodes=>lnodes,:rnodes=>rnodes,:allnodes=>allnodes)
end

# initialize a bipartite graph where {left,right,control} ʌ {interactions} = ∅
# for edge [1,2], direction is 1->2.  define: [left,interaction], [interaction,right], [ctrl,interaction]
function initBpGraph(df::DataFrame,pw::Symbol,nestedParams,
					 participantLocRef::Symbol,participantType::Symbol,participantRefs::Symbol,participantRefTypes::Symbol,
					 partPred::Symbol,intRef::Symbol,intType::Symbol,
					 ctrlLocRef::Symbol,ctrlRxn::Symbol,ctrlRxnType::Symbol,ctrlRxnDir::Symbol,
					 ctrlEntityType::Symbol,ctrlEntityRefs::Symbol,ctrlEntityRefTypes::Symbol)
	len = size(df,1)
	# collect verts
	inds=[];rxns=[];participants=[];preds=[];ctrlRxns=[];ctrlParticipants=[]
	# criteria for each nested entity
	simpleF = r->occursin("Unification",r)
	nestedF = r->!occursin("Unification",r)
	# criteria for participant edge direction
	edgeSym = partPred
	symOut = "http://www.biopax.org/release/biopax-level3.owl#left"
	symIn = "http://www.biopax.org/release/biopax-level3.owl#right"
	for i in 1:len
		p_nested = getNested(df,nestedParams,i,
						simpleF,nestedF,
						participantRefTypes,participantRefs)
		for ii in 1:length(p_nested)
			p_i=p_nested[ii]
			if !ismissing(df[!,ctrlRxn][i])
				c_nested = getNested(df,nestedParams,i,
								simpleF,nestedF,
								ctrlEntityRefTypes,ctrlEntityRefs)
				for iii in 1:length(c_nested)
					c_i = c_nested[iii]
					push!(ctrlParticipants,c_i)
					push!(ctrlRxns,df[!,ctrlRxn][i])
					push!(inds,i)
					push!(participants,p_i)
					push!(rxns,df[!,intRef][i])
					push!(preds,df[!,edgeSym][i])
				end
			else
				push!(ctrlParticipants,missing)
				push!(ctrlRxns,missing)
				push!(inds,i)
				push!(participants,p_i)
				push!(rxns,df[!,intRef][i])
				push!(preds,df[!,edgeSym][i])
			end
		end
	end


	allverts = sortUnique(rxns,participants,ctrlRxns,ctrlParticipants)
	g = MetaDiGraph(length(allverts))

	for i in 1:length(inds)
		ind = inds[i]
		# participants
		p_ind = findfirst(n->n==participants[i],allverts)
		set_props!(g, p_ind, Dict(:reactome=>participants[i],
								  :entityType=>df[!,participantType][inds[i]],
								  :location=>df[!,participantLocRef][inds[i]]))

		rx_ind = findfirst(n->n==rxns[i],allverts)
		set_props!(g,rx_ind,Dict(:reactome=>rxns[i],
								  :entityType=>df[!,intType][inds[i]]))

		if preds[i] == symIn
			add_edge!(g,p_ind,rx_ind)
		elseif preds[i] == symOut
			add_edge!(g,rx_ind,p_ind)
		else
			throw("only left/right reactions are supported")
		end

		# optional ctrlRef vertex and edge properties
		if !ismissing(ctrlRxns[i])
			if ismissing(df[!,ctrlRxnDir][inds[i]])
				throw("ctrl rxn dir missing")
			end
			ct = ctrlParticipants[i]
			ct_rx = ctrlRxns[i]
			ct_ind = findfirst(n->n==ct,allverts)
			ct_rx_ind = findfirst(n->n==ct_rx,allverts)
			set_props!(g, ct_ind, Dict(:reactome=>ctrlParticipants[i],
									   :entityType=>df[!,ctrlEntityType][inds[i]],
									   :location=>df[!,ctrlLocRef][inds[i]]))

			set_props!(g, ct_rx_ind, Dict(:reactome=>ctrlRxns[i],
									   	  :entityType=>df[!,ctrlRxnType][inds[i]],
										  :controlType=>df[!,ctrlRxnDir][inds[i]]))

			add_edge!(g,ct_ind,ct_rx_ind)
			add_edge!(g,ct_rx_ind,rx_ind)
		end
	end
	Dict(:graph=>g,:vertices=>allverts)
end

function getAllPaths(paths,startnodes,endnodes)
	foundpaths = []
	for startnode in startnodes
		for endnode in endnodes # start/end pairs
			for p in paths[startnode]
				if length(p) > 0
					if p[1] == startnode && p[end] == endnode
						push!(foundpaths,p)
					end
				end
			end
		end
	end
	foundpaths
end

function collectVertex!(refereceVec::Vector,
					   vertexRef::String)
   ind = findfirst(n->n==vertexRef,refereceVec)
   if isnothing(ind)
	   push!(referenceVec,vertexRef)
	   ind = length(referenceVec)
   end
   ind
end
