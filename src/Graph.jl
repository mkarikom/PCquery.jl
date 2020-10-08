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
function initBpGraph(df::DataFrame,
					pw::Symbol,
					lRef::Symbol,rRef::Symbol,ctrlRef::Symbol,rxn::Symbol,
					lRefType::Symbol,rRefType::Symbol,ctrlRefType::Symbol,
					lType::Symbol,rType::Symbol,ctrlType::Symbol,rxnType::Symbol,
					lLoc::Symbol,rLoc::Symbol,ctrlLoc::Symbol)
	len = size(df,1)
	lnodes = []
	rnodes = []
	allpairs = sort(collect(skipmissing(unique(vcat(
				df[!,lRef],df[!,rRef],df[!,rxn],df[!,ctrlRef])))))

	for i in 1:len
		# lRef vertex
		l = df[!,lRef][i]
		l_ind = findfirst(n->n==l,allnodes)
		set_props!(g, l_ind, Dict(:type=>df[!,lType][i], :location =>df[!,lLoc][i]))

		# rRef vertex
		r = df[!,rRef][i]
		r_ind = findfirst(n->n==r,allnodes)
		set_props!(g, r_ind, Dict(:type=>df[!,rType][i], :location =>df[!,rLoc][i]))

		# l/r edge properies
		rx = df[!,rxn][i]
		rx_ind = findfirst(n->n==rx,allnodes)
		add_edge!(g,l_ind,rx_ind)
		add_edge!(g,rx_ind,r_ind)
		set_prop!(g, Edge(l_ind, rx_ind), :interaction, df[!,rxnType][i])
		set_prop!(g, Edge(rx_ind, r_ind), :interaction, df[!,rxnType][i])

		# optional ctrlRef vertex and edge properties
		if !ismissing(df[!,ctrlRef][i])
			ct = df[!,ctrlRef][i]
			ct_ind = findfirst(n->n==ct,allnodes)
			set_props!(g, ct_ind, Dict(:type=>df[!,ctrlType][i], :location =>df[!,ctrlLoc][i]))
			add_edge!(g,ct_ind,rx_ind)
			set_prop!(g, Edge(ct_ind, rx_ind), :interaction, df[!,ctrlType][i])
		end
	end


	g = MetaDiGraph(length(allnodes))
	return Dict(:g=>g,:allnodes=>allnodes)
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
