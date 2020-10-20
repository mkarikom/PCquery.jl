# addEdge!(g,p_ind,rx_ind,int_df[i,:],true)
function addEdge!(g::AbstractMetaGraph,source::Int,dest::Int,
					eProps::Dict,debug::Bool)
	hasedge = has_edge(g, source, dest)
	if hasedge
		if debug
			println("\n\nedge $source -> $dest already exists")
		end
	else
		added = add_edge!(g,source,dest)
		set_props!(g,source,dest,eProps)
		@assert added "edge $source -> $dest error"
	end
end

function filterSubgraph(g::AbstractMetaGraph,props::Vector,values)
	lst = []
	for p in 1:length(props)
		for v in values[p]
			bcv = sort(unique(vcat(map(e->[src(e),dst(e)],collect(filter_edges(g,props[p],v)))...)))
			push!(lst,bcv...)
		end
	end
	sg,vmap = induced_subgraph(g,Int.(unique(lst)))
end

function findVertex(g::AbstractMetaGraph,vProps::Dict)
	inds = []
	for v in 1:nv(g)
		if props(g,v) == vProps
			push!(inds,v)
		end
	end
	reduce(vcat,inds)
end

# filter the verts, update a dict and return a simple vector of verts
function filterVertices!(g::AbstractMetaGraph,filter::Dict,resultDict::Dict)
	resultVec = Vector{eltype(collect(vertices(g)))}()
	for k in keys(filter)
		resultDict[k] = Dict()
		for val in get(filter,k,"")
			resultDict[k][val] = []
			for v in 1:nv(g)
				vp = props(g,v)
				if haskey(vp,k)
					if props(g,v)[k] == val
						push!(resultVec,v)
						push!(resultDict[k][val],v)
					end
				end
			end
		end
	end
end

# filter the edges, update a dict and return a simple vector of edges
function filterEdges!(g::AbstractMetaGraph,filter::Dict,resultDict::Dict)
	resultVec = Vector{eltype(collect(edges(g)))}()
	for k in keys(filter)
		resultDict[k] = Dict()
		for val in get(filter,k,"")
			resultDict[k][val] = []
			for e in edges(g)
				ep = props(g,e)
				if haskey(ep,k)
					if props(g,e)[k] == val
						push!(resultVec,e)
						push!(resultDict[k][val],e)
					end
				end
			end
		end
	end
	resultVec
end

# as above but evaluate a function on a key
function filterVertices!(g::AbstractMetaGraph,k::Symbol,f::Function)
	resultVec = Vector{eltype(collect(vertices(g)))}()
	for v in vertices(g)
		vp = props(g,v)
		if haskey(vp,k)
			if f(props(g,v)[k])
				push!(resultVec,v)
			end
		end
	end
	resultVec
end

# as above but evaluate a function on a key
function filterEdges!(g::AbstractMetaGraph,k::Symbol,f::Function)
	resultVec = Vector{eltype(collect(edges(g)))}()
	for e in edges(g)
		ep = props(g,e)
		if haskey(ep,k)
			if f(props(g,e)[k])
				push!(resultVec,e)
			end
		end
	end
	resultVec
end

# as above but also return a mapping of the value
function filterVertices!(g::AbstractMetaGraph,k::Symbol,f::Function,m::Function)
	resultVec = Vector{eltype(collect(vertices(g)))}()
	valVec = []
	for v in vertices(g)
		vp = props(g,v)
		if haskey(vp,k)
			if f(props(g,v)[k])
				push!(resultVec,v)
				push!(valVec,m(props(g,v)[k]))
			end
		end
	end
	(ind=resultVec,val=valVec)
end

# as above but also return a mapping of the value
function filterEdges!(g::AbstractMetaGraph,k::Symbol,f::Function,m::Function)
	resultVec = Vector{eltype(collect(edges(g)))}()
	valVec = []
	for e in edges(g)
		ep = props(g,e)
		if haskey(ep,k)
			if f(props(g,e)[k])
				push!(resultVec,e)
				push!(valVec,m(props(g,e)[k]))
			end
		end
	end
	# resultVec,valVec
	(ind=resultVec,val=valVec)
end
