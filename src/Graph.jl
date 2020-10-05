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
	# Ex.
	# keyvals[1] = Dict(:lr=>:left,
	# 					:key=>:llocref,
	# 				    :val=>"http://pathwaycommons.org/pc11/#UnificationXref_gene_ontology_GO_0005829")


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
