function annotatePathway!(dbParams::Dict,graph::AbstractMetaGraph)
	# get the uniprot ids of all proteins in the graph
	verts = filterVertices(graph,:entIdDb,p->occursin("uniprot",p))
	uniprot_ids = map(v->props(graph,v)[:entId],verts)

	# get the ortholog info for all upids
	ogs = selectOrthologs(unique(uniprot_ids),dbParams)

	# update the graph
	for ind in 1:length(verts)
		vi = verts[ind]
		uid = uniprot_ids[ind]
		dstr = unique(ogs[!,:orthoDist])
		dists = Dict()
		for d in 1:length(dstr)
			di = dstr[d]
			ogd = @from i in ogs begin
				  @where (i.entId,i.orthoDist) == (uid,di)
				  @select i
				  @collect DataFrame
			  	end
			genes = []
			for g in unique(ogd[!,:orthoDBId])
				ogg = @from i in ogd begin
					  @where (i.orthoDBId) == (g)
					  @select i
					  @collect DataFrame
				  	end
				geneFeat = [:orthoGroup,:orthoDist,
							:orthoSpecies,:orthoDBId,
							:orthoHGNC,:orthoEnsembleGeneId]
				gene = filter(v->!ismissing(v.second),
								Dict{Symbol,Any}(Symbol.(names(ogg[1,geneFeat])).=>
									values(ogg[1,geneFeat])))
				upids = []
				for op in unique(ogg.orthoEntId)
					println("update ortho distance $di for vertex $vi, gene $g, protein $op")
					push!(upids,op)
				end
				gene[:orthoEntIds] = Tuple(upids)
				push!(genes,gene)
			end
			og = Dict(
					:orthoSpecies=>reduce(vcat,unique(ogd.orthoSpecies)),
					:orthoGroup=>reduce(vcat,unique(ogd.orthoGroup)),
					:members=>Tuple(genes))
			# push!(dists,og)
			dists[parse(Int64,di)]=og
		end

		msg = set_props!(graph,verts[ind],Dict(:orthoDist=>dists))
		@assert msg == true "failed to add orthodb for vertex upid"
	end
	ogs
end
