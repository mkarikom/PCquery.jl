function annotateGraphFcn!(fcnParams::Dict,g::AbstractMetaGraph)
	simpleVerts = filterVertices(g,:entIdDb,p->occursin("uniprot",p))

	uniprot_ids = map(v->props(g,v)[:entId],simpleVerts)
	entFilter = delimitValues(map(r->string(split(r,"/")[end]),unique(uniprot_ids)),"unipage:")
	fcnParams[:entFilter] = entFilter
	ann = getNextProtFcn(fcnParams)
	if size(ann,1) > 0
		# get the shorthand and add to the result
		terms = fcnParams[:annotationMapAncTerm](ann)
		insertcols!(ann, 3, :goAncShort=>terms)

		ctrlFcnVert = [collect(filter_vertices(g,:ctrlEntityEntId,up)) for up in map(iri->split(iri,"/")[end],ann.uniprot)]
		fcnVert = [collect(filter_vertices(g,:entId,up)) for up in map(iri->split(iri,"/")[end],ann.uniprot)]
		insertcols!(ann, 3, :ctrlFcnVert=>ctrlFcnVert)
		insertcols!(ann, 3, :fcnVert=>fcnVert)

		a = @from i in ann begin
			@group i by {i.goAncShort} into g
			@select {goAncShort=reduce(vcat,unique(g.goAncShort)),verts=vcat(g.fcnVert...)}
			@collect DataFrame
		end

		for r in 1:length(a[!,:goAncShort])
			for v in a[r,:verts]
				annot = fcnParams[:annotationKey]
				npval = a[r,:goAncShort]
				println("setting $annot = $npval for vertex $v")
				set_prop!(g, v, annot, npval)
			end
		end
		return a
	else
		println(string("did not find any functional annotations matching filter: ",fcnParams[:goFilter]))
	end
	nothing
end



# add gene names and ensids for all the upids
function annotateGraphP!(fcnParams::Dict,g::AbstractMetaGraph)
	prot_v = filterVertices(g,:participantType,i->i=="http://www.biopax.org/release/biopax-level3.owl#Protein",a->get(a,:entId,missing))
	if length(prot_v.ind) > 0
		fcnParams[:unipage]=delimitValues(unique(prot_v.val),"unipage:")
		ann = getNextProtP(fcnParams)
		if size(ann,1) > 0
			upids = map(x->split(x,"/")[end],ann[!,:uniprot])
			insertcols!(ann, 1, :upid=>upids)

			for v in 1:length(prot_v.ind)
				println("ind v $v")
				ind = findfirst(x->x==prot_v.val[v],ann[!,:upid])
				ens = split(ann[ind,:gene],"/")[end]
				geneName = ann[ind,:gname]
				status = set_props!(g, prot_v.ind[v], Dict(:ensId=>ens,
												  :gname=>geneName))
				if status
					println("success: ensId=$ens, gene name=$geneName")
				else
					println("failed to set :ensId or :gname on ",prot_v.ind[v])
				end
			end
			println("protein annotation complete")
			return ann
		else
			println("did not find any protein annotations")
		end
	end
	nothing
end

# add gene names for all the ensids
function annotateGraphG!(fcnParams::Dict,g::AbstractMetaGraph)
	gene_v = filterVertices(g,:participantType,i->i=="http://www.biopax.org/release/biopax-level3.owl#Dna",a->get(a,:entId,missing))
	if length(gene_v.ind) > 0
		fcnParams[:gene]=delimitValues(unique(gene_v.val),"gene:")
		ann = getNextProtG(fcnParams)
		if size(ann,1) > 0
			ensids = map(x->split(x,"/")[end],ann[!,:gene])
			insertcols!(ann, 1, :ensid=>ensids)

			for v in 1:length(gene_v.ind)
				ind = findfirst(x->x==gene_v.val[v],ann[!,:ensid])
				geneName = ann[ind,:gname]
				status = set_props!(g, gene_v.ind[v], Dict(:gname=>geneName))
				if status
					println("success: geneName=$geneName")
				else
					println("failed to set :gname on ",prot_v.ind[v])
				end
			end
			println("gene annotation complete")
			return ann
		else
			println("did not find any gene annotations")
		end
	end
	nothing
end
