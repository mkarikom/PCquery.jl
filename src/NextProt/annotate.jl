function annotateGraphFcn!(fcnParams::Dict,g::AbstractMetaGraph;verbose=false)
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
				verbose ? println("setting $annot = $npval for vertex $v") : nothing
				set_prop!(g, v, annot, npval)
			end
		end
		return a
	else
		verbose ? println(string("did not find any functional annotations matching filter: ",fcnParams[:goFilter])) : nothing
	end
	nothing
end



# add gene names and ensids for all the upids
function annotateGraphP!(fcnParams::Dict,g::AbstractMetaGraph;verbose=false)
	filt = [v->!haskey(v,:members),
			v->haskey(v,:participantType) && v[:participantType] == "http://www.biopax.org/release/biopax-level3.owl#Protein"]
	v_ind = filterVertices(g,filt)
	v_val = map(v->get(props(g,v),:entId,missing),v_ind)
	if length(v_ind) > 0
		fcnParams[:unipage]=delimitValues(unique(v_val),"unipage:")
		ann = getNextProtP(fcnParams)
		if size(ann,1) > 0
			upids = map(x->split(x,"/")[end],ann[!,:uniprot])
			insertcols!(ann, 1, :upid=>upids)

			for v in 1:length(v_ind)
				if v==133
					verbose ? println("stop") : nothing
					verbose ? println("ind v ",v_ind[v]) : nothing
				end
				verbose ? println("ind v ",v_ind[v]) : nothing
				ind = findfirst(x->x==v_val[v],ann[!,:upid])
				ens = split(ann[ind,:gene],"/")[end]
				geneName = ann[ind,:gname]
				status = set_props!(g, v_ind[v], Dict(:ensId=>ens,
												  :gname=>geneName))
				if status
					verbose ? println("success: ensId=$ens, gene name=$geneName") : nothing
				else
					verbose ? println("failed to set :ensId or :gname on ",v_ind[v]) : nothing
				end
			end
			verbose ? println("protein annotation complete") : nothing
			return ann
		else
			verbose ? println("did not find any protein annotations") : nothing
		end
	end
	nothing
end

# add gene names for all the ensids
function annotateGraphG!(fcnParams::Dict,g::AbstractMetaGraph;verbose=false)
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
					verbose ? println("success: geneName=$geneName") : nothing
				else
					verbose ? println("failed to set :gname on ",gene_v.ind[v]) : nothing
				end
			end
			verbose ? println("gene annotation complete") : nothing
			return ann
		else
			verbose ? println("did not find any gene annotations") : nothing
		end
	end
	nothing
end
