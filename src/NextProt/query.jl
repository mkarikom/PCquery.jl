function getNextProt(fcnParams::Dict)
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/NextProt/rq")
    str = open(f->read(f, String), string(rqDir,"/","getNextProt_gdb.rq"));
    turtle = Mustache.render(str,
            Dict{Any,Any}("goFilter"=>fcnParams[:goFilter],
						  "entityFilter"=>fcnParams[:entFilter]))
    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
              "Accept" => "application/sparql-results+xml"]
    fmt = "application/sparql-results+xml"
    resp = PCquery.requestTTL(fcnParams[:port],fcnParams[:protocol],fmt,
                        fcnParams[:host],fcnParams[:path],fcnParams[:method],
                        header,turtle)

    ann = parseSparqlResponse(resp)
end

function annotateGraph!(fcnParams::Dict,g::AbstractMetaGraph)
	simpleVerts = filterVertices(g,:entIdDb,p->occursin("uniprot",p))
	uniprot_ids = map(v->props(g,v)[:entId],simpleVerts)
	entFilter = delimitValues(map(r->string(split(r,"/")[end]),unique(uniprot_ids)),"unipage:")
	fcnParams[:entFilter] = entFilter
	ann = getNextProt(fcnParams)

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

	for r in 1:length(a[!,:goAncShort])
		for v in a[r,:verts]
			set_prop!(g, v, fcnParams[:annotationKey], a[r,:goAncShort])
		end
	end
	a
end

function getNextProtGene(fcnParams::Dict)
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/NextProt/rq")
    str = open(f->read(f, String), string(rqDir,"/","getNextProtGene_gdb.rq"));
    turtle = Mustache.render(str,
            Dict{Any,Any}("unipage"=>fcnParams[:unipage]))
    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
              "Accept" => "application/sparql-results+xml"]
    fmt = "application/sparql-results+xml"
    resp = PCquery.requestTTL(fcnParams[:port],fcnParams[:protocol],fmt,
                        fcnParams[:host],fcnParams[:path],fcnParams[:method],
                        header,turtle)

    ann = parseSparqlResponse(resp)
end

# add gene ids to all the proteins
function annotateGraphGene!(fcnParams::Dict,g::AbstractMetaGraph)
	prot_v = filterVertices(g,:participantType,i->i=="http://www.biopax.org/release/biopax-level3.owl#Protein",a->get(a,:entId,missing))
	fcnParams[:unipage]=delimitValues(unique(prot_v.val),"unipage:")
	ann = getNextProtGene(fcnParams::Dict)

	# get the shorthand and add to the result
	upids = map(x->split(x,"/")[end],ann[!,:uniprot])
	insertcols!(ann, 1, :upid=>upids)

	for v in 1:length(prot_v.ind)
		ind = findfirst(x->x==prot_v.val[v],ann[!,:upid])
		ens = split(ann[ind,:gene],"/")[end]
		set_props!(g, prot_v.ind[v], Dict(:ensId=>ens))
		setval = props(g,prot_v.ind[v])[:ensId]
		println("success: $setval")
	end
	ann
end
# # find the simple lig rec
# lrverts = filterVertices(g,:roleLR,v->true,p->get(p,:roleLR,missing)) # the vertices representing complexes
# # find the complexes which include a ligand, receptor or both, excluding simple lig rec
# cxlr = PCquery.getCxLR(g,lrverts.ind)

# traverse inward edges from a start node and get values of target features
# return a dict with a vector of verts for each key
function getLRtree(g::AbstractMetaGraph,dir::Symbol,
					start::Int,targetFeatures::Dict)
	# get the subgraph
	vertFound = Dict()
	ng = bfs_tree(g, start, dir=dir)
	connected_v = findall(x->x > 0,degree(ng))
	pc_sg,vmap = induced_subgraph(g,connected_v)
	for k in keys(targetFeatures)
		k_fnd = []
		for v in targetFeatures[k]
			res = filterVertices(pc_sg,k,p->p==v,p->get(p,k,missing)) # the vertices representing complexes
			vertFound[v]=sort(res.ind)
		end
	end
	(v=vertFound,g=pc_sg)
end

function getTxGraphs(g::AbstractMetaGraph,dir::Symbol,txList::DataFrame,targetFeatures::Dict)
	# pull the protein ind (should already have ENSG id for zebrafish)
	paths = []
	for i in 1:size(txList,1)
		protInd = txList[i,:protInd]
		lrt = getLRtree(g,:in,protInd,targetFeatures)

		# find the product protInd in the new lrt graph
		inds = []
		for v in 1:nv(lrt.g)
			prp = props(lrt.g,v)
			if prp == props(g,protInd)
				push!(inds,v)
			end
		end
		orig = reduce(vcat,inds)
		p = dijkstra_shortest_paths(reverse(lrt.g), orig);
		pathverts = []
		allverts = Vector{eltype(collect(vertices(g)))}()
		dstverts = []
		for ft in targetFeatures[:roleLR]
			dst = filterVertices(lrt.g,:roleLR,f->f==ft)
			ep = enumerate_paths(p,dst)
			push!(pathverts,(ft,unique(reduce(vcat,ep))))
			push!(allverts,reduce(vcat,ep)...)
			push!(dstverts,(ft,unique(dst)))
		end
		allverts = unique(allverts)

		# construct the lig rec targ graph
		lrtg,vmap = induced_subgraph(lrt.g,allverts)
		recind = filterVertices(lrtg,:roleLR,f->f=="receptor",p->get(p,:ensId,missing))
		ligind = filterVertices(lrtg,:roleLR,f->f=="ligand",p->get(p,:ensId,missing))
		targgene = props(g,txList[i,:geneInd])[:entId]
		targind = filterVertices(lrtg,:ensId,f->f==targgene)
		recgenes = unique(recind.val)
		liggenes = unique(ligind.val)
		push!(paths,(recind=recind,
					 ligind=ligind,
					 targind=targind,
					 recgene=recgenes,
					 liggene=liggenes,
					 targgene=targgene,
					 graph=lrtg))
	end
	paths
end

function plotTxGraphs(g,paths,drnm,cols,lw)
	for i in 1:length(paths)
		i = 1
		recinds = paths[i].recind.ind
		liginds = paths[i].ligind.ind
		targind = paths[i].targind
		gg = paths[i].graph
		layout=(args...)->spring_layout(args...; C=20)
		#nodesize = log.([LightGraphs.degree(gg, v) for v in LightGraphs.vertices(gg)])
		nodesize=fill(1,nv(gg))
		#alphas = nodesize/maximum(nodesize)
		alphas=fill(1,nv(gg))

		nodefillc = fill(cols[3],nv(gg))
		nodefillc[targind] .= cols[4]
		nodefillc[liginds] .= cols[1]
		nodefillc[recinds] .= cols[2]

		nodelabel = fill("",nv(gg))
		nodelabel[recinds] .= string.("REC:",paths[i].recind.val)
		nodelabel[liginds] .= string.("LIG:",paths[i].ligind.val)
		nodelabel[targind] .= string.("TARG:",paths[i].targgene)

		gp = gplot(gg, layout=layout,
			  arrowlengthfrac=0.02,
			  arrowangleoffset = π/8,
			  nodesize=nodesize, nodefillc=nodefillc, nodelabel=nodelabel)
		draw(PDF(string(drnm,"path_$i","_",paths[i].targgene,"_path.pdf"), lw, lw), gp)
	end
end
