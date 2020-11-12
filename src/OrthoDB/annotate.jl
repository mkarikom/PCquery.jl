function annotatePathway!(graph::AbstractMetaGraph,dbParams::Dict)

	df = DataFrame()
	verts = []
	ids = []

	dbParams[:entpfx]="http://purl.uniprot.org/uniprot/"
	dbParams[:entpfxs]="uniprot:"
	un_verts = filterVertices(graph,:entIdDb,p->occursin("uniprot",p))
	un_ids = map(v->props(graph,v)[:entId],un_verts)
	selectOrthologs!(unique(un_ids),dbParams,df)

	ids = vcat(ids,un_ids)
	verts = vcat(verts,un_verts)

	dbParams[:entpfx]="http://rdf.ebi.ac.uk/resource/ensembl/"
	dbParams[:entpfxs]="ensembl:"
	en_verts = filterVertices(graph,:entIdDb,p->occursin("ensembl",p))
	en_ids = map(v->props(graph,v)[:entId],en_verts)
	selectOrthologs!(unique(en_ids),dbParams,df)
	ids = vcat(ids,en_ids)
	verts = vcat(verts,en_verts)

	addOrthoDist!(graph,verts,ids,df)

	(df=df,ids=ids,verts=verts)
end

function addOrthoDist!(graph::AbstractMetaGraph,verts,eids,ogs)
	# update the graph
	dstr = unique(ogs[!,:orthoDist])
	for ind in 1:length(verts)
		vi = verts[ind]
		uid = eids[ind]
		dists = Dict()
		for d in 1:length(dstr)
			di = dstr[d]
			ogd = @from i in ogs begin
				  @where (i.entId,i.orthoDist) == (uid,di)
				  @select i
				  @collect DataFrame
				end
			if size(ogd)[1] > 0
				println("adding dist ind $ind vert $vi uid $uid dist $di, ortholog count is ",size(ogd)[1])
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
						# println("update ortho distance $di for vertex $vi, gene $g, protein $op")
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
			else
				println("no orthologs reported for ind $ind vert $vi uid $uid at dist $di")
			end
		end

		if length(dists) > 0
			msg = set_props!(graph,verts[ind],Dict(:orthoDist=>dists))
			@assert msg == true "failed to add orthodb for vertex upid"
		else
			println("no orthologs reported for ind $ind vert $vi uid $uid")
		end
	end
end

# add gene expression for the given key:
# all orthoDBId's at all levels
function addExpression(graph::AbstractMetaGraph,
						df::DataFrame,
						key::Symbol,
						modality::Symbol,
						barcode::String)
	# get vertices with orthologs
	g = deepcopy(graph)
	sv = filterVertices(g,:orthoDist,v->true)
	# iterate over vertices
	for v in sv
		dists = get(props(g,v),:orthoDist,"")
		dks = keys(dists)
		for k in dks
			mems = get(get(dists,k,""),:members,"")
			for mem in mems
				# replace the accession value eg for :orthoHGNC the value is "wwp2"
				# with a dict including the accession value and the expression value
				ind =  findfirst(x->x==get(mem,key,""),names(df))
				indbc = findfirst(x->x==barcode,df[:barcode])
				if !isnothing(ind) && !isnothing(indbc)
					mem[modality] = Dict(:accession=>get(mem,key,""),
										 :barcode=>barcode,
										 :value=>df[indbc,ind])
					# println("added gene:$ind, barcode:$indbc")
				else
					# println("cannot add gene:$ind, barcode $indbc")
				end

			end
		end
	end
	g
end

# add gene expression for the given key:
# all orthoDBId's at all levels
function addExpression(graph::AbstractMetaGraph,
						dfrow::DataFrameRow,
						key::Symbol,
						modality::Symbol)
	# get vertices with orthologs
	barcode = dfrow[:barcode]
	g = deepcopy(graph)
	sv = filterVertices(g,:orthoDist,v->true)
	# iterate over vertices
	for v in sv
		dists = get(props(g,v),:orthoDist,"")
		dks = keys(dists)
		for k in dks
			mems = get(get(dists,k,""),:members,"")
			for mem in mems
				# replace the accession value eg for :orthoHGNC the value is "wwp2"
				# with a dict including the accession value and the expression value
				gnn = get(mem,key,"")
				ind = findfirst(x->x==gnn,names(dfrow))
				if !isnothing(ind)
					mem[modality] = Dict(:accession=>get(mem,key,""),
										 :barcode=>dfrow[:barcode],
										 :value=>dfrow[ind])
					# println("added $gnn = feature:$ind, barcode:$barcode")
				else
					# println("not found $gnn, barcode $barcode")
				end
			end
		end
	end
	g
end
