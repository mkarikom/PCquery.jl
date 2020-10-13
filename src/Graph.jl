

# initialize a bipartite graph where {left,right,control} ʌ {interactions} = ∅
# for edge [1,2], direction is 1->2.  define: [left,interaction], [interaction,right], [ctrl,interaction]
function initBpGraph(df::DataFrame,nestedParams,dataKeys,ctrlKeys)
	len = size(df,1)
	# collect verts
	allKeys = vcat(dataKeys,ctrlKeys)
	colnames = [:ind,vcat(map(k->nestedParams[k],allKeys)...)...]
	coltypes = fill(Union{Missing,String},length(colnames))
	coltypes[1] = Union{Missing,Int64}
	int_df = DataFrame(coltypes,colnames)

	# criteria for participant edge direction
	symOut = "http://www.biopax.org/release/biopax-level3.owl#right"
	symIn = "http://www.biopax.org/release/biopax-level3.owl#left"
	for i in 1:len

		p_nested = getNested(nestedParams,df[i,:participant])
		for ii in 1:size(p_nested,1)
			if hasCol(p_nested,nestedParams[:simpleEntity][1])
				pTerms = vcat(map(k->nestedParams[k],dataKeys[[3,4]])...)
			else
				pTerms =vcat(map(k->nestedParams[k],dataKeys[[3]])...)
			end
			if !ismissing(df[i,nestedParams[:ctrlInteraction][1]])
				baseTerms = vcat(map(k->nestedParams[k],dataKeys[[1,2]])...)
				c_nested = getNested(nestedParams,df[i,:ctrlEntity])
				for iii in 1:size(c_nested,1)
					if hasCol(c_nested,nestedParams[:simpleEntity][1])
						cTerms =vcat(map(k->nestedParams[k],dataKeys[[3,4]])...)
						cColNames =vcat(map(k->nestedParams[k],ctrlKeys[[1,2]])...)
					else
						cTerms =vcat(map(k->nestedParams[k],dataKeys[[3]])...)
						cColNames =vcat(map(k->nestedParams[k],ctrlKeys[[1]])...)
					end
					data = tuple(i,collect(df[i,baseTerms])...,
									collect(p_nested[ii,pTerms])...,
									collect(c_nested[iii,cTerms])...)
					push!(int_df,initRow(
							colnames,
							[:ind,vcat(baseTerms,pTerms,cColNames)...],
							data))
				end
			else
				baseTerms = vcat(map(k->nestedParams[k],dataKeys[[1]])...)
				data = tuple(i,collect(df[i,baseTerms])...,
								collect(p_nested[ii,pTerms])...)
				push!(int_df,initRow(
						colnames,
						[:ind,vcat(baseTerms,pTerms)...],
						data))
			end
		end
	end
	allverts = sortUnique(int_df[!,:participantRef],int_df[!,:interaction],
							int_df[!,:ctrlEntityRef],int_df[!,:ctrlRxn])
	g = MetaDiGraph(length(allverts))

	for i in 1:size(int_df,1)
		# participants
		ptype = split(int_df[i,:participantType],"#")[2]
		rxtype = split(int_df[i,:intType],"#")[2]
		dxn = int_df[i,:partPred][end-4:end]

		p_ind = findfirst(n->n==int_df[i,:participantRef],allverts)
		if !ismissing(int_df[i,nestedParams[:simpleEntity][1]])
			set_props!(g,p_ind,
				Dict(nestedParams[:simpleEntity] .=>
					collect(int_df[i,nestedParams[:simpleEntity]])))
		end
		set_props!(g, p_ind, Dict(:unification=>int_df[i,:participantRef],
								  :entityType=>int_df[i,:participantType],
								  :location=>int_df[i,:participantLocRef]))

		rx_ind = findfirst(n->n==int_df[i,:interaction],allverts)
		set_props!(g,rx_ind,Dict(:unification=>int_df[i,:interaction],
								  :entityType=>int_df[i,:intType]))


	    cond = int_df[i,:participantType] == "http://www.biopax.org/release/biopax-level3.owl#Dna"

		if int_df[i,:partPred] == symIn
			hasedge = has_edge(g, p_ind, rx_ind)
			if hasedge
				# println("\n\nedge $p_ind -> $rx_ind already exists")
				# debugedge(int_df,i,g,p_ind,rx_ind,cond)
			else
				added = add_edge!(g,p_ind,rx_ind)
			end
			# added = add_edge!(g,p_ind,rx_ind)
			# @assert added "edge $p_ind -> $rx_ind was not added"
		elseif int_df[i,:partPred] == symOut
			hasedge = has_edge(g, rx_ind, p_ind)
			if hasedge
				# println("\n\nedge $rx_ind -> $p_ind  already exists")
				# debugedge(int_df,i,g,p_ind,rx_ind,cond)
			else
				added = add_edge!(g,rx_ind,p_ind)
			end
			# added = add_edge!(g,rx_ind,p_ind)
			# @assert added "edge $rx_ind -> $p_ind was not added"
		else
			throw("only left/right reactions are supported")
		end

		# if int_df[i,:participantRef] == "http://pathwaycommons.org/pc11/#UnificationXref_reactome_R-HSA-2127254"
		# 	println("found dna rxn")
		# 	debugedge(int_df,i,g,p_ind,rx_ind)
		# end


		# optional ctrlRef vertex and edge properties
		if !ismissing(int_df[i,:ctrlRxn])
			ct_ind = findfirst(n->n==int_df[i,:ctrlEntityRef],allverts)
			ct_rx_ind = findfirst(n->n==int_df[i,:ctrlRxn],allverts)
			if !ismissing(int_df[i,nestedParams[:ctrlSimpleEntity][1]])
				set_props!(g,ct_ind,
					Dict(nestedParams[:ctrlSimpleEntity] .=>
						collect(int_df[i,nestedParams[:ctrlSimpleEntity]])))
			end
			set_props!(g, ct_ind, Dict(:unification=>int_df[i,:ctrlEntityRef],
									   :entityType=>int_df[i,:ctrlEntityType],
									   :location=>int_df[i,:ctrlEntityLocRef]))

			set_props!(g, ct_rx_ind, Dict(:unification=>int_df[i,:ctrlRxn],
									   	  :entityType=>int_df[i,:ctrlRxnType],
										  :controlType=>int_df[i,:ctrlRxnDir]))

			add_edge!(g,ct_ind,ct_rx_ind)
			add_edge!(g,ct_rx_ind,rx_ind)
		end
	end
	Dict(:graph=>g,:vertices=>allverts,:simple=>int_df)
end

function debugedge(int_df,i,g,p_ind,rx_ind,cond)
	if cond
	rxtype = int_df[i,:intType]
	ptype = int_df[i,:participantType]
	oldRxType = props(g,rx_ind)[:entityType]
	oldPType = props(g,p_ind)[:entityType]
	rxref = int_df[i,:interaction]
	pref = int_df[i,:participantRef]
	rdir = int_df[i,:partPred]
	pinn = inneighbors(g,p_ind)
	pout = outneighbors(g,p_ind)
	rinn = inneighbors(g,rx_ind)
	rout = outneighbors(g,rx_ind)
	println("\nparticipant direction is $rdir")
	println("participant $p_ind ref is $pref")
	println("rx $rx_ind ref is $rxref")
	println("participant $p_ind in neighbors are $pinn")
	println("participant $p_ind out neighbors are $pout")
	println("rx $rx_ind in neighbors are $rinn")
	println("rx $rx_ind out neighbors are $rout")
	println("attempting to change p from $oldPType to $ptype")
	println("attempting to change rx from $oldRxType to $rxtype\n\n")
end
end
