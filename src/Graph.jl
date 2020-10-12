

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
	symOut = "http://www.biopax.org/release/biopax-level3.owl#left"
	symIn = "http://www.biopax.org/release/biopax-level3.owl#right"
	for i in 1:len

		p_nested = getNested(nestedParams,df[i,:participant])
		if hasCol(p_nested,nestedParams[:simpleEntity][1])
			pTerms = vcat(map(k->nestedParams[k],dataKeys[[3,4]])...)
		else
			pTerms =vcat(map(k->nestedParams[k],dataKeys[[3]])...)
		end
		for ii in 1:size(p_nested,1)
			if !ismissing(df[i,nestedParams[:ctrlInteraction][1]])
				baseTerms = vcat(map(k->nestedParams[k],dataKeys[[1,2]])...)
				c_nested = getNested(nestedParams,df[i,:ctrlEntity])
				if hasCol(c_nested,nestedParams[:simpleEntity][1])
					cTerms =vcat(map(k->nestedParams[k],dataKeys[[3,4]])...)
					cColNames =vcat(map(k->nestedParams[k],ctrlKeys[[1,2]])...)
				else
					cTerms =vcat(map(k->nestedParams[k],dataKeys[[3]])...)
					cColNames =vcat(map(k->nestedParams[k],ctrlKeys[[1]])...)
				end
				for iii in 1:size(c_nested,1)
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
		p_ind = findfirst(n->n==int_df[i,:participantRef],allverts)
		set_props!(g, p_ind, Dict(:unification=>int_df[i,:participantRef],
								  :entityType=>int_df[i,:participantType],
								  :location=>int_df[i,:participantLocRef]))

		rx_ind = findfirst(n->n==int_df[i,:interaction],allverts)
		set_props!(g,rx_ind,Dict(:unification=>int_df[i,:interaction],
								  :entityType=>int_df[i,:intType]))

		if int_df[i,:partPred] == symIn
			add_edge!(g,p_ind,rx_ind)
		elseif int_df[i,:partPred] == symOut
			add_edge!(g,rx_ind,p_ind)
		else
			throw("only left/right reactions are supported")
		end

		# optional ctrlRef vertex and edge properties
		if !ismissing(int_df[i,:ctrlRxn])
			ct_ind = findfirst(n->n==int_df[i,:ctrlEntityRef],allverts)
			ct_rx_ind = findfirst(n->n==int_df[i,:ctrlRxn],allverts)
			set_props!(g, ct_ind, Dict(:unification=>int_df[i,:ctrlEntityRef],
									   :entityType=>int_df[i,:ctrlEntityType],
									   :location=>int_df[i,:ctrlEntityLocRef]))

			set_props!(g, ct_rx_ind, Dict(:unification=>int_df[i,:ctrlRxn],
									   	  :entityType=>int_df[i,:ctrlEntityRef],
										  :controlType=>int_df[i,:ctrlRxnType]))

			add_edge!(g,ct_ind,ct_rx_ind)
			add_edge!(g,ct_rx_ind,rx_ind)
		end
	end
	Dict(:graph=>g,:vertices=>allverts,:simple=>int_df)
end
