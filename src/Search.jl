# get a vector of Union{Missing,Type} for the df constructor so that rows with missing values can be added
function makeUnions(types::Vector{DataType})
    tvec = Union[]
    for t in 1:length(types)
        push!(tvec,Union{Missing,types[t]})
    end
    tvec
end

# get the row-dictionary to add to the df
function initRow(allkeys::Vector{Symbol},dkeys::Vector{Symbol},data::Vector)
    tuples = []
    for k in 1:length(allkeys)
        ind = findfirst(x->x==allkeys[k],dkeys)
        if isnothing(ind)
            push!(tuples,(allkeys[k],missing))
        else
            push!(tuples,(allkeys[k],data[ind]))
        end
    end
    Dict(tuples)
end


function resultData(xmlRoot::XMLElement)
    ces = collect(child_elements(xmlRoot));
    # get the column names
    dkeys = []
    dtypes = []
    for c in ces
        hit = collect(child_elements(c))
        dkey = name.(hit)
        dtype = typeof.(content.(hit))
        push!(dkeys,dkey)
        push!(dtypes,dtype)
    end
    dkeys = reduce(vcat,dkeys)
    dtypes = reduce(vcat,dtypes)
    ind = [findfirst(x->x==k,dkeys) for k in unique(dkeys)]

    colnames = Symbol.(dkeys[ind])
    coltypes = makeUnions(dtypes[ind])
    df = DataFrame(coltypes,colnames)
    for c in ces
        hit = collect(child_elements(c))
        keys = Symbol.(name.(hit))
        data = content.(hit)
        # tdata = initRow(colnames,keys,data)
        tdata = initRow(colnames,keys,data)
        push!(df,tdata)
    end
    df
end

function getNames(node)
    dkeys = []
    for c in child_nodes(node)  # c is an instance of XMLNode
    if is_elementnode(c)
    if name(XMLElement(c)) == "results"
    for res in child_nodes(c)
    for cres in child_nodes(res)  # c is an instance of XMLNode
    if is_elementnode(cres)
    e = XMLElement(cres)  # this makes an XMLElement instance
    for a in attributes(e)  # a is an instance of XMLAttr
        n = name(a)
        v = value(a)
        cont = content(e)
        push!(dkeys,v)
        # println("$v = $cont")
    end end end end end end end
    unique(dkeys)
end

function getXMLres(node)
    dkeys = getNames(node)
    colnames = Symbol.(dkeys)
    coltypes = fill(Union{Missing,String},length(colnames))
    df = DataFrame(coltypes,colnames)
    counter = 1
    for c in child_nodes(node)  # c is an instance of XMLNode
    if is_elementnode(c)
    if name(XMLElement(c)) == "results"
    for res in child_nodes(c)
        counter += 1
        # println(counter)
        node_data = []
        node_names = []
        for cres in child_nodes(res)  # c is an instance of XMLNode
            if is_elementnode(cres)
                e = XMLElement(cres)  # this makes an XMLElement instance
                for a in attributes(e)  # a is an instance of XMLAttr
                    n = name(a)
                    v = value(a)
                    cont = content(e)
                    push!(node_data,Base.strip(cont))
                    push!(node_names,Base.strip(v))
                    # println("$v = $cont")
                end
            end
        end
        if length(node_data) > 0
            push!(df,initRow(colnames,Symbol.(node_names),node_data))
        end
    end end end end
    df
end

function indexFeature(flist::Vector,fname::Symbol,
                      df::DataFrame,offset::Vector)
    inds = []
    for r in 1:size(df,1)
        len = length(df[!,fname][r])
        urng = len-offset[1]:len-offset[2]
        ind = findall(x->x==df[!,fname][r][urng],flist)[]
        push!(inds,ind)
    end
    inds
end

# get left,right using the value of surveilance rows to identify:
# each row will always have a "left" and a "right" (for the head and tail of a pathway these carry special importance)
# the head of the pathway as a "left" where the surveilance column is "lig_rec_complex"
# the tail of the pathway as a "right" where the surveilance column is "transcription complex"
function nodeLabels(df::DataFrame,left::Symbol,right::Symbol,
					keyvals...)
	# Ex.
	# keyvals[1] = Dict(:lr=>:left,
	#					:key=>:llocref,
	# 				    :val=>"http://pathwaycommons.org/pc11/#UnificationXref_gene_ontology_GO_0005829")
	len = size(df,1)
	nodes = vcat(df[!,left],df[!,right])
	vals = []
	lr = []

	for n in 1:length(nodes)
		for kv in 1:length(keyvals)
			if keyvals[kv][:lr] == lr[n]
				push!(lr,lr[n])
				push!(vals,keyvals[kv][:val])
			end
		end

	end
	DataFrame(:nodes=>nodes,:vals=>vals,:lr=>lr)

	#
	# allkeys = vcat(df[!,leftkey],df[!,rightkey])
	# lr = vcat(fill(leftkey,len),fill(rightkey,len))
	# uniqueind = unique(i -> allkeys[i], 1:length(allkeys))
	#
	# allkeys = allkeys[uniqueind]
	# lr = lr[uniqueind]
	#
	# sind = sortperm(allkeys)
	# allkeys = allkeys[sind]
	# lr = lr[sind]
	#
	# dflr = []
	# for k in 1:length(allkeys)
	# 	if lr[k] == leftkey
	# 		findfirst(x->x==leftval,df[!,leftkey])
	# 	end
	# 	findfirst(x->, df[!,lr[k]])
	# end
	# unique(i -> x[i], 1:length(x))
	#
	# labels = Int64(zeros(length(nodes)))
	#
	# for v in 1:length(featvals)
	# 	for n in 1:length(nodes)
	# 		findfirst
	# 	end
	# end
end
