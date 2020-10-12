

# get the row-dictionary to add to the df
function initRow(allkeys::Vector{Symbol},dkeys::Vector{Symbol},data::Tuple)
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
    if length(dkeys) > 0
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
                push!(df,initRow(colnames,Symbol.(node_names),tuple(node_data...)))
            end
        end end end end
    else
        df = Array{Any,1}()
    end
    df
end

function getXMLres!(node,df)
    dkeys = getNames(node)
    colnames = Symbol.(dkeys)
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
            push!(df,initRow(colnames,Symbol.(node_names),tuple(node_data...)))
        end
    end end end end
end
