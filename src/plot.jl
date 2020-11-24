function plotDagLRT(targPath::NamedTuple,gParams::Dict)
    g = targPath.graph
    nodelabels = Vector{String}(undef,0)
    nodestyles = Dict{Int64,String}()
    nlab = ""
    for v in 1:nv(g)
        if any(in.(v,targPath.recinds.ind))
            println("adding rec $v")
            nlab = join(["R",v,props(g,v)[gParams[:lrt]]],":")
            push!(nodestyles,v=>"fill=green!50")
        elseif any(in.(v,targPath.liginds.ind))
            println("adding lig $v")
            nlab = join(["L",v,props(g,v)[gParams[:lrt]]],":")
            push!(nodestyles,v=>"fill=yellow!50")
        elseif any(in.(v,targPath.targinds))
            println("adding targ $v")
            nlab = join(["T",v,props(g,v)[gParams[:lrt]]],":")
            push!(nodestyles,v=>"fill=red!50")
        else
            if haskey(props(g,v),:displayName)
                println("adding p $v")
                nlab = join(["P",v,props(g,v)[:displayName]],":")
                push!(nodestyles,v=>"fill=white")
            elseif haskey(props(g,v),:displayNameIntxn)
                println("adding p $v")
                nlab = join(["P",v,props(g,v)[:displayNameIntxn]],":")
                push!(nodestyles,v=>"fill=white")
            elseif haskey(props(g,v),:ctrlRxn)
                ident = split(props(g,v)[:ctrlRxn],"/")[end]
                nlab = join(["P",v,ident],":")
            else
                println("adding p $v")
                nlab = join(["P",v],":")
                push!(nodestyles,v=>"fill=white")
            end
        end
        if length(nlab) > gParams[:nodelabelmax]
            nlab = string(nlab[1:gParams[:nodelabelmax]],"...")
            push!(nodelabels, nlab)
        else
            push!(nodelabels, nlab)
        end
        println(nlab)
    end

    gp = TikzGraphs.plot(DiGraph(g),gParams[:lt],nodelabels,
                    node_style="draw,rounded corners",
                    node_styles=nodestyles,
                    edge_style="black,line width=0.25em",
                    options=gParams[:opt])

    if split(gParams[:fname],".")[end] == "tex"
        TikzPictures.save(TEX(gParams[:fname]), gp)
    elseif split(gParams[:fname],".")[end] == "svg"
        TikzPictures.save(SVG(gParams[:fname]), gp)
    elseif split(gParams[:fname],".")[end] == "pdf"
        TikzPictures.save(PDF(gParams[:fname]), gp)
    else
        println("unknown output format: ",split(gParams[:fname],".")[end])
        throw(error())
    end
    gp
end

function plotDag(g::AbstractGraph,gParams::Dict)
    nodelabels = Vector{String}(undef,0)
    nodestyles = Dict{Int64,String}()
    # apply default styles
    for v in 1:nv(g)
        if haskey(props(g,v),:displayName)
            println("adding $v")
            nlab = join([v,props(g,v)[:displayName]],":")
            push!(nodestyles,v=>"fill=white")
        elseif haskey(props(g,v),:displayNameIntxn)
            println("adding $v")
            nlab = join([v,props(g,v)[:displayNameIntxn]],":")
            push!(nodestyles,v=>"fill=white")
        elseif haskey(props(g,v),:entId)
            println("adding $v")
            nlab = join([v,props(g,v)[:entId]],":")
            push!(nodestyles,v=>"fill=white")
        else
            println("adding $v")
            nlab = "$v"
            push!(nodestyles,v=>"fill=white")
        end
        if length(nlab) > gParams[:nodelabelmax]
            nlab = string(nlab[1:gParams[:nodelabelmax]],"...")
            push!(nodelabels, nlab)
        else
            push!(nodelabels, nlab)
        end
    end

    # process custom key styles and filter
    if haskey(gParams,:keystyles)
        for s in 1:length(gParams[:keystyles])
            sty = collect(gParams[:keystyles])[s]
            for v in 1:nv(g)
                if haskey(props(g,v),sty[1])
                    nodestyles[v] = sty[2]
                end
            end
        end
    end

    if haskey(gParams,:filterstyles)
        for s in 1:length(gParams[:filterstyles])
            sty = collect(gParams[:filterstyles][s])
            verts = filterVertices(g,sty[1][1],sty[1][2][2])
            for v in verts
                nodestyles[v] = sty[1][2][1]
            end
        end
    end

    gp = TikzGraphs.plot(DiGraph(g),gParams[:lt],nodelabels,
                    node_style="draw,rounded corners",
                    node_styles=nodestyles,
                    edge_style="black,line width=0.25em",
                    options=gParams[:opt])

    if split(gParams[:fname],".")[end] == "tex"
        TikzPictures.save(TEX(gParams[:fname]), gp)
    elseif split(gParams[:fname],".")[end] == "svg"
        TikzPictures.save(SVG(gParams[:fname]), gp)
    elseif split(gParams[:fname],".")[end] == "pdf"
        TikzPictures.save(PDF(gParams[:fname]), gp)
    else
        println("unknown output format: ",split(gParams[:fname],".")[end])
        throw(error())
    end
    gp
end

# create a colorscale for the plot based on max value orthologs across all nodes
function plotDagExp(g::AbstractGraph,gParams::Dict)
    nodelabels = Vector{String}(undef,0)
    nodequants = Vector{Union{Missing,Float64}}(undef,0)
    nodestyles = Dict{Int64,String}()
    if haskey(gParams,:colorscheme)
        cscheme = gParams[:colorscheme]
    else
        cscheme = ColorSchemes.plasma
    end


    # get the max value of all expressed orthologs at each node
    for v in 1:nv(g)
        quantmax = missing
        if haskey(props(g,v),:orthoDist)
            vorth = props(g,v)[:orthoDist][gParams[:orthodist]][:members]
            if length(vorth) > 0
                println("length of vorth is ",length(vorth))
                # for o in 1:length(vorth)
                #     if haskey(vorth[o],:rnaSeq)
                #     end
                # end
                if length(collect(skipmissing([haskey(vorth[o],:rnaSeq) ? vorth[o][:rnaSeq][:value] : missing for o in 1:length(vorth)]))) > 0
                    quantmax = reduce(max,skipmissing([haskey(vorth[o],:rnaSeq) ? vorth[o][:rnaSeq][:value] : missing for o in 1:length(vorth)]))
                    println(quantmax)
                else
                    println("no expression of ortholog")
                end
            else
                println("vorth is empty")
            end
        else
            println("orthoDist missing")
        end
        push!(nodequants,quantmax)
    end
    # derive a log colorscale for the expression values
    scaledquants = copy(nodequants)
    if gParams[:scale] == :log
        logquants = log.(collect(skipmissing(nodequants)).+1)
        α = 1/maximum(logquants)
        for v in 1:nv(g)
            if !ismissing(scaledquants[v])
                scaledquants[v] = log(scaledquants[v])*α
            end
        end
    elseif gParams[:scale] == :ident
        α = 1/maximum(collect(skipmissing(scaledquants)))
        for v in 1:nv(g)
            if !ismissing(scaledquants[v])
                scaledquants[v] = scaledquants[v]*α
            end
        end
    end


    for v in 1:nv(g)
        # set node styles as expression, when possible
        if !ismissing(scaledquants[v])
            scaled = get(cscheme,scaledquants[v])
            rgbstr = string("fill={rgb:red,",scaled.r,
                                 ";green,",scaled.g,
                                 ";blue,",scaled.b,"}")
            push!(nodestyles,v=>rgbstr)
        else
            push!(nodestyles,v=>"fill=white")
        end

        # set node labels
        if haskey(props(g,v),:displayName)
            println("adding $v")
            nlab = join([v,props(g,v)[:displayName]],":")
        elseif haskey(props(g,v),:displayNameIntxn)
            println("adding $v")
            nlab = join([v,props(g,v)[:displayNameIntxn]],":")
        elseif haskey(props(g,v),:entId)
            println("adding $v")
            nlab = join([v,props(g,v)[:entId]],":")
        else
            println("adding $v")
            nlab = "$v"
        end
        if length(nlab) > gParams[:nodelabelmax]
            nlab = string(nlab[1:gParams[:nodelabelmax]],"...")
            push!(nodelabels, nlab)
        else
            push!(nodelabels, nlab)
        end
    end

    # process custom key styles and filter
    if haskey(gParams,:keystyles)
        for s in 1:length(gParams[:keystyles])
            sty = collect(gParams[:keystyles])[s]
            for v in 1:nv(g)
                if haskey(props(g,v),sty[1])
                    nodestyles[v] = sty[2]
                end
            end
        end
    end

    if haskey(gParams,:filterstyles)
        for s in 1:length(gParams[:filterstyles])
            sty = collect(gParams[:filterstyles][s])
            verts = filterVertices(g,sty[1][1],sty[1][2][2])
            for v in verts
                nodestyles[v] = sty[1][2][1]
            end
        end
    end

    gp = TikzGraphs.plot(DiGraph(g),gParams[:lt],nodelabels,
                    node_style="draw,rounded corners",
                    node_styles=nodestyles,
                    edge_style="black,line width=0.25em",
                    options=gParams[:opt])

    if split(gParams[:fname],".")[end] == "tex"
        TikzPictures.save(TEX(gParams[:fname]), gp)
    elseif split(gParams[:fname],".")[end] == "svg"
        TikzPictures.save(SVG(gParams[:fname]), gp)
    elseif split(gParams[:fname],".")[end] == "pdf"
        TikzPictures.save(PDF(gParams[:fname]), gp)
    else
        println("unknown output format: ",split(gParams[:fname],".")[end])
        throw(error())
    end
    gp
end
