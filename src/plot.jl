function plotDagLRT(targPath::NamedTuple,gParams::Dict)
    g = targPath.graph
    nodelabels = Vector{String}(undef,0)
    nodestyles = Dict{Int64,String}()
    for v in 1:nv(g)
        if any(in.(v,targPath.recinds.ind))
            println("adding rec $v")
            push!(nodelabels, join(["R",v,props(g,v)[gParams[:lrt]]],":"))
            push!(nodestyles,v=>"fill=green!50")
        elseif any(in.(v,targPath.liginds.ind))
            println("adding lig $v")
            push!(nodelabels, join(["L",v,props(g,v)[gParams[:lrt]]],":"))
            push!(nodestyles,v=>"fill=yellow!50")
        elseif any(in.(v,targPath.targinds))
            println("adding targ $v")
            push!(nodelabels, join(["T",v,props(g,v)[gParams[:lrt]]],":"))
            push!(nodestyles,v=>"fill=red!50")
        else
            if haskey(props(g,v),:displayName)
                println("adding p $v")
                # push!(nodelabels,string(v))
                push!(nodelabels, join(["P",v,props(g,v)[:displayName]],":"))
                push!(nodestyles,v=>"fill=white")
            elseif haskey(props(g,v),:displayNameIntxn)
                println("adding p $v")
                # push!(nodelabels,string(v))
                push!(nodelabels, join(["P",v,props(g,v)[:displayNameIntxn]],":"))
                push!(nodestyles,v=>"fill=white")
            elseif haskey(props(g,v),:ctrlRxn)
                ident = split(props(g,v)[:ctrlRxn],"/")[end]
                push!(nodelabels, join(["P",v,ident],":"))
            else
                println("adding p $v")
                push!(nodelabels,string(v))
                # push!(nodelabels, join(["P",v,props(g,v)[pw]],":"))
                push!(nodestyles,v=>"fill=white")
            end
        end
    end

    gp = TikzGraphs.plot(DiGraph(g),gParams[:lt],nodelabels,
                    node_style="draw,rounded corners",
                    node_styles=nodestyles,
                    edge_style="black,line width=0.25em",
                    options="scale=2, font=\\huge\\sf")

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

function plotDag(targPath::NamedTuple,gParams::Dict)
    g = targPath.graph
    nodelabels = Vector{String}(undef,0)
    nodestyles = Dict{Int64,String}()
    for v in 1:nv(g)
        if haskey(props(g,v),:displayName)
            println("adding $v")
            push!(nodelabels, join([v,props(g,v)[:displayName]],":"))
            push!(nodestyles,v=>"fill=white")
        elseif haskey(props(g,v),:displayNameIntxn)
            println("adding $v")
            push!(nodelabels, join([v,props(g,v)[:displayNameIntxn]],":"))
            push!(nodestyles,v=>"fill=white")
        elseif haskey(props(g,v),:entId)
            println("adding $v")
            push!(nodelabels, join([v,props(g,v)[:entId]],":"))
            push!(nodestyles,v=>"fill=white")
        else
            println("adding $v")
            push!(nodelabels,"$v")
            push!(nodestyles,v=>"fill=white")
        end
    end

    gp = TikzGraphs.plot(DiGraph(g),gParams[:lt],nodelabels,
                    node_style="draw,rounded corners",
                    node_styles=nodestyles,
                    edge_style="black,line width=0.25em",
                    options="scale=2, font=\\huge\\sf")

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
