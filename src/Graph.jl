function getEdges(edgelist,edgelabels)
    list = Int.(zeros(length(edgelist)))
    for n in 1:length(edgelabels)
        list[findall(x->x==edgelabels[n],edgelist)] .= n
    end
    list
end

function loadGraph(left,right,allnodes)
    g = SimpleDiGraph(length(allnodes))
    for i in 1:length(left)
        add_edge!(g, left[i], right[i]);
    end
    g
end
