function subStrVec(vals::Vector,pos::Vector)
    output = []
    for i in vals
        push!(output,i[pos...])
    end
end
