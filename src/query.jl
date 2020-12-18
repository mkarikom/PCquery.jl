function delimitValues(iris,domain::String,wrap;prependhash=false)
	# some pathway commons rdf xml dumps specify "#" after base uri, add this manually to pre-stored or GET requested URI from pathway commons site, where there are no "#" in the URI
	if prependhash
		iris = map(x->joinpath(split(x,"/")[1:end-1]...,string("#",split(x,"/")[end])),iris)
	end
	if length(wrap) > 0
		xs = string.(wrap[1],domain,iris,wrap[2])
	else
		xs = string.(domain,iris)
	end
	if typeof(iris) <: Vector
		xs = join(xs," ")
		xs = "{$xs}"
		xs = "$xs"
	else
		xs = "{$xs}"
		xs = "$xs"
	end
end

function delimitValues(iris::Vector,domain::String;prependhash=false)
	if prependhash
		iris = map(x->joinpath(split(x,"/")[1:end-1]...,string("#",split(x,"/")[end])),iris)
	end
	xs = string.(domain,iris)
	xs = join(xs," ")
	xs = "{$xs}"
	xs = "$xs"
end

function urlForm(domain::String,path::String,value::String)
	xs = "<$domain$path$value>"
end

function summarizeVals(df::DataFrame, scol::Symbol, fcol::Symbol)
	svals = unique(df[!,scol])
	featvals = []
	for i in 1:length(svals)
	      vals = []
	      for j in 1:size(df,1)
			  if df[!,scol][j] == svals[i]
				  push!(vals,df[!,fcol][j])
			  end
	      end
	      push!(featvals,unique(vals))
	end
	DataFrame(scol=>svals,fcol=>featvals)
end

function requestTTLold(port::Int64,fmt::String,host::String,path::String,
					accept::String,qstr::String)
	scheme = "http"
	port = string(port)
	fmt = fmt
	#host = "nie2.math.uci.edu"
	host = host
	path = path
	turtle = qstr

	query = Dict("timeout"=>"0",
	             "query"=>turtle,
	             "debug"=>"on",
	             "format"=>fmt)
	uri = HTTP.URI(scheme=scheme,host=host,port=port,path=path,query=query)
	resp = HTTP.request(:GET,uri,["Accept" => accept]);
	#resp = HTTP.request(:GET,uri);
	resp_str = String(resp.body);
	xdoc = LightXML.parse_string(resp_str)
	xroot = LightXML.root(xdoc)
	results = PCquery.getXMLres(xroot)
end


function requestTTL(port::Int64,
					scheme::String,
					fmt::String,
					host::String,
					path::String,
					method::Symbol,
					header::Array{Pair{String,String},1},
					qstr::String)
	port = string(port)
	query = Dict("query"=>qstr)
	uri = HTTP.URI(scheme=scheme,host=host,port=port,path=path)
	resp = HTTP.request(method,uri,header,HTTP.escapeuri(query));
end

function requestTTLget(scheme::String,fmt::String,host::String,path::String,
					accept::String,qstr::String)
	query = Dict("timeout"=>"0",
	             "query"=>qstr,
	             "debug"=>"on",
	             "format"=>fmt)
	uri = HTTP.URI(scheme=scheme,host=host,path=path,query=query)
	resp = HTTP.request(:GET,uri,["Accept" => accept]);
	#resp = HTTP.request(:GET,uri);
	resp_str = String(resp.body);
	xdoc = LightXML.parse_string(resp_str)
	xroot = LightXML.root(xdoc)
	results = PCquery.getXMLres(xroot)
end


function requestTTL(fmt::String,
					scheme::String,
					host::String,
					path::String,
					method::Symbol,
					header::Array{Pair{String,String},1},
					qstr::String)
	fmt = fmt
	#host = "nie2.math.uci.edu"
	host = host
	path = path
	turtle = qstr

	query = Dict("query"=>turtle)

	uri = HTTP.URI(scheme=scheme,host=host,path=path)
	# resp = HTTP.request(:GET,uri,["Accept" => accept]);
	resp = HTTP.request(method,uri,header,HTTP.escapeuri(query));
	#resp = HTTP.request(:GET,uri);
end

function parseSparqlResponse(resp::HTTP.Messages.Response)
	resp_str = String(resp.body);
	xdoc = LightXML.parse_string(resp_str)
	xroot = LightXML.root(xdoc)
	results = PCquery.getXMLres(xroot)
end

function parseSparqlResponse!(resp::HTTP.Messages.Response,df)
	resp_str = String(resp.body);
	xdoc = LightXML.parse_string(resp_str)
	xroot = LightXML.root(xdoc)
	results = PCquery.getXMLres!(xroot,df)
end



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

function getNames(node;verbose=false)
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
        verbose ? println("$v = $cont") : nothing
    end end end end end end end
    unique(dkeys)
end

function getXMLres(node;verbose=false)
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
            verbose ? println(counter) : nothing
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
                        verbose ? println("$v = $cont") : nothing
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

function getXMLres!(node,df;verbose=false)
    dkeys = getNames(node)
    colnames = Symbol.(dkeys)
    counter = 1
    for c in child_nodes(node)  # c is an instance of XMLNode
    if is_elementnode(c)
    if name(XMLElement(c)) == "results"
    for res in child_nodes(c)
        counter += 1
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
                    verbose ? println("$v = $cont") : nothing
                end
            end
        end
        if length(node_data) > 0
            push!(df,initRow(colnames,Symbol.(node_names),tuple(node_data...)))
        end
    end end end end
end

# convenience get request to pc top pathways
function topPathsGet(qParams::Dict)
	rParams = Dict(
			:scheme=>"https",
			:host=>"www.pathwaycommons.org",
			:path=>"/pc2/top_pathways.json",
			:query=>qParams)
	uri = URIs.URI(
			scheme=rParams[:scheme],
			host=rParams[:host],
			path=rParams[:path],
			query=rParams[:query])
	resp = HTTP.request("GET",uri)
	parsed = JSON.parse(String(resp.body))
	hits = parsed["searchHit"]
	df = DataFrame(
			typeof.(values(hits[1])),
			Symbol.(names(hits[1])),0)
	for i in 1:length(hits)
		push!(df,hits[i])
	end
	unique(df)
end

function searchPathsGet(params::Dict)
	qParams = copy(params)
	namefilter = nothing
	if haskey(qParams,:filter)
		namefilter = pop!(qParams,:filter)
	end
	rParams = Dict(
			:scheme=>"https",
			:host=>"www.pathwaycommons.org",
			:path=>"/pc2/search.json",
			:query=>qParams)
	uri = URIs.URI(
			scheme=rParams[:scheme],
			host=rParams[:host],
			path=rParams[:path],
			query=rParams[:query])
	resp = HTTP.request("GET",uri)
	parsed = JSON.parse(String(resp.body))
	hits = parsed["searchHit"]
	df = DataFrame(
			typeof.(values(hits[1])),
			Symbol.(names(hits[1])),0)
	if !isnothing(namefilter)
		for i in 1:length(hits)
			for j in 1:length(namefilter)
				if all(map(t->occursin(lowercase(t),lowercase.(hits[i]["name"])),namefilter[j]))
					push!(df,hits[i])
				end
			end
		end
	else
		for i in 1:length(hits)
			push!(df,hits[i])
		end
	end
	unique(df)
end
