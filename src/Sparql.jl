
function delimitValues(iris::Vector{String},varname::String)
	xs = join(iris," ")
	xs = "{$xs}"
	xs = "$varname $xs"
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
					fmt::String,
					host::String,
					path::String,
					method::Symbol,
					header::Array{Pair{String,String},1},
					qstr::String)
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
	uri = HTTP.URI(scheme=scheme,host=host,port=port,path=path)
	# resp = HTTP.request(:GET,uri,["Accept" => accept]);
	resp = HTTP.request(method,uri,header,HTTP.URI(query));
	#resp = HTTP.request(:GET,uri);
	resp_str = String(resp.body);
	xdoc = LightXML.parse_string(resp_str)
	xroot = LightXML.root(xdoc)
	results = PCquery.getXMLres(xroot)
end
