function searchPathways(dbParams::Dict;verbose=false)
	if dbParams[:pnbool]==|
		pnamefilter = lowercase(string("\"",join(string.("(?=.*",dbParams[:pnames],")"),"|"),"\""))
	elseif dbParams[:pnbool] == &
		pnamefilter = lowercase(string("\"",join(string.("(?=.*",dbParams[:pnames],")")),"\""))
	end
	if dbParams[:pdsbool]==|
		pdsnamefilter = lowercase(string("\"",join(string.("(?=.*",dbParams[:pdsnames],")"),"|"),"\""))
	elseif dbParams[:pdsbool] == &
		pdsnamefilter = lowercase(string("\"",join(string.("(?=.*",dbParams[:pdsnames],")")),"\""))
	end
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/PathwayCommons/rq")
    str = open(f->read(f, String), string(rqDir,"/",dbParams[:pwrqfile]));
	cnames = [:pw,:pdn,:ds,:dsdn]
	ctypes = [fill(Vector{Union{Missing,String}},4)...]
	df_pairs = DataFrame(cnames .=> [type[] for type in ctypes])
	lengthdf = 10000
	offset = 0
	while lengthdf > 0
		turtle = Mustache.render(str,Dict{Any,Any}("pnamefilter"=>pnamefilter,
		     									   "pdsnamefilter"=>pdsnamefilter,
												   "resultLim"=>dbParams[:resultLim],
												   "resultOffset"=>string(offset)))

	    # compose the header and execute query
	    header = ["Content-Type" => "application/x-www-form-urlencoded",
	    		  "Accept" => "application/sparql-results+xml"]
		# add default graph uri if provided
		haskey(dbParams,:defaultgraphuri) ? push!(header,"default-graph-uri"=>dbParams[:defaultgraphuri]) : nothing
		# local dump of
		resp = requestTTL(dbParams[:port],dbParams[:protocol],dbParams[:format],
	                        dbParams[:host],dbParams[:path],dbParams[:method],
	                        header,turtle)
		new_pairs = parseSparqlResponse(resp)
		lengthdf = size(new_pairs)[1]
		if lengthdf > 0
			df_pairs = vcat(df_pairs,new_pairs;cols=:union)
			oldoffset = offset
			offset = offset + size(new_pairs)[1]
			verbose ? println("processed results $oldoffset-$offset") : false
		end
	end
	df_pairs
end


function getPathways(dbParams::Dict;verbose=false)
	if haskey(dbParams,:prependhash) && dbParams[:prependhash]
		val = delimitValues(dbParams[:refs],"","<>";prependhash=dbParams[:prependhash])
	else
		val = delimitValues(dbParams[:refs],"","<>")
	end
    # val = delimitValues(dbParams[:refs],"","<>")
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/PathwayCommons/rq")
    str = open(f->read(f, String), string(rqDir,"/",dbParams[:pwrqfile]));
	cnames = [:pw,:displayNamePw,:interaction,:intDisplayName,:intType,:partPred,:participant,:particiapantType,:ctrlRxn,:ctrlRxnDisplayName,:ctrlRxnType,:ctrlRxnDir,:ctrlEntity,:ctrlEntityType,:ctrlEntityLocRef,:intPubTitles,:intPubDbs,:intPubIds]
	ctypes = [fill(Vector{Union{Missing,String}},15)...,fill(Vector{Union{Missing,Vector{String}}},3)...]
	df_pairs = DataFrame(cnames .=> [type[] for type in ctypes])
	lengthdf = 10000
	offset = 0
	while lengthdf > 0
		turtle = Mustache.render(str,Dict{Any,Any}("pathVals"=>val,
		     									   "resultLim"=>dbParams[:resultLim],
												   "resultOffset"=>string(offset),
												   "fromgraph"=>dbParams[:fromgraph]))

	    # compose the header and execute query
	    header = ["Content-Type" => "application/x-www-form-urlencoded",
	    		  "Accept" => "application/sparql-results+xml"]
		# add default graph uri if provided
		# haskey(dbParams,:defaultgraphuri) ? header["default-graph-uri"] = dbParams[:defaultgraphuri] : nothing # old style
		haskey(dbParams,:defaultgraphuri) ? push!(header,"default-graph-uri"=>dbParams[:defaultgraphuri]) : nothing

		# local dump of
		resp = requestTTL(dbParams[:port],dbParams[:protocol],dbParams[:format],
	                        dbParams[:host],dbParams[:path],dbParams[:method],
	                        header,turtle)
		new_pairs = parseSparqlResponse(resp)
		lengthdf = size(new_pairs)[1]
		if lengthdf > 0
			df_pairs = vcat(df_pairs,new_pairs;cols=:union)
			oldoffset = offset
			offset = offset + size(new_pairs)[1]
			verbose ? println("processed results $oldoffset-$offset") : false
		end
	end
	df_pairs
end

function getStatic(dbParams::Dict)
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/PathwayCommons/rq")
    str = open(f->read(f, String), string(rqDir,"/",dbParams[:rqFile]));

	cnames = [:pw1,:pw2]
	ctypes = fill(Union{Missing,String},2)
	members = DataFrame(ctypes,cnames,0)

    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
              "Accept" => "application/sparql-results+xml"]
	resp = requestTTL(dbParams[:port],dbParams[:protocol],dbParams[:format],
						dbParams[:host],dbParams[:path],dbParams[:method],
						header,str)

	temp = parseSparqlResponse(resp)

	if size(temp,1) > 0
		members = vcat(members,temp;cols=:union)
	else
		members = temp
	end
	members
end

function getParticipant(dbParams::Dict,entity::String)
    val = delimitValues([entity],"","<>")
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/PathwayCommons/rq")
    str = open(f->read(f, String), string(rqDir,"/",dbParams[:rqFile]));

	cnames = [:superEntity,:participant,:participantType,:participantLocRef,:displayName,
			  :dbname,:id,:participantRef,
			  :participantEntRef,:standardName,:participantEntRefType,:entId,:entIdDb]
	ctypes = fill(Union{Missing,String},13)
	members = DataFrame(cnames .=> [type[] for type in ctypes])
    turtle = Mustache.render(str,Dict{Any,Any}("ent"=>val,
												"fromgraph"=>dbParams[:fromgraph]))
    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
              "Accept" => "application/sparql-results+xml"]
	resp = requestTTL(dbParams[:port],dbParams[:protocol],dbParams[:format],
						dbParams[:host],dbParams[:path],dbParams[:method],
						header,turtle)

	temp = parseSparqlResponse(resp)

	if size(temp,1) > 0
		members = vcat(members,temp;cols=:union)
	else
		members = temp
	end
	members
end

function getComplex(dbParams::Dict,cxref::String)
    val = delimitValues([cxref],"","<>")
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/PathwayCommons/rq")
	str = open(f->read(f, String), string(rqDir,"/","getComplex_gdb.rq"));

	cnames = [:cx,:comp,:compType,:participantLocRef,:cxLocRef]
	ctypes = fill(Union{Missing,String},5)
	cxEnts = DataFrame(cnames .=> [type[] for type in ctypes])

	turtle = Mustache.render(str,Dict{Any,Any}("cx"=>val,
											   "fromgraph"=>dbParams[:fromgraph]))

    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
    		  "Accept" => "application/sparql-results+xml"]
	resp = requestTTL(dbParams[:port],dbParams[:protocol],dbParams[:format],
						dbParams[:host],dbParams[:path],dbParams[:method],
						header,turtle)

	temp = parseSparqlResponse(resp)
	if size(temp,1) > 0
		cxEnts = vcat(cxEnts,temp;cols=:union)
	else
		cxEnts = temp
	end
	cxEnts
end
