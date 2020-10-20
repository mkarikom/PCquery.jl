# get a vector of simple physical entity dictionaries
function getNested(dbParams::Dict,entity::String,
					entNames::Vector,entNamesSimple::Vector)
	reflist = []
	dbParams[:rqFile] = "getParticipantNested_gdb.rq"
	resp = getParticipant(dbParams::Dict,entity::String)
	if size(resp,1) == 0 # simple ref was provided
		dbParams[:rqFile] = "getParticipantSimple_gdb.rq"
		resp = getParticipant(dbParams::Dict,entity::String)
	else # simple refs were retrieved from nested
		indMissing = findall(n->ismissing(n),collect(resp[1,:]))
		@assert length(indMissing) == 0 "all nested participants must be resolved"
	end
	# loop over refs and recurse any complexes
	for i in 1:size(resp,1)
		simpleNames = [dbParams[:physicalEntity]...,dbParams[:simpleEntity]...]
		inds = indexin(simpleNames,Symbol.(names(resp)))
		ind_found = findall(r->!isnothing(r),inds)
		ind_notfound = findall(r->isnothing(r),inds)
		# determine if resp[i,:] is a complex or a non-complex
		if length(ind_found) < length(simpleNames)
			println("adding complex row $i")
			cxref =  resp[i,dbParams[:physicalEntity][4]]
			members = decomposeComplex(dbParams,cxref)
            cxdict = Dict([entNames...,:members] .=> [collect(resp[i,dbParams[:physicalEntity]])...,members])
			push!(reflist,cxdict)
		else
			println("adding participant row $i")
			push!(reflist,Dict([entNames...,entNamesSimple...] .=>
				collect(resp[i,[dbParams[:physicalEntity]...,dbParams[:simpleEntity]...]])))
		end
	end
	# @assert size(reflist,1) > 0 "require at least one simple ref"
    reflist
end

# decompose all complexes in the interaction list
function decomposeComplex(dbParams::Dict,cxref::String)
	cxRefs = getComplex(dbParams,cxref)
	members = []
	for m in 1:size(cxRefs,1)
        # need to include getNested here, temporarily omit
        nonNested = getNested(dbParams,cxRefs[m,:comp],
                                dbParams[:physicalEntity],
                                dbParams[:simpleEntity])
        for nn in 1:length(nonNested)
            push!(members,nonNested[nn])
        end    
	end
	members
end

function getComplex(dbParams::Dict,cxref::String)
    val = delimitValues([cxref],"","<>")
    srcDir = join(split(pathof(PCquery),"/")[1:end-1],"/")
    rqDir = string(srcDir,"/rq")
	str = open(f->read(f, String), "$rqDir/getComplex_gdb.rq");
    turtle = Mustache.render(str,Dict{Any,Any}("cx"=>val))

    # compose the header and execute query
    header = ["Content-Type" => "application/x-www-form-urlencoded",
    		  "Accept" => "application/sparql-results+xml"]
    fmt = "application/sparql-results+xml"
    resp = PCquery.requestTTL(dbParams[:port],dbParams[:protocol],fmt,
                        dbParams[:host],dbParams[:path],dbParams[:method],
                        header,turtle)
    cxEnts = PCquery.parseSparqlResponse(resp)
end
