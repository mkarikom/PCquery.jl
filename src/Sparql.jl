function delimitedAS(paths::Vector{String},varname::String)
	xs = map(x->"""('$x')""",paths)
	xs = join(xs," ")
	xs = "{$xs}"
	xs = "($varname) $xs"
end
