@inline issafe(c::Char) = c == '-' ||
                          c == '.' ||
                          c == '_' ||
                          c == '+' || # added '+' to the list of characters to ignore
                          c == '&' ||
                          c == '=' ||
                          c == '*' ||
                          (isascii(c) && (isletter(c) || isnumeric(c)))

utf8_chars(str::AbstractString) = (Char(c) for c in IOExtras.bytes(str))

"percent-encode a string, dict, or pair for a uri"
function escapeuri end

escapeuri(c::Char) = string('%', uppercase(string(Int(c), base=16, pad=2)))
escapeuri(str::AbstractString, safe::Function=issafe) =
  join(safe(c) ? c : escapeuri(c) for c in utf8_chars(str))

escapeuri(bytes::Vector{UInt8}) = bytes
escapeuri(v::Number) = escapeuri(string(v))
escapeuri(v::Symbol) = escapeuri(string(v))

escapeuri(key, value) = string(escapeuri(key), "=", escapeuri(value))
escapeuri(key, values::Vector) = escapeuri(key => v for v in values)
escapeuri(query) = join((escapeuri(k, v) for (k,v) in query), "&")
escapeuri(nt::NamedTuple) = escapeuri(pairs(nt))
