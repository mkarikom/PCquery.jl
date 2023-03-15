using PCquery
using HTTP
using Test

@testset "escape search uri" begin
    @test HTTP.escapeuri("i+am sam.") == "i%2Bam%20sam."
    @test PCquery.escapeuri("i+am sam.") == "i+am%20sam."
end
