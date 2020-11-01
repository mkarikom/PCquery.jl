using PCquery
using Test

@testset "escape search uri" begin
    @test HTTP.escapeuri("i+am sam.") == "i%2Bam%20sam."
    @test escapeuri("i+am sam.") == "i+am%20sam."
end
