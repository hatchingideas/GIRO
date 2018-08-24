module TestmzML

include(realpath(joinpath(@__DIR__, "..", "..", "src", "mzML", "mzML.jl")))

using GIRO.mzML, Base.Test

function testmzml()

    println(realpath(joinpath(@__DIR__, "..", "data")))

    MD = mzMLData(realpath(joinpath(@__DIR__, "..", "data")), "spectrum.xml")

    @testset begin

        @test length(MD.Spectrum) == 3

    end

end

export testmzml

end
