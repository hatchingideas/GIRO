module TestmzML

using Revise

using GIRO.mzML, Test

function testmzml()

    println(realpath(joinpath(@__DIR__, "..", "data")))

    MD = mzMLData(realpath(joinpath(@__DIR__, "..", "data")), "spectrum.xml")

    @testset begin

        @test length(MD.Spectrum) == 3

    end

end

testmzml()

export testmzml

end
