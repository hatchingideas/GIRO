println("Starting to run unit test.")

include(joinpath("TestmzML", "TestmzML.jl"))
using TestmzML
a = testmzml()

include(joinpath("TestImageRepresentation", "TestImageRepresentation.jl"))
using TestImageRepresentation
testimagerepresentation()

println("All unit tests passed successfully.")
