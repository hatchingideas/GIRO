function write_rtadj_trafoxml(FilePath :: String, DictRTAdj :: Dict)

    RTLoc = [k for k in keys(DictRTAdj)]

    NumRTLoc = length(RTLoc)

    open(FilePath, "w") do f

        println(f, """<?xml version="1.0" encoding="UTF-8"?>""")
        println(f, """<TrafoXML version="1.0" xsi:noNamespaceSchemaLocation="http://open-ms.sourceforge.net/schemas/TrafoXML_1_0.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">""")
	    println(f, """\t<Transformation name="b_spline">""")
        println(f, """\t\t<Param  type="float" name="wavelength" value="0"/>""")
		println(f, """\t\t<Param  type="int" name="num_nodes" value="5"/>""")
		println(f, """\t\t<Param  type="string" name="extrapolate" value="linear"/>""")
		println(f, """\t\t<Param  type="int" name="boundary_condition" value="2"/>""")
        println(f, """\t\t<Pairs count="$NumRTLoc">""")
        [println(f, """\t\t\t<Pair from="$(RTLoc[i])" to="$(DictRTAdj[RTLoc[i]])"/>""") for i in 1:NumRTLoc]
        println(f, """\t\t</Pairs>""")
        println(f, """\t</Transformation>""")
        println(f, """</TrafoXML>""")

    end

    nothing

end
