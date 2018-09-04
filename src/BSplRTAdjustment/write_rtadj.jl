function write_rtadj(FilePath :: String,  RTA :: RTAdjRec)

    # interpolate to recover the deformation at scan_start_time:
    RTAdj = get_rt_adj_vec(ResLevelDeformFieldParamVec[i])[StartIdx, EndIdx] * RTRes

    RTAdjInterp = RT#RTVec[i]

    D = Dict("Retention_Time" => RT,
             "Adjustment_in_Second" => RTAdj)





             ResLevelDeformFieldParamVec

             using Plots

             p = plot(1:2048, get_rt_adj_vec(ResLevelDeformFieldParamVec[3]))

             display(p)

             println("Writing out retention time adjustments in csv files:")
             for i in 1:NumImg


                 writecsv(joinpath(FileDir, "$(FileName[i][1:end-4])csv"), D)

             end

             println("GIRO successfully finished. Output csv files written into directory $FileDir.")








    writecsv(FilePath, D)

    nothing

end
