function get_rt_adj_vec(RT_AR :: RTAdjRec)

    (RT_AR.AdjMat*RT_AR.BsplCP)[RT_AR.StartIdx : RT_AR.EndIdx, :]

end
