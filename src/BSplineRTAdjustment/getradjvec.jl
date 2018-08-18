function getadjvec(RT_AR :: RTAdjRec)

    (RT_AR.AdjMat*RT_AR.BsplCP)[RT_AR.StartIdx : RT_AR.EndIdx, :]

end
