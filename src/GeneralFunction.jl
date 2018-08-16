function flatmap(f :: Function, v...)

    m = map(f, v...)

    r = vcat(m...)

end
