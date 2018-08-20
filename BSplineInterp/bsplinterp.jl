function bsplinterp(x :: Vector{Float64}, Knots :: Vector{Float64}, SplOrder = 4 :: Int64)

    NumKnots = length(Knots)

    BSplInterp = zeros(Float64, length(x), SplOrder)

    for i in 0:SplOrder

        NumCP = NumKnots - i - 1

        for j = 1:NumCP

            if i == 0

                BSplInterp[ (x .>= Knots[j]) .& (x .< Knots[j+1]) , j] = 1

            else

                BSplInterp[(x .>= Knots[j]) .& (x .< Knots[i+j]) , j] = BSplInterp[(x .>= Knots[j]) .& (x .< Knots[i+j]), j] .* (x[(x .>= Knots[j]) .& (x .< Knots[i+j])] - Knots[j]) / (Knots[i+j] - Knots[j]);
                BSplInterp[(x .>= Knots[j+1]) .& (x .< Knots[i+j+1]), j] = BSplInterp[(x .>= Knots[j+1]) .& (x .< Knots[i+j+1]), j] + BSplInterp[(x .>= Knots[j+1]) .& (x .< Knots[i+j+1]), j+1] .* (Knots[i+j+1] - x[(x .>= Knots[j+1]) .& (x .< Knots[i+j+1])]) / (Knots[i+j+1] - Knots[j+1]);

            end

        end

    end

    BSplInterp[:,1]

end
