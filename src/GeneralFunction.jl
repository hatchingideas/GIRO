function flatmap(f :: Function, v...)

    m = map(f, v...)

    r = vcat(m...)

end

function anscombe(x)

    sqrt.(x .+ .375)

end

function readinspecifiedlines(f :: IO, StartMark :: String, EndMark :: String)

    Buffer = IOBuffer()

    ReadFlag = false

    ContainStartMark = false

    while !eof(f)

        Str = readline(f)

        if contains(Str, StartMark)

            ReadFlag = true

            ContainStartMark = true

        end

        if (ReadFlag == true) && contains(Str, EndMark)

            write(Buffer, Str)

            return String(Buffer)

        else

            nothing

        end

        if ReadFlag == true

            write(Buffer, Str, "\n")

        end

    end

    eof(f) || !ContainStartMark ? nothing : throw(ErrorException("Marks Not Found!"))

end
