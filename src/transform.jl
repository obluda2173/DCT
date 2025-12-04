
alpha(u, N) = return (u == 0) ? sqrt(1/N) : sqrt(2/N)

function dct_1d(f)
    N = length(f)
    C = zeros(Float64, N)

    for u in 0:(N-1)
        sigma = 0.0

        for x in 0:(N-1)
            cos_term = cos((pi * (2*x + 1) * u)/(2*N))
            sigma += f[x+1] * cos_term
        end

        C[u+1] = alpha(u, N) * sigma
    end

    return C
end


function idct_d1(C)
    N = length(C)
    f = zeros(Float64, N)

    for k in 0:(N-1)
        sigma = 0.0

        for u in 0:(N-1)
            cos_term = cos((pi * (2*k + 1) * u)/(2*N))
            sigma += alpha(u, N) * C[u+1] * cos_term
        end

        f[u+1] = sigma
    end

    return f
end

