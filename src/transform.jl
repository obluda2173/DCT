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

        f[k+1] = sigma
    end

    return f
end


function dct_2d(M)
    h, w = size(M)

    row_result = zeros(Float64, h, w)
    final_result = zeros(Float64, h, w)

    for i in 1:h
        row_result[i, :] = dct_1d(M[i, :])
    end

    for j in 1:w
        final_result[:, j] = dct_1d(row_result[:, j])
    end

    return final_result
end


function idct_2d(C)
    h, w = size(C)

    col_result = zeros(Float64, h, w)
    final_pixels = zeros(Float64, h, w)

    for j in 1:w
        col_result[:, j] = idct_d1(C[:, j])
    end

    for i in 1:h
        final_pixels[i, :] = idct_d1(col_result[i, :])
    end

    return final_pixels
end
