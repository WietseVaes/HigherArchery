include("bin_tree.jl");
include("ID.jl");
function HBS_ID(A, k)
    Ares, Ut, y, z, Vt, B = HBS(A, k)
    N = size(A, 2)
    levels = bin_tree(N, k)
    L = length(levels)

    Usamp = Dict()
    Vsamp = Dict()
    Bskel = Dict()

    for l in 2:(L - 1)
        local_level = levels[l]
        for tau1 in 1:length(local_level)
            Ytmp = Ut[l][tau1] * y[l][tau1]
            Ztmp = Vt[l][tau1] * z[l][tau1]

            _, Jin = ID(Ytmp')
            _, Jout = ID(Ztmp')

            Usamp[l] = get(Usamp, l, Dict())
            Vsamp[l] = get(Vsamp, l, Dict())

            Usamp[l][tau1] = Ut[l][tau1][Jin[1:k], :]
            Vsamp[l][tau1] = Vt[l][tau1][Jout[1:k], :]
        end
    end

    for l in (L - 1):-1:3
        local_level = levels[l]
        for tau1 in 1:2:(length(local_level) - 1)
            alpha = tau1
            beta = tau1 + 1

            Uzeroab = zeros(size(Usamp[l][alpha], 1), size(Usamp[l][beta], 2))
            Uzeroba = zeros(size(Usamp[l][beta], 1), size(Usamp[l][alpha], 2))

            Vzeroab = zeros(size(Vsamp[l][alpha], 1), size(Vsamp[l][beta], 2))
            Vzeroba = zeros(size(Vsamp[l][beta], 1), size(Vsamp[l][alpha], 2))

            Utmp = [Usamp[l][alpha] Uzeroab; Uzeroba Usamp[l][beta]] * Ut[l - 1][div(tau1 + 1, 2)]
            Vtmp = [Vsamp[l][alpha] Vzeroab; Vzeroba Vsamp[l][beta]] * Vt[l - 1][div(tau1 + 1, 2)]

            Ytmp = Utmp * y[l - 1][div(tau1 + 1, 2)]
            Ztmp = Vtmp * z[l - 1][div(tau1 + 1, 2)]

            Tin, Jin = ID(Ytmp')
            Tout, Jout = ID(Ztmp')

            Usamp[l - 1] = get(Usamp, l - 1, Dict())
            Vsamp[l - 1] = get(Vsamp, l - 1, Dict())

            Usamp[l - 1][div(tau1 + 1, 2)] = Utmp[Jin[1:k], :]
            Vsamp[l - 1][div(tau1 + 1, 2)] = Vtmp[Jout[1:k], :]

            Bskel[l] = get(Bskel, l, Dict())

            Bskel[l][alpha] = Ut[l][alpha] * B[l][alpha] * (Vt[l][beta])'
            Bskel[l][beta] = Ut[l][beta] * B[l][beta] * (Vt[l][alpha])'
        end
    end

    alpha = 1
    beta = 2

    Bskel[2] = get(Bskel, 2, Dict())
    Bskel[2][alpha] = Ut[2][alpha] * B[2][alpha] * (Vt[2][beta])'
    Bskel[2][beta] = Ut[2][beta] * B[2][beta] * (Vt[2][alpha])'

    return Ares
end
