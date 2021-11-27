using CPUTime
include("utils.jl")
include("lsp.jl")
include("vrp.jl")


function Heuris(path, m)
    data = Read_data(path)

    if occursin("A_", path)
        type = "A"
    else
        type = "B"
    end

    SC = Initialize_SC(data, type)
    Nc, opt = Lsp_solv_heu(data, SC)
    vrp = []

    for p in 1:data["l"]
        rvrp, rval = Vrp_local(data, type, Nc, m, p)
        push!(vrp, rvrp)
    end

    tt = []
    tvrp = []
    indx = -1

    for i in 1:15
        SC = Update_SC(data, type, vrp)
        Nc, opt = Lsp_solv_heu(data, SC)
        vrp = []
        ri = 0

        for p in 1:data["l"]
            rvrp, rval = Vrp_local(data, type, Nc, m, p)
            ri = ri + rval
            push!(vrp, rvrp)
        end

        if -1 in vrp
            #println("solution non realisable")
            return -1
        end

        ropt = opt + ri
        push!(tt, ropt)
        push!(tvrp, vrp)

        if length(tt) >= 3
            if tt[end-2] == tt[end]
                indx = i
                break
            end
        end
    end

    #println("i = ", indx)
    #println(tt)

    indm = argmin(tt)
    pvrp = tvrp[indm]

    for i in 1:length(pvrp)
        Draw_vrp(path, data, pvrp[i], i, "heu")
    end

    return tt[indm]
end


function Read_inst_pdi_heu(s, m)
    ldir = readdir("./PRP_instances/")

    for i in 1:length(ldir)
        if occursin(s, ldir[i])
            x = split(ldir[i], ".")
            t0 = CPUtime_us()
            tmin = Heuris(ldir[i], m)
            t1 = CPUtime_us()

            open(string("./Results/PDI_heu_", x[1], "/", "PDI_heu_", x[1], ".log"), "w") do f
                write(f, string("Execution time : ", (t1 - t0) / 1e6, " s\n"))
                write(f, string("Objective value : ", tmin, "\n"))
            end
        end
    end
end
