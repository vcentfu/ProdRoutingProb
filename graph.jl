#ENV["GKSwstype"] = "nul"
using Plots


function collect_log(type, met, sizes = ["15", "50", "100"], carac = ["o", "t", "g"])
    ldir = readdir(string("./Results/", met, "/"))
    ldir = [ldirt for ldirt in ldir if occursin(type, ldirt)]
    ldir = sort(ldir)
    log = Dict()
    is = 0
    ic = 0

    for s in sizes
        log[s] = Dict()

        for c in carac
            log[s][c] = []
        end
    end

    i = 0

    for path in ldir
        for s in sizes
            if occursin(string("_", s, "_"), path)
                is = s
                break
            end
        end

        if i % 5 == 0
            for c in carac
                push!(log[is][c], [])
            end

            i = 0
        end

        open(string("./Results/", met, "/", path, "/", path, ".log"), "r") do f
            for line in readlines(f)
                tline = split(line, " ")

                for c in carac
                    if tline[1] == c
                        ic = c
                        break
                    end
                end

                if ic != "g"
                    push!(log[is][ic][end], parse(Float64, tline[5]))
                else
                    push!(log[is][ic][end], parse(Float64, tline[5][1:end - 1]))
                end
            end
        end

        i = i + 1
    end

    return log
end


function Draw_perf_graph_A(ins_size, file_name, root_name = "Graphs")
    if ins_size <= 50
        log_heu = collect_log("A", "heu")
        log_exact = collect_log("A", "exact")
        lv = [1.5 + i - 1 for i in 1:length(log_exact[string(ins_size)]["o"])]


        x = Any[round(0.6 + 0.2 * j + i - 1.2, digits = 1) for i in 1:length(log_exact[string(ins_size)]["o"]) for j in 1:length(log_exact[string(ins_size)]["o"][1])]
        o = []
        g = []
        a = []
        t1 = []
        t2 = []

        for z in zip(log_exact[string(ins_size)]["o"], log_exact[string(ins_size)]["g"], log_heu[string(ins_size)]["o"])
            bo, bg, ba = z

            for y in zip(bo, bg, ba)
                io, ig, ia = y
                push!(o, io)
                push!(g, io * (1 - ig / 100))
                push!(a, ia)
            end
        end

        t = zeros(length(x), 3)

        for i in 1:length(x)
            t[i, 1] = o[i]
            t[i, 2] = g[i]
            t[i, 3] = a[i]
        end

        if !ispath(string("./", root_name))
            mkdir(string("./", root_name))
        end

        scatter(x, t, shape=[:+ :hline :x], xlabel = "Classe", ylabel = "Valeur Objectif", label = ["branch and cut avec config. CPLEX" "Borne relaxation" "Méthode approchée"])
        plot!(lv, seriestype="vline", label="")
        plot!(size=(800,800))
        savefig(string("./", root_name, "/perf_", file_name, ".png"))

        for z in zip(log_exact[string(ins_size)]["t"],log_heu[string(ins_size)]["t"])
            s1, s2 = z

            for y in zip(s1, s2)
                y1, y2 = y
                push!(t1, y1)
                push!(t2, y2)
            end
        end

        t = zeros(length(x), 2)

        for i in 1:length(x)
            t[i, 1] = t1[i]
            t[i, 2] = t2[i]
        end

        scatter(x, t, xlabel = "Classe", ylabel = "Temps d'exécution", label = ["branch and cut avec config. CPLEX" "Méthode approchée"])
        plot!(lv, seriestype="vline", label="")
        plot!(size=(800,800))
        savefig(string("./", root_name, "/time_", file_name, ".png"))
    else
        log_heu = collect_log("A", "heu")
        lv = [1.5 + i - 1 for i in 1:length(log_heu[string(ins_size)]["o"])]

        x = Any[round(0.6 + 0.2 * j + i - 1.2, digits = 1) for i in 1:length(log_heu[string(ins_size)]["o"]) for j in 1:length(log_heu[string(ins_size)]["o"][1])]
        o = []
        t1 = []

        for z in zip(log_heu[string(ins_size)]["o"], log_heu[string(ins_size)]["t"])
            bo, bt = z

            for y in zip(bo, bt)
                io, it= y
                push!(o, io)
                push!(t1, it)
            end
        end

        if !ispath(string("./", root_name))
            mkdir(string("./", root_name))
        end

        scatter(x, o, xlabel = "Classe", ylabel = "Valeur Objectif", label = "Méthode approchée")
        plot!(lv, seriestype="vline", label="")
        plot!(size=(800,800))
        savefig(string("./", root_name, "/perf_", file_name, ".png"))

        scatter(x, t1, xlabel = "Classe", ylabel = "Temps d'exécution", label = "Méthode approchée")
        plot!(lv, seriestype="vline", label="")
        plot!(size=(800,800))
        savefig(string("./", root_name, "/time_", file_name, ".png"))
    end

    return
end


function draw_all_graphs_A()
    for s in [15, 50, 100]
        Draw_perf_graph_A(s, string("graph_", s))
    end

    return
end
