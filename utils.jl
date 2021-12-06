using LinearAlgebra
using LightGraphs, SimpleWeightedGraphs, Cairo, Compose, Fontconfig, Colors
using GraphPlot

const MCOLORS = [colorant"red", colorant"blue", colorant"green", colorant"yellow", colorant"orange", colorant"purple", colorant"tan", colorant"pink", colorant"cyan", colorant"silver"]


function Read_data(path)
    data = Dict()
    nend = false

    open(string("./PRP_instances/", path)) do f
        for line in readlines(f)[2:end]
            tline = split(line, " ")

            if length(tline) == 2 && !nend
                data[tline[1]] = Int(parse(Float64, tline[2]))
            elseif "h" in tline
                if tline[1] == "0"
                    data["x"] = []
                    data["y"] = []
                    data["h"] = []
                    data["L"] = []
                    data["L0"] = []
                end

                push!(data["x"], Int(parse(Float64, tline[2])))
                push!(data["y"], Int(parse(Float64, tline[3])))
                push!(data["h"], Int(parse(Float64, tline[6])))
                push!(data["L"], Int(parse(Float64, tline[8])))
                push!(data["L0"], Int(parse(Float64, tline[10])))
            elseif tline[1] == "d"
                nend = true
                data["d"] = []
                continue
            elseif nend
                d = map(x -> parse(Int, x), tline[2:end - 1])
                push!(data["d"], d)
            end
        end
    end

    d = data["d"]
    n = data["n"]
    l = data["l"]
    md = zeros(Int64, n, l)

    for i in 1:n
        for t in 1:l
            md[i, t] = d[i][t]
        end
    end

    data["d"] = md

    return data
end


function Distance_A(data, i, j)
    x = data["x"]
    y = data["y"]

    return floor(sqrt(((x[i] - x[j]) ^ 2) + ((y[i] - y[j]) ^ 2)) + 1 / 2)
end


function Distance_B(data, i, j)
    x = data["x"]
    y = data["y"]
    mc = data["mc"]

    return mc * sqrt(((x[i] - x[j]) ^ 2) + ((y[i] - y[j]) ^ 2))
end


function Initialize_SC(data, type)
    n = data["n"]
    l = data["l"]
    res = zeros(Float64, n, l)

    for i in 1:n
        for t in 1:l
            if type == "A"
                res[i, t] = Distance_A(data, 1, i + 1) + Distance_A(data, i + 1, 1)
            elseif type == "B"
                res[i, t] = Distance_B(data, 1, i + 1) + Distance_B(data, i + 1, 1)
            end
        end
    end

    return res
end


function Update_SC(data, type, vrp)
    n = data["n"]
    l = data["l"]
    res = zeros(Float64, n, l)

    for i in 1:n
        for t in 1:l
            e = -1
            s = -1

            for (l, k) in vrp[t]
                if l == i
                    s = k
                elseif k == i
                    e = l
                end
            end

            if e == -1 && s == -1
                continue
            end

            if type == "A"
                res[i, t] = Distance_A(data, e + 1, i + 1) + Distance_A(data, i + 1, s + 1) - Distance_A(data, e + 1, s + 1)
            elseif type == "B"
                res[i, t] = Distance_B(data, e + 1, i + 1) + Distance_B(data, i + 1, s + 1) - Distance_B(data, e + 1, s + 1)
            end
        end
    end

    return res
end


function Draw_vrp(path, data, vrp, tp, opts)
    n = data["n"] + 1
    g = DiGraph(n)
    xvrp = copy(vrp)
    t = []
    f = true
    r = -1
    Ncg = []

    for (i, j) in vrp
        if i + 1 != 1 && !(i + 1 in Ncg)
            push!(Ncg, i + 1)
        elseif j + 1 != 1 && !(j + 1 in Ncg)
            push!(Ncg, j + 1)
        end
    end

    while f
        f = false
        mt = [0]

        for (i, j) in xvrp
            if i == 0
                push!(mt, j)
                deleteat!(xvrp, findall(x -> x == (i, j), xvrp))
                f = true
                break
            end
        end

        if f
            while mt[end] != 0
                for (i, j) in xvrp
                    if i == mt[end]
                        push!(mt, j)
                        break
                    end
                end
            end

        push!(t, mt)
        end
    end

    #println(t)

    for tr in t
        for i in 1:length(tr) - 1
            add_edge!(g, tr[i] + 1, tr[i + 1] + 1)
        end
    end

    ecolors = []

    for e in edges(g)
        ent = src(e)
        sor = dst(e)

        for i in 1:length(t)
            for j in 1:length(t[i]) - 1
                if t[i][j] + 1 == ent && t[i][j + 1] + 1 == sor
                    push!(ecolors, MCOLORS[i])
                end
            end
        end
    end

    tx = [x for x in data["x"]]
    ty = [y for y in data["y"]]
    nodec = [colorant"turquoise" for i in 1:n + 1]
    nodec[1] = colorant"yellow"

    for k in Ncg
        nodec[k] = colorant"orange"
    end

    if !ispath("./Results")
        mkdir("./Results")
    end

    testname = split(path, ".")
    file = string("PDI_", opts)

    if !ispath(string("./Results/", file, "_", testname[1]))
        mkdir(string("./Results/", file,  "_", testname[1]))
    end

    draw(PDF(string("./Results/", file, "_", testname[1], "/", "VRP_", opts, "_", testname[1], "_", "p", tp, ".pdf"), 16cm, 16cm), gplot(g, tx, ty, nodelabel = 1:nv(g), edgestrokec = ecolors, nodefillc = nodec))

    return
end


function Partition_s2(l)
    res = []

    for k in l
        ti = length(res)

        for i in 1:ti
            ts = deepcopy(res[i])
            push!(ts, k)
            push!(res, ts)
        end

        push!(res, [k])
    end

    push!(res, [])

    fres = []

    for tres in res
        if length(tres) >= 2
            push!(fres, tres)
        end
    end

    return fres
end
