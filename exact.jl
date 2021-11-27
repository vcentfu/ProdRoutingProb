using CPUTime
using JuMP
using CPLEX

const OPTIMAL = JuMP.MathOptInterface.OPTIMAL
const INFEASIBLE = JuMP.MathOptInterface.INFEASIBLE
const UNBOUNDED = JuMP.MathOptInterface.DUAL_INFEASIBLE

include("utils.jl")


function Pdi_exact_ne(data, type, m)
    mod = Model(CPLEX.Optimizer)
    u = data["u"]
    f = data["f"]
    h = data["h"]
    n = data["n"]
    l = data["l"]
    c = zeros(n + 1, n + 1)

    for i in 1:n + 1
        for j in 1:n + 1
            if type == "A"
                c[i, j] = Distance_A(data, i, j)
            elseif type == "B"
                c[i, j] = Distance_B(data, i, j)
            end
        end
    end

    d = data["d"]
    C = data["C"]
    Q = data["Q"]
    L = data["L"]
    L0 = data["L0"]

    @variable(mod, p[1:l] >= 0)
    @variable(mod, I[1:n + 1, 1:l] >= 0)
    @variable(mod, y[1:l], Bin)
    @variable(mod, z0[1:l] >= 0, Int)
    @variable(mod, z[1:n, 1:l], Bin)
    @variable(mod, x[1:n + 1, 1:n + 1, 1:l], Bin)

    for i in 1:n + 1
        delete(mod, x[i, i, :])
    end

    @variable(mod, q[1:n, 1:l] >= 0)
    @variable(mod, w[1:n, 1:l] >= 0)

    @objective(mod, Min, sum([u * p[t] + f * y[t] + sum([h[i] * I[i, t] for i in 1:n + 1]) + sum([c[i, j] * x[i, j, t] for i in 1:n + 1 for j in 1:n + 1 if i != j]) for t in 1:l]))

    @constraint(mod, c2z, L0[1] + p[1] == sum(q[1:n, 1]) + I[1, 1])
    @constraint(mod, c2[t in 1:l - 1], I[1, t] + p[t + 1] == sum(q[1:n, t + 1]) + I[1, t + 1])

    @constraint(mod, c3z[i in 1:n], L0[i + 1] + q[i, 1] == d[i, 1] + I[i + 1, 1])
    @constraint(mod, c3[i in 1:n, t in 1:l - 1], I[i + 1, t] + q[i, t + 1] == d[i, t + 1] + I[i + 1, t + 1])

    M = [min(C, sum(d[:, t:l])) for t in 1:l]
    @constraint(mod, c4[t in 1:l], p[t] <=  M[t] * y[t])

    @constraint(mod, c5[t in 1:l], I[1, t] <= L[1])

    @constraint(mod, c6z[i in 1:n], L0[i + 1] + q[i, 1] <= L[i + 1])
    @constraint(mod, c6[i in 1:n, t in 1:l - 1], I[i + 1, t] + q[i, t + 1] <= L[i + 1])

    Mt = zeros(n, l)

    for i in 1:n
        for t in 1:l
            Mt[i, t] = min(L[i + 1], Q, sum(d[i, t:l]))
        end
    end

    @constraint(mod, c7[i in 1:n, t in 1:l], q[i, t] <= Mt[i, t] * z[i, t])

    @constraint(mod, c8[i in 1:n, t in 1:l], sum([x[i + 1, j, t] for j in 1:n + 1 if j != i + 1]) == z[i, t])

    @constraint(mod, c9z[t in 1:l], sum([x[j, 1, t] for j in 1:n + 1 if j != 1]) + sum([x[1, j, t] for j in 1:n + 1 if j != 1]) == 2 * z0[t])
    @constraint(mod, c9[i in 1:n, t in 1:l], sum([x[j, i + 1, t] for j in 1:n + 1 if j != i + 1]) + sum([x[i + 1, j, t] for j in 1:n + 1 if j != i + 1]) == 2 * z[i, t])

    @constraint(mod, c10[t in 1:l], z0[t] <= m)

    @constraint(mod, c11[i in 1:n, j in 1:n, t in 1:l], w[i, t] - w[j, t] >= q[i, t] - Mt[i, t] * (1 - x[i + 1, j + 1, t]))

    for i in 1:n
        delete(mod, c11[i, i, :])
    end

    @constraint(mod, c12[i in 1:n, t in 1:l], w[i, t] <= Q * z[i, t])

    println("Affichage du modèle avant résolution")
    print(mod)
    println()

    println("Résolution")
    optimize!(mod)
    println()

    println("Récupération et affichage \"à la main\" d'informations précises")
    status = termination_status(mod)

    if status == INFEASIBLE
        println("Le problème n'est pas réalisable")
    else
        println("Affichage de tous les détails de la solution avec la commande solution_summary")
        println(solution_summary(mod, verbose = true))
        println()

        res = []

        for t in 1:l
            resp = []

            for i in 1:r + 1
                for j in 1:r + 1
                    if i != j && abs(value(x[i, j, t]) - 1) < EPS
                        push!(resp, (inds[i], inds[j]))
                    end
                end
            end

            push!(res, resp)
        end
    end

    return res, JuMP.objective_value(mod)
end


function Exact_ne(path, m)
    data = Read_data(path)

    if occursin("A_", path)
        type = "A"
    else
        type = "B"
    end

    res, robj = Pdi_exact_ne(data, type, m)

    for i in 1:length(res)
        Draw_vrp(path, data, res[i], i, "exact")
    end

    return robj
end


function Read_inst_pdi_exact_ne(s, m)
    ldir = readdir("./PRP_instances/")

    for i in 1:length(ldir)
        if occursin(s, ldir[i])
            x = split(ldir[i], ".")
            t0 = CPUtime_us()
            tmin = Exact_ne(ldir[i], m)
            t1 = CPUtime_us()

            open(string("./Results/PDI_exact_", x[1], "/", "PDI_exact_", x[1], ".log"), "w") do f
                write(f, string("Execution time : ", (t1 - t0) / 1e6, " s\n"))
                write(f, string("Objective value : ", tmin, "\n"))
            end
        end
    end
end


function Pdi_exact_e(data, type, m)
    mod = Model(CPLEX.Optimizer)
    u = data["u"]
    f = data["f"]
    h = data["h"]
    n = data["n"]
    l = data["l"]
    c = zeros(n + 1, n + 1)

    for i in 1:n + 1
        for j in 1:n + 1
            if type == "A"
                c[i, j] = Distance_A(data, i, j)
            elseif type == "B"
                c[i, j] = Distance_B(data, i, j)
            end
        end
    end

    d = data["d"]
    C = data["C"]
    Q = data["Q"]
    L = data["L"]
    L0 = data["L0"]

    @variable(mod, p[1:l] >= 0)
    @variable(mod, I[1:n + 1, 1:l] >= 0)
    @variable(mod, y[1:l], Bin)
    @variable(mod, z[1:n + 1, 1:m, 1:l], Bin)
    @variable(mod, x[1:n + 1, 1:n + 1, 1:m, 1:l], Bin)

    for i in 1:n + 1
        for k in 1:m
            delete(mod, x[i, i, k, :])
        end
    end

    @variable(mod, q[1:n, 1:m, 1:l] >= 0)

    @objective(mod, Min, sum([u * p[t] + f * y[t] + sum([h[i] * I[i, t] for i in 1:n + 1]) + sum([c[i, j] * sum(x[i, j, k, t] for k in 1:m) for i in 1:n + 1 for j in 1:n + 1 if i != j]) for t in 1:l]))

    @constraint(mod, c21z, L0[1] + p[1] == sum(q[1:n, 1:m, 1]) + I[1, 1])
    @constraint(mod, c21[t in 1:l - 1], I[1, t] + p[t + 1] == sum(q[1:n, 1:m, t + 1]) + I[1, t + 1])

    @constraint(mod, c22z[i in 1:n], L0[i + 1] + sum(q[i, 1:m, 1]) == d[i, 1] + I[i + 1, 1])
    @constraint(mod, c22[i in 1:n, t in 1:l - 1], I[i + 1, t] + sum(q[i, 1:m, t + 1]) == d[i, t + 1] + I[i + 1, t + 1])

    M = [min(C, sum(d[:, t:l])) for t in 1:l]
    @constraint(mod, c23[t in 1:l], p[t] <=  M[t] * y[t])

    @constraint(mod, c24[t in 1:l], I[1, t] <= L[1])

    @constraint(mod, c25z[i in 1:n], L0[i + 1] + sum(q[i, 1:m, 1]) <= L[i + 1])
    @constraint(mod, c25[i in 1:n, t in 1:l - 1], I[i + 1, t] + sum(q[i, 1:m, t + 1]) <= L[i + 1])

    Mt = zeros(n, l)

    for i in 1:n
        for t in 1:l
            Mt[i, t] = min(L[i + 1], Q, sum(d[i, t:l]))
        end
    end

    @constraint(mod, c26[k in 1:m, i in 1:n, t in 1:l], q[i, k, t] <= Mt[i, t] * z[i, k, t])

    @constraint(mod, c27[i in 1:n, t in 1:l], sum(z[i + 1, 1:m, t]) <= 1)

    @constraint(mod, c28[k in 1:m, i in 1:n + 1, t in 1:l], sum([x[j, i, k, t] for j in 1:n + 1 if j != i]) + sum([x[i, j, k, t] for j in 1:n + 1 if j != i]) == 2 * z[i, k, t])

    Nc = [i for i in 1:n + 1]
    PNc = Partition_s2(Nc)

    @constraint(mod, c29[i in 1:length(PNc), k in 1:m, t in 1:l], sum([x[PNc[i][ii], PNc[i][ij], k, t] for ii in 1:length(PNc[i]) for ij in 1:length(PNc[i]) if ii != ij]) <= length(PNc[i]) - 1)

    @constraint(mod, c30[k in 1:m, t in 1:l], sum(q[:, k, t]) <= Q * z[1, k, t])

    #println("Affichage du modèle avant résolution")
    #print(mod)
    #println()

    #println("Résolution")
    optimize!(mod)
    #println()

    #println("Récupération et affichage \"à la main\" d'informations précises")
    status = termination_status(mod)

    if status == INFEASIBLE
        #println("Le problème n'est pas réalisable")
    else
        #println("Affichage de tous les détails de la solution avec la commande solution_summary")
        println(solution_summary(mod, verbose = true))
        #println()

        res = []

        for t in 1:l
            resp = []

            for k in 1:m
                for i in 1:r + 1
                    for j in 1:r + 1
                        if i != j && abs(value(x[i, j, k, t]) - 1) < EPS
                            push!(resp, (inds[i], inds[j]))
                        end
                    end
                end
            end

            push!(res, resp)
        end
    end

    return res, JuMP.objective_value(mod)
end


function Exact_e(path, m)
    data = Read_data(path)

    if occursin("A_", path)
        type = "A"
    else
        type = "B"
    end

    res, robj = Pdi_exact_e(data, type, m)

    for i in 1:length(res)
        Draw_vrp(path, data, res[i], i, "exact2")
    end

    return robj
end


function Read_inst_pdi_exact_e(s, m)
    ldir = readdir("./PRP_instances/")

    for i in 1:length(ldir)
        if occursin(s, ldir[i])
            x = split(ldir[i], ".")
            t0 = CPUtime_us()
            tmin = Exact_e(ldir[i], m)
            t1 = CPUtime_us()

            open(string("./Results/PDI_exact_", x[1], "/", "PDI_exact_", x[1], ".log"), "w") do f
                write(f, string("Execution time : ", (t1 - t0) / 1e6, " s\n"))
                write(f, string("Objective value : ", tmin, "\n"))
            end
        end
    end
end
