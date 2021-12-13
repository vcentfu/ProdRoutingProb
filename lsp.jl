using JuMP
using CPLEX


const OPTIMAL = JuMP.MathOptInterface.OPTIMAL
const INFEASIBLE = JuMP.MathOptInterface.INFEASIBLE
const UNBOUNDED = JuMP.MathOptInterface.DUAL_INFEASIBLE
const M = Int(100000)


function Lsp_solv(data)

    """ Dict{String : Float} ->

        data : Les donnees extraits d'une instance PRP (cf. function read_data dans utils.jl).

        Retourne la solution optimale du LSP pour data.
    """

    m = Model(CPLEX.Optimizer)
    #set_silent(m)
    l = data["l"]
    n = data["n"]
    u = data["u"]
    f = data["f"]
    h = data["h"]
    d = data["d"]

    @variable(m, p[1:l] >= 0)
    @variable(m, y[1:l], Bin)
    @variable(m, I[1:n + 1, 1:l] >= 0)
    @variable(m, q[1:n, 1:l] >= 0)

    @objective(m, Min, sum([u * p[t] + f * y[t] + sum([h[i] * I[i + 1, t] for i in 1:n]) for t in 1:l]))

    @constraint(m, caz, data["L0"][1] + p[1] == sum(q[1:n, 1]) + I[1, 1])
    @constraint(m, ca[t in 1:l - 1], I[1, t] + p[t + 1] == sum(q[1:n, t + 1]) + I[1, t + 1])

    @constraint(m, cbz[i in 1:n], data["L0"][i + 1] + q[i, 1] == d[i, 1] + I[i + 1, 1])
    @constraint(m, cb[i in 1:n, t in 1:l - 1], I[i + 1, t] + q[i, t + 1] == d[i, t + 1] + I[i + 1, t + 1])

    @constraint(m, cc[t in 1:l], p[t] <= M * y[t])

    @constraint(m, cd[t in 1:l - 1], I[1, t] <= data["L"][1])

    @constraint(m, cez[i in 1:n], data["L0"][i + 1] + q[i, 1] <= data["L"][i + 1])
    @constraint(m, ce[i in 1:n, t in 1:l - 1], I[i + 1, t] + q[i, t + 1] <= data["L"][i + 1])

    println("Affichage du modèle avant résolution")
    println(m)
    println()

    println("Résolution")
    optimize!(m)
    println()

    println("Récupération et affichage \"à la main\" d'informations précises")
    status = termination_status(m)

    if status == INFEASIBLE
        println("Le problème n'est pas réalisable")
    else
        println("Affichage de tous les détails de la solution avec la commande solution_summary")
        println(solution_summary(m, verbose = true))
        println()
    end

    return
end


function Lsp_solv_heu(data, SC)

    """ Dict{String : Float} * Matrix{Float} -> Matrix{Float} * Cplex.Objective_value

        data : Les donnees extraits d'une instance PRP (cf. function read_data dans utils.jl).

        Retourne la solution optimale privee des couts heuristiques du LSP modifie pour data ainsi que la matrice des revendeurs a livrer.
    """

    m = Model(CPLEX.Optimizer)
    set_silent(m)
    l = data["l"]
    n = data["n"]
    u = data["u"]
    f = data["f"]
    h = data["h"]
    d = data["d"]

    @variable(m, p[1:l] >= 0)
    @variable(m, y[1:l], Bin)
    @variable(m, I[1:n + 1, 1:l] >= 0)
    @variable(m, q[1:n, 1:l] >= 0)
    @variable(m, z[1:n, 1:l], Bin)

    @objective(m, Min, sum([u * p[t] + f * y[t] + sum([h[i] * I[i, t] for i in 1:n + 1]) for t in 1:l]) + sum([SC[i, t] * z[i, t] for i in 1:n for t in 1:l]))

    @constraint(m, caz, data["L0"][1] + p[1] == sum(q[1:n, 1]) + I[1, 1])
    @constraint(m, ca[t in 1:l - 1], I[1, t] + p[t + 1] == sum(q[1:n, t + 1]) + I[1, t + 1])

    @constraint(m, cbz[i in 1:n], data["L0"][i + 1] + q[i, 1] == d[i, 1] + I[i + 1, 1])
    @constraint(m, cb[i in 1:n, t in 1:l - 1], I[i + 1, t] + q[i, t + 1] == d[i, t + 1] + I[i + 1, t + 1])

    @constraint(m, cc[t in 1:l], p[t] <= M * y[t])

    @constraint(m, cd[t in 1:l - 1], I[1, t] <= data["L"][1])

    @constraint(m, cez[i in 1:n], data["L0"][i + 1] + q[i, 1] <= data["L"][i + 1])
    @constraint(m, ce[i in 1:n, t in 1:l - 1], I[i + 1, t] + q[i, t + 1] <= data["L"][i + 1])

    @constraint(m, cf[i in 1:n, t in 1:l], q[i, t] <= M * z[i, t])

    #println("Affichage du modèle avant résolution")
    #println(m)
    #println()

    #println("Résolution")
    optimize!(m)
    #println()

    #println("Récupération et affichage \"à la main\" d'informations précises")
    status = termination_status(m)

    if status == INFEASIBLE
        #println("Le problème n'est pas réalisable")
    else
        #println("Affichage de tous les détails de la solution avec la commande solution_summary")
        #println(solution_summary(m, verbose = true))
        #println()
    end

    return JuMP.value.(z), JuMP.objective_value(m) - sum([SC[i, t] * JuMP.value.(z[i, t]) for i in 1:n for t in 1:l])
end
