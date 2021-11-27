using JuMP
using CPLEX
include("utils.jl")


const EPS = 1e-10
const OPTIMAL = JuMP.MathOptInterface.OPTIMAL
const INFEASIBLE = JuMP.MathOptInterface.INFEASIBLE
const UNBOUNDED = JuMP.MathOptInterface.DUAL_INFEASIBLE;


function Vrp_local(data, type, Nc, m, t)
    mod = Model(CPLEX.Optimizer)
    Q = data["Q"]
    d = data["d"][:, t]
    inds = append!([0], [i for i in 1:length(Nc[:, t]) if Nc[i, t] == 1])
    r = length(inds) - 1

    if r == 0
        #println("Pas de tournée")
        return [], 0
    end

    #println(inds)

    @variable(mod, 0 <= w[1:r] <= Q)
    @variable(mod, x[1:r + 1, 1: r + 1], Bin)

    for i in 1:r + 1
        delete(mod, x[i, i])
    end

    if type == "A"
        @objective(mod, Min, sum([Distance_A(data, inds[i] + 1, inds[j] + 1) * x[i, j] for i in 1:r + 1 for j in 1:r + 1]))
    elseif type == "B"
        @objective(mod, Min, sum([Distance_B(data, inds[i] + 1, inds[j] + 1) * x[i, j] for i in 1:r + 1 for j in 1:r + 1]))
    end

    @constraint(mod, ca, sum(x[1, 2:r + 1]) <= m)
    @constraint(mod, cb, sum(x[2:r + 1, 1]) <= m)
    @constraint(mod, cc[i in 1:r], sum([x[i + 1, j] for j in 1:r + 1 if j != i + 1]) == 1)
    @constraint(mod, cd[j in 1:r], sum([x[i, j + 1] for i in 1:r + 1 if i != j + 1]) == 1)
    @constraint(mod, ce[i in 1:r, j in 1:r], w[i] - w[j] >= d[inds[2:r + 1][i]] - (Q + d[inds[2:r + 1][i]]) * (1 - x[i + 1, j + 1]))

    for i in 1:r
        delete(mod, ce[i, i])
    end

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
        #println(solution_summary(mod, verbose = true))
        #println()

        res = []

        for i in 1:r + 1
            for j in 1:r + 1
                if i != j && abs(value(x[i, j]) - 1) < EPS
                    push!(res, (inds[i], inds[j]))
                end
            end
        end

        return res, JuMP.objective_value(mod)
    end

    return -1, -1
end
