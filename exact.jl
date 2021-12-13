using CPUTime
using JuMP
using CPLEX
using Random

const OPTIMAL = JuMP.MathOptInterface.OPTIMAL
const INFEASIBLE = JuMP.MathOptInterface.INFEASIBLE
const UNBOUNDED = JuMP.MathOptInterface.DUAL_INFEASIBLE
const EPS = 1e-10

include("utils.jl")


function Pdi_exact_ne(data, type, m)
    mod = Model(CPLEX.Optimizer)
    #set_silent(mod)
    if data["n"] <= 14
        set_optimizer_attribute(mod, "CPX_PARAM_TILIM" , 60 * 5)
    end

    set_optimizer_attribute(mod, "CPX_PARAM_EPINT" , 1e-10)

    if data["n"] > 14
        set_optimizer_attribute(mod, "CPX_PARAM_INTSOLLIM" , 5)
    end
    #println(time_limit_sec(mod))
    u = data["u"]
    f = data["f"]
    h = data["h"]
    n = data["n"]
    l = data["l"]
    c = zeros(n + 1, n + 1)
    nbaddc = 0
    nbaddc2 = 0

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

    @variable(mod, p[1:l] >= 0)  # OK
    @variable(mod, I[1:n + 1, 1:l] >= 0)  #  OK
    @variable(mod, y[1:l], Bin)           # OK
    @variable(mod, z0[1:l] >= 0, Int)     # OK
    @variable(mod, z[1:n, 1:l], Bin)      # OK
    @variable(mod, x[1:n + 1, 1:n + 1, 1:l], Bin)    # OK

    for i in 1:n + 1
        delete(mod, x[i, i, :])
    end

    @variable(mod, q[1:n, 1:l] >= 0)   # OK
    #@variable(mod, w[1:n + 1, 1:l] >= 0)   # OK

    @objective(mod, Min, sum([u * p[t] + f * y[t] + sum([h[i] * I[i, t] for i in 1:n + 1]) + sum([c[i, j] * x[i, j, t] for i in 1:n + 1 for j in 1:n + 1 if i != j]) for t in 1:l]))  # OK

    @constraint(mod, c2z, L0[1] + p[1] == sum(q[1:n, 1]) + I[1, 1])  # OK
    @constraint(mod, c2[t in 1:l - 1], I[1, t] + p[t + 1] == sum(q[1:n, t + 1]) + I[1, t + 1])  # OK

    @constraint(mod, c3z[i in 1:n], L0[i + 1] + q[i, 1] == d[i, 1] + I[i + 1, 1])  # OK
    @constraint(mod, c3[i in 1:n, t in 1:l - 1], I[i + 1, t] + q[i, t + 1] == d[i, t + 1] + I[i + 1, t + 1])  # OK

    M = [min(C, sum(d[:, t:l])) for t in 1:l]
    @constraint(mod, c4[t in 1:l], p[t] <=  M[t] * y[t])  # OK

    @constraint(mod, c5[t in 1:l], I[1, t] <= L[1])  # OK

    @constraint(mod, c6z[i in 1:n], L0[i + 1] + q[i, 1] <= L[i + 1])  # OK
    @constraint(mod, c6[i in 1:n, t in 1:l - 1], I[i + 1, t] + q[i, t + 1] <= L[i + 1])  #  OK

    Mt = zeros(n, l)

    for i in 1:n
        for t in 1:l
            Mt[i, t] = min(L[i + 1], Q, sum(d[i, t:l]))
        end
    end

    @constraint(mod, c7[i in 1:n, t in 1:l], q[i, t] <= Mt[i, t] * z[i, t])  # OK
    @constraint(mod, c7bis[t in 1:l], sum([q[i, t] for i in 1:n]) <= Q * z0[t])

    @constraint(mod, c8[i in 1:n, t in 1:l], sum([x[i + 1, j, t] for j in 1:n + 1 if j != i + 1]) == z[i, t]) # OK

    @constraint(mod, c9z[t in 1:l], sum([x[j, 1, t] for j in 1:n + 1 if j != 1]) + sum([x[1, j, t] for j in 1:n + 1 if j != 1]) == 2 * z0[t])  # OK
    @constraint(mod, c9[i in 1:n, t in 1:l], sum([x[j, i + 1, t] for j in 1:n + 1 if j != i + 1]) + sum([x[i + 1, j, t] for j in 1:n + 1 if j != i + 1]) == 2 * z[i, t]) # OK

    @constraint(mod, c10[t in 1:l], z0[t] <= m) # OK

    """
    @constraint(mod, c11z[j in 1:n, t in 1:l], w[1, t] - w[j + 1, t] >= 0)
    @constraint(mod, c11z2[i in 1:n, t in 1:l], w[i + 1, t] - w[1, t] >= q[i, t] - Mt[i, t] * (1 - x[i + 1, 1, t]))
    @constraint(mod, c11[i in 1:n, j in 1:n, t in 1:l], w[i + 1, t] - w[j + 1, t] >= q[i, t] - Mt[i, t] * (1 - x[i + 1, j + 1, t]))

    for i in 1:n
        delete(mod, c11[i, i, :])
    end

    @constraint(mod, c12[i in 1:n, t in 1:l], w[i + 1, t] <= Q * z[i, t])
    """

    function lazySep(cb_data)
        #println("test")
        xsep = zeros(n + 1, n + 1, l)

        for i in 1:n + 1
            for j in 1:n + 1
                for t in 1:l
                    if i != j
                        if callback_value(cb_data, x[i, j, t]) > 0.999
                            xsep[i, j, t] = 1
                        end
                    end
                end
            end
        end

        zsep = zeros(n, l)
        qsep = zeros(n, l)

        for i in 1:n
            for t in 1:l
                if callback_value(cb_data, z[i, t]) > 0.999
                    zsep[i, t] = 1
                end

                if callback_value(cb_data, q[i, t]) > 0.999
                    qsep[i, t] = callback_value(cb_data, q[i, t])
                end
            end
        end

        function FCCs(S, t)
            if sum([xsep[i, j + 1, t] for i in 1:n + 1 for j in S if !(i - 1 in S)]) < sum(qsep[S, t]) / Q
                return true
            end

            return false
        end

        function GFSECs(S, t)
            if sum([xsep[i + 1, j + 1, t] for i in S for j in S if i != j]) > length(S) - sum(qsep[S, t]) / Q
                return true
            end

            return false
        end

        res = []

        for t in 1:l
            resp = []

            for i in 1:n + 1
                for j in 1:n + 1
                    if i != j && abs(callback_value(cb_data, x[i, j, t]) - 1) < EPS
                        push!(resp, (i - 1, j - 1))
                    end
                end
            end

            push!(res, resp)
        end


        for t in 1:l
            R, lS = Detect_subtour(res[t])
            lS2 = [S[1:end - 1] for S in lS]


            for S in lS
                S = S[1:end - 1]

                if GFSECs(S, t)
                    con = @build_constraint(sum([x[i + 1, j + 1, t] * xsep[i + 1, j + 1, t] for i in S for j in S if i != j]) <= length(S) - sum(q[S, t]) / Q)
                    MOI.submit(mod, MOI.LazyConstraint(cb_data), con)
                    nbaddc = nbaddc + 1
                end
            end 

            for S in lS2
                if FCCs(S, t)
                    con = @build_constraint(sum([x[i, j + 1, t] for i in 1:n + 1 for j in S if !(i - 1 in S)]) >= sum(q[S, t]) / Q)
                    MOI.submit(mod, MOI.LazyConstraint(cb_data), con)
                    nbaddc2 = nbaddc2 + 1
                end
            end
        end
    end

    MOI.set(mod, MOI.LazyConstraintCallback(), lazySep)

    #println("Affichage du modèle avant résolution")
    print(mod)
    #println()

    #println("Résolution")
    optimize!(mod)
    #println()

    #println("Récupération et affichage \"à la main\" d'informations précises")
    status = termination_status(mod)

    if status == INFEASIBLE
        println("Le problème n'est pas réalisable")
    else
        #println("Affichage de tous les détails de la solution avec la commande solution_summary")
        println(solution_summary(mod, verbose = true))
        #println()

        res = []

        for t in 1:l
            resp = []

            for i in 1:n + 1
                for j in 1:n + 1
                    if i != j && abs(value(x[i, j, t]) - 1) < EPS
                        push!(resp, (i - 1, j - 1))
                    end
                end
            end

            push!(res, resp)
        end
    end

    return res, JuMP.objective_value(mod), JuMP.relative_gap(mod), nbaddc, nbaddc2
end


function Exact_ne(path, m, root_name = "Results")
    data = Read_data(path)

    if occursin("A_", path)
        type = "A"
    else
        type = "B"
    end

    res, robj, gap, nbaddc, nbaddc2 = Pdi_exact_ne(data, type, m)

    for i in 1:length(res)
        Draw_vrp(path, data, res[i], i, "exact", root_name)
    end

    return robj, gap
end


function Read_inst_pdi_exact_ne(s, m)
    ldir = readdir("./PRP_instances/")

    for i in 1:length(ldir)
        if occursin(s, ldir[i])
            x = split(ldir[i], ".")
            t0 = CPUtime_us()
            tmin, gap = Exact_ne(ldir[i], m)
            t1 = CPUtime_us()

            println(string("./Results/exact/PDI_exact_", x[1], "/", "PDI_exact_", x[1], ".log"))
            open(string("./Results/exact/PDI_exact_", x[1], "/", "PDI_exact_", x[1], ".log"), "w") do f
                write(f, string("t Execution time : ", (t1 - t0) / 1e6, " s\n"))
                write(f, string("o Objective value : ", tmin, "\n"))
                write(f, string("g Relative gap : ", round(gap * 100, digits = 2), "%\n"))
            end

            println(ldir[i], " checked")
        end
    end
end
