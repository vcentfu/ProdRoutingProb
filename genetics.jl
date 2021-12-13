using JuMP
using CPLEX
using Random
include("utils.jl")
include("lsp.jl")


const EPS = 1e-10
const OPTIMAL = JuMP.MathOptInterface.OPTIMAL
const INFEASIBLE = JuMP.MathOptInterface.INFEASIBLE
const UNBOUNDED = JuMP.MathOptInterface.DUAL_INFEASIBLE;


function Init_vrp(Nc, t)
    vals  = [i for i in 1:length(Nc[:, t]) if Nc[i, t] == 1]
    inds = shuffle!([i for i in 1:length(vals)])

    return inds, vals
end


function Meta_vrp(res)
    return [k for k in res if k != 0]
end


function Inci_vrp(res)
    inci_res = [0 for i in 1:maximum(res)]
    tour = 1

    for k in res[2:end]
        if k != 0
            inci_res[k] = tour
        else
            tour = tour + 1
        end
    end

    return inci_res
end


function Classic_crossover(p1, p2)
    e1 = copy(p1)
    e2 = copy(p2)
    cs = div(length(p1), 3)
    ss = div(length(p1), 2)

    for i in 1:cs
        e1[i] = p2[i]
        e2[i] = p1[i]
    end

    for i in (cs + ss + 1):length(p1)
        e1[i] = p2[i]
        e2[i] = p1[i]
    end

    return e1, e2
end


function Order_crossover(p1, p2)
    e1 = copy(p1)
    e2 = copy(p2)
    cs = div(length(p1), 3)
    ss = div(length(p1), 2)
    o1 = append!([p2[i] for i in (cs + ss + 1):length(p2)],[p2[i] for i in 1:(cs + ss)])
    o2 = append!([p1[i] for i in (cs + ss + 1):length(p2)],[p1[i] for i in 1:(cs + ss)])
    s1 = [p1[i]  for i in (cs + 1):(cs + ss)]
    s2 = [p2[i]  for i in (cs + 1):(cs + ss)]
    c = 0
    i = 0

    while c != length(p1) - ss
        i = i + 1

        if o1[i] in s1
            continue
        else
            e1[((cs + ss) + c) % length(p1) + 1] = o1[i]
            c = c + 1
        end
    end

    c = 0
    i = 0

    while c != length(p2) - ss
        i = i + 1

        if o2[i] in s2
            continue
        else
            e2[((cs + ss) + c) % length(p2) + 1] = o2[i]
            c = c + 1
        end
    end

    return e1, e2
end


function Mate(i1)
    res = copy(i1)
    r = rand([i for i in 1:length(i1)])
    rv = [i for i in 1: length(i1) if !(i in [r])]

    if length(rv) != 0
        s = rand(rv)
        temp = i1[r]
        res[r] = i1[s]
        res[s] = temp
    end

    return res
end


function Seperation(data, i1, v1, m)
    res = [0]
    c = 0
    cm = m - 1

    for i in 1:length(i1)
        if c + data["d"][v1[i1[i]]] > data["Q"] && cm > 0
            push!(res, 0)
            c = 0
            cm = cm - 1
        end

        if cm == 0
            #println([i for i in 1:length(res) if res[i] == 0])
            println("pas assez de voitures")
            return []
        end

        push!(res, i1[i])
        c = c + data["d"][v1[i1[i]]]
    end

    push!(res, 0)

    return res
end


function Distance_vrp(data, type, i1, v1)
    res = 0

    for i in 2:length(i1)
        if type == "A"
            if i1[i - 1] == 0
                res = res + Distance_A(data, 1, v1[i1[i]] + 1)
            elseif i1[i] == 0
                res = res + Distance_A(data, v1[i1[i - 1]] + 1, 1)
            else
                res = res + Distance_A(data, v1[i1[i - 1]] + 1, v1[i1[i]] + 1)
            end
        end

        if type == "B"
            if i1[i - 1] == 0
                res = res + Distance_B(data, 1, v1[i1[i]] + 1)
            elseif i1[i] == 0
                res = res + Distance_B(data, v1[i1[i - 1]] + 1, 1)
            else
                res = res + Distance_B(data, v1[i1[i - 1]] + 1, v1[i1[i]] + 1)
            end
        end
    end

    return res
end


function Computes_r(vrp)
    r = length([k for k in vrp if k == 0]) - 1

    return r
end



function Fitness(data, type, vrp, v1, rmin, dmin)
    r = Computes_r(vrp)
    dvrp = Distance_vrp(data, type, vrp, v1)

    return r - rmin + min(dvrp, 2 * dmin) / dmin
end


function Launch_GA(data, type, Nc, m, t, lamb, mu, ngen)
    v1 = Init_vrp(Nc, t)[2]
    population = [Init_vrp(Nc, t)[1] for i in 1:lamb]
    population = [ind for ind in population if ind != []]
    population = [Seperation(data, population[i], v1, m) for i in 1:length(population)]

    if length(population) == 0
        println("Pas assez de voitures, pas besoin de voitures")
        return [], 0
    end

    pr = []
    pd = []

    for j in 1:ngen
        offspring = []
        for i in 1:div(lamb, 2)

            p1 = rand(population)
            p2 = rand(population)

            t1 = Meta_vrp(p1)
            t2 = Meta_vrp(p2)

            while length(t1) == 0 || length(t2) == 0
                p1 = rand(population)
                p2 = rand(population)

                t1 = Meta_vrp(p1)
                t2 = Meta_vrp(p2)
            end

            e1, e2 = Order_crossover(t1, t2)

            if rand() < 0.5
                e1 = Mate(e1)
            end

            if rand() < 0.5
                e2 = Mate(e2)
            end

            e1 = Seperation(data, e1, v1, m)
            e2 = Seperation(data, e2, v1, m)
            push!(offspring, e1)
            push!(offspring, e2)
        end

        offspring = append!(population, offspring)
        pr = [Computes_r(offspring[i]) for i in 1:length(offspring)]
        pd = [Distance_vrp(data, type, offspring[i], v1) for i in 1:length(offspring)]
        imin = argmin(pd)
        rmin = pr[imin]
        dmin = pd[imin]
        pFitness = [Fitness(data, type, offspring[i], v1, rmin, dmin) for i in 1:length(offspring)]
        best = sort([(pFitness[i], i) for i in 1:length(offspring) if pFitness[i] != 0])
        population = []

        for i in 1:mu
            push!(population, offspring[best[i][2]])

            if i == length(best)
                break
            end
        end

        population = [ind for ind in population if ind != []]

    end

    pd = [Distance_vrp(data, type, population[i], v1) for i in 1:length(population)]

    res = []

    for i in 2:length(population[1])
        push!(res, (population[1][i - 1], population[1][i]))
    end

    return res, pd[1]
end



"""
data = Read_data("A_014_ABS96_15_5.prp")
SC = Initialize_SC(data, "A")
Nc, = Lsp_solv_heu(data, SC)
t1, te1 = Init_vrp(Nc, 2)
t2, = Init_vrp(Nc, 2)
t1 = Meta_vrp(t1)
t2 = Meta_vrp(t2)
e1, e2 = Order_crossover(t1, t2)
e3 = Mate(e1)
r4 = Seperation(data, e3, te1, 10)
Distance_vrp(data, "A", r4, te1)
x1, x2 = Launch_GA(data, "A", Nc, 10, 2, 200, 100, 200)"""
