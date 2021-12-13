include("heuris.jl")
include("exact.jl")


function Provide_full_results(s, m)

    """ String * Int -> Void

        s : Une sous chaine des noms d'instances a executer
        m : Le nombre maximal de vehicules

        Lance la resolution heuritique en deux phases et la resolution exacte (branch and cut) pour toute les instances ayant pour sous chaine s.
    """

    ldir = readdir("./PRP_instances/")

    for i in 1:length(ldir)
        if occursin(s, ldir[i])
            x = split(ldir[i], ".")
            t0 = CPUtime_us()
            tmin = Heuris(ldir[i], m, "Example_results")
            t1 = CPUtime_us()

            open(string("./Example_results/heu/PDI_heu_", x[1], "/", "PDI_heu_", x[1], ".log"), "w") do f
                write(f, string("t Execution time : ", (t1 - t0) / 1e6, " s\n"))
                write(f, string("o Objective value : ", tmin, "\n"))
            end

            println(ldir[i], " heuristic checked")
        end

        if occursin(s, ldir[i])
            x = split(ldir[i], ".")
            t0 = CPUtime_us()
            tmin, gap = Exact_ne(ldir[i], m, "Example_results")
            t1 = CPUtime_us()

            open(string("./Example_results/exact/PDI_exact_", x[1], "/", "PDI_exact_", x[1], ".log"), "w") do f
                write(f, string("t Execution time : ", (t1 - t0) / 1e6, " s\n"))
                write(f, string("o Objective value : ", tmin, "\n"))
                write(f, string("g Relative gap : ", round(gap * 100, digits = 2), "%\n"))
            end

            println(ldir[i], " optimal checked")
        end
    end
end


function Provide_class_results(s, m)

    """ String * Int -> Void

        s : Une sous chaine des noms d'instances a executer
        m : Le nombre maximal de vehicules

        Lance la resolution heuritique en deux phases et la resolution exacte (branch and cut) pour toute les instances par classe ayant pour sous chaine s.
    """

    ldir = readdir("./PRP_instances/")

    for i in 1:length(ldir)
        if occursin(s, ldir[i]) && (occursin("ABS12", ldir[i]) || occursin("ABS36", ldir[i]) || occursin("ABS60", ldir[i]) || occursin("ABS84", ldir[i]))
            x = split(ldir[i], ".")
            t0 = CPUtime_us()
            tmin = Heuris(ldir[i], m)
            t1 = CPUtime_us()

            open(string("./Results/heu/PDI_heu_", x[1], "/", "PDI_heu_", x[1], ".log"), "w") do f
                write(f, string("t Execution time : ", (t1 - t0) / 1e6, " s\n"))
                write(f, string("o Objective value : ", tmin, "\n"))
            end

            println(ldir[i], " heuristic checked")
        end

        if occursin(s, ldir[i]) && (occursin("ABS12", ldir[i]) || occursin("ABS36", ldir[i]) || occursin("ABS60", ldir[i]) || occursin("ABS84", ldir[i]))
            x = split(ldir[i], ".")
            t0 = CPUtime_us()
            tmin, gap = Exact_ne(ldir[i], m)
            t1 = CPUtime_us()

            open(string("./Results/exact/PDI_exact_", x[1], "/", "PDI_exact_", x[1], ".log"), "w") do f
                write(f, string("t Execution time : ", (t1 - t0) / 1e6, " s\n"))
                write(f, string("o Objective value : ", tmin, "\n"))
                write(f, string("g Relative gap : ", round(gap * 100, digits = 2), "%\n"))
            end

            println(ldir[i], " optimal checked")
        end
    end
end
