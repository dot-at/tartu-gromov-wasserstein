models = Dict{String,Function}()

using LightGraphs



models["ERDOS-RENY"] =
    function erdos_reny(args::Array{String,1}) :: Int
        local usg::String =
"""
Erdos-Reny
----------
Parameters (in this order):
  n        number of nodes
  p        edge probability
"""
        if isempty(args) || length(args)==1 && uppercase(args[2])=="HELP" || length(args) > 2
            println(usg)
            return 0
        end
        if !(  typeof(parse(args[2])) <: Integer && typeof(parse(args[3])) <: Real  )
            println(usg)
            return 0
        end
        local n::Integer = parse(args[2])
        local p::Real    = parse(args[3])
        


        return 0;
    end

models["STOCHASTIC_BLOCK_MODEL"] =
    function stochastic_block_model(args::Array{String,1}) :: Int
        if isempty(args) || length(args)==1 && uppercase(args[2])=="HELP"
            println(
"""
Stochastic Block Model
----------------------
Parameters (in this order):
  k        number of blocks
  p_0      edge probability inside blocks
  p_1      edge probability between blocks
  b_1      size of block #1
  ...
  b_k      size of block #k
"""
                    )
            return 0
        end
        # work

        return 0
        ;
    end

models["RANDOM_REGULAR"] =
    function random_regular(args::Array{String,1}) :: Int
        if isempty(args) || length(args)==1 && uppercase(args[2])=="HELP"
            println(
"""
Random Regular
--------------
Parameters (in this order):
  n        number of nodes
  d        degree
"""
                    )
            return 0
        end
        return 0;
    end

models["STATIC_FITNESS"] =
    function static_fitness(args::Array{String,1}) :: Int
        if isempty(args) || length(args)==1 && uppercase(args[2])=="HELP"
            println(
"""
Static Fitness
--------------
Parameters (in this order):
  ?        ?
"""
                    )
            return 0
        end
        return 0;
    end

models["SCALE_FREE"] =
    function scale_free(args::Array{String,1}) :: Int
        if isempty(args) || length(args)==1 && uppercase(args[2])=="HELP"
            println(
"""
Scale Free
----------
Parameters (in this order):
  ?        ?
"""
                    )
            return 0
        end
        return 0;
    end

models["BARABASI-ALBERT"] =
    function barabasi_albert(args::Array{String,1}) :: Int
        if isempty(args) || length(args)==1 && uppercase(args[2])=="HELP"
            println(
"""
Barabasi-Albert
---------------
Parameters (in this order):
  ?        ?
"""
                    )
            return 0
        end
        return 0;
    end

models["WATTS-STROGATZ"] =
    function watts_strogatz(args::Array{String,1}) :: Int
        if isempty(args) || length(args)==1 && uppercase(args[2])=="HELP"
            println(
"""
Watts-Strogatz
--------------
Parameters (in this order):
  ?        ?
"""
                    )
            return 0
        end
        return 0;
    end


function model_default(args::Array{String,1}) :: Int
    println("""
I don't recognize the model ``$(args[1])''
Try ``generate_random help'' ?
"""
            )
    return 1
    ;
end

usage = """
generate_random  ---
generate metric (measure) space from random network.
Usage:
  generate_random [help]            - this help text
  generate_random model [help]      - help for model
  generate_random model parameters  - do it.

``model'' is one of the following:
""" * string(keys(models))

function do_parameters(args::Array{String,1}) :: Int
    if isempty(args) || length(args)==1 && uppercase(args[1])=="HELP"
        println(usage)
        return 0
    else
        f = getindex(models, uppercase(args[1]), model_default)
        return f( uppercase(args[1]) )
    end
end

# main()

return do_parameters(ARGS)

;
#EOF
