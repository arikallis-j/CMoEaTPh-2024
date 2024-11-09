using Comonicon

include("src/func.jl")
include("src/manage.jl")

# #registration of packages
# include("src/data/gt/gt_fit.jl")

# fits = Dict(
#     "lfc" => fit_lfc, 
#     "gt" => fit_gt,
#     "be" => fit_be,
#     "fe" => fit_fe,
# )

const greeting = 
"""
# USAGE
The program for computing of the vychmaty
"""

const theoretical_notes = """
# THEORY 
* [data 13.03.2024]
Theory real...
* [data 31.10.2020]
"""



# realisations of schemes:
"""
Program CLI.

# Options
- `-e, --equation=<str>`: equation
    * equation = I.1

- `-n, --num=<str>`: num
    * num = 1 
    
- `-x, --x_n=<int>`: x
- `-t, --t_n=<int>`: t
- `--a=<int>`: a
- `--q=<int>`: q

# Flags 
- `-p, --print`: flag
"""
@cast function calc(;equation::String="I.1", num::Int=1, x_n::Int=100, t_n::Int=100, a::Float64 = 2.0, q::Float64 = 0.25, print::Bool=false)
    args = Dict(
        "x" => x_n,
        "t" => t_n,
        "a" => a,
        "q" => q,
    )
    calculate(equation, num, args) 
    if print
        println("drawing...")
    end
end


# """
# Program CLI.

# # Options
# - `-s, --str=<str>`: str
#     * str = abc
#     * str = efg
    
# - `--int=<int>`: int
# - `--flt=<int>`: flt

# # Flags 
# - `-f, --flag`: flag
# """
# @cast function calc(;str::String="lfc", int::Int=0, flt::Float64 = 0.0, flag::Bool=false)
#     args = Dict(
#         "str" => str,
#         "int" => int,
#         "flt" => flt,
#         "flag" => flag,
#     )
#     calculate(args) 
# end

"""
Program CLI.

# Flags 
- `--theory`: printing the theoretical notes
"""
@cast function guide(;theory::Bool=false)
    if theory
        println(theoretical_notes)
    else
        println(greeting)
    end
end

@main