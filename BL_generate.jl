##########################################################################
# Compressible Blasius boundary layer with variable properties
# by Ludovico Fossa, 2024
#
# the code computes the self-similar solution to the 
# ordinary differential equations that govern a compressible flat-plate
# boundary layer with variable viscosity and thermal conductivity
# The code uses a standard block-elimination algorithm (see Cebeci, 2002)
# See the README and the attached .pdf for more technical details 
# 
# IMPORTANT
# to run the code open BL_main.jl and press CTRL+F5
# you can edit the parameters in the BL_input.jl file 
####################################################
################ input routines ####################
####################################################

println("BL_generate.jl reads the input")

include("BL_input.jl")

## OUTPUT
# generate the fields
f = Array{Float64, 1}(undef, N_max)
u = Array{Float64, 1}(undef, N_max)
v = Array{Float64, 1}(undef, N_max)
g = Array{Float64, 1}(undef, N_max)
p = Array{Float64, 1}(undef, N_max)
# generate the elements
b = Array{Float64, 1}(undef, N_max)
d = Array{Float64, 1}(undef, N_max)
e = Array{Float64, 1}(undef, N_max)
# arrays in the block elimination algorithm
w = Array{Float64, 2}(undef, 5, N_max)
delta = Array{Float64, 2}(undef, 5, N_max)
DeltaMatrix = Array{Float64, 2}(undef, 5, 5)
invDeltaMatrix = Array{Float64, 3}(undef, 5, 5, N_max)
GammaMatrix = Array{Float64, 2}(undef, 5, 5)

# discretization outputs
h = eta_max/N_max

# physical outputs
T_w = TwTad_ratio*(1+(gamma-1.0)/2.0*sqrt(Prandtl)*Mach^2)
println("BL_generate.jl T_w = $(@sprintf("%.2f", T_w))")