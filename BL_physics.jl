##########################################################################
# Compressible Blasius boundary layer with variable properties
# by Ludovico Fossa, 2024
#
# the code computes the self-similar solution to the 
# ordinary differential equations that govern a compressible flat-plate
# boundary layer with variable viscosity and thermal conductivity
# The code uses a standard block-elimination algorithm (see Cebeci, 2002)
# See the README for more details 
# 
# IMPORTANT
# to run the code open BL_main.jl and press CTRL+F5
# you can edit the parameters in the BL_input.jl file 
module  BL_physics
    export  mu,kcond
    #----------------------------------------------------------------------------------------------------------
    function mu(T) # shear viscosity
        T^(3/2)*(1+Main.chi_mu)/(T+Main.chi_mu)
    end
    #----------------------------------------------------------------------------------------------------------
    function kcond(T) # thermal conductivity
        T^(3/2)*(1+Main.chi_k)/(T+Main.chi_k)
    end
    #----------------------------------------------------------------------------------------------------------
end # end module
using .BL_physics
