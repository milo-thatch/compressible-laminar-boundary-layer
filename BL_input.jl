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
##########################################################################
############################ INPUT FILE ##################################
##########################################################################

# discretization
eta_max     = 10.0          # maximum value of eta
N_max       = 1000          # discretization
tol         = 1e-12         # tolerance

# physical parameters
gamma       = 1.4
Prandtl     = 0.72
Mach        = 3.0
TwTad_ratio = 1.1
chi_mu      = 0.43          # non-dimensional sutherland constant
chi_k       = 0.66          # non-dimensional sutherland constant (conductivity)
