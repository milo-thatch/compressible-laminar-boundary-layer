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
########################################################################
############# generate initial guess for F and T #######################

println("BL_initial_guess.jl defines the initial guess")
println("BL_initial_guess.jl TwTad_ratio = $(@sprintf("%.2f", TwTad_ratio))")
if TwTad_ratio == 1.0
    alpha0 = 0.0
    alpha1 = 1.0
    println("BL_initial_guess.jl wall is ADIABATIC, alpha0 = $(@sprintf("%.2f", alpha0)), alpha1 = $(@sprintf("%.2f", alpha1))")
else
    alpha0 = 1.0
    alpha1 = 0.0
    println("BL_initial_guess.jl wall is ISOTHERMAL, alpha0 = $(@sprintf("%.2f", alpha0)), alpha1 = $(@sprintf("%.2f", alpha1))")
end
eta_ext = 5.0; # boundary layer thickness 
N_ext = trunc(Int, eta_ext/eta_max*N_max); # eta_max point
println("BL_initial_guess.jl eta_ext     = $(@sprintf("%.2f", eta_ext))")
println("BL_initial_guess.jl eta_max     = $(@sprintf("%.2f", eta_max))")
println("BL_initial_guess.jl N_max       = $N_max")
println("BL_initial_guess.jl N_ext       = $N_ext")

# initial guess boundary layer of thickness eta_max (Cebeci, 2002)
for j in 1:N_ext
    f[j] = eta_ext/4.0*((j-1)*h/eta_ext)^2*(3.0-((j-1)*h/eta_ext)^2/2.0)
    u[j] = 1.0/2.0*((j-1)*h/eta_ext)*(3.0-((j-1)*h/eta_ext)^2)
    v[j] = 3.0/(2.0*eta_ext)*(1.0-((j-1)*h/eta_ext)^2)
    if(alpha0==1.0)
        g[j] = T_w + (T_w-1.0)*j*h/eta_ext*((j-1)*h/eta_ext-2.0) 
        p[j] = 2*(T_w-1.0)/eta_ext*((j-1)*h/eta_ext-1.0)
    else
        g[j] = 1.0
        p[j] = 0.0
    end
end

# free-stream region
for j in (N_ext+1):N_max
    f[j] = eta_ext/4.0*(1.0)^2*(3.0-(1.0)^2/2.0) + (h*(j-1) - eta_ext)
    u[j] = 1.0
    v[j] = 0.0
    g[j] = 1.0
    p[j] = 0.0
end

# initial distributions of the viscous, conduction and dissipation factors (Cebeci, 2002)
for i in 1:N_max
    b[i] = mu(g[i])/g[i]
    d[i] = (gamma-1.0)*Mach*b[i]
    e[i] = 1.0/Prandtl*mu(g[i])/g[i]
end

# print initial conditions
fname =   "initial_conditions.dat"
f2 = open(fname, "w")
    println(f2,"eta     f       u       v       g       p        b       d       e")
    for j in 1:N_max
        println(f2,"$(@sprintf("%.5f", (j-1)*h)) $(@sprintf("%.5f", f[j])) $(@sprintf("%.5f", u[j])) $(@sprintf("%.5f", v[j])) $(@sprintf("%.5f", g[j])) $(@sprintf("%.5f", p[j])) $(@sprintf("%.5f", b[j])) $(@sprintf("%.5f", d[j])) $(@sprintf("%.5f", e[j]))")
    end
close(f2)