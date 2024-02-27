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
################################ BEGIN ###################################
using Printf
using Base.Filesystem
# clear all .dat files
script_dir = @__DIR__
all_files = readdir(script_dir)
dat_files = filter(file -> endswith(file, ".dat"), all_files)
for dat_file in dat_files
    file_path = joinpath(script_dir, dat_file)
    rm(file_path)
    println("File $file_path deleted successfully.")
end

# preliminary ops: read input and generate initial guess
# define global variables
global alpha0, alpha1, h, N_max, f, u, v, g, p, b, d, e
global GammaMatrix, DeltaMatrix
global error, k

include("BL_physics.jl")
include("BL_generate.jl")
include("BL_initial_guess.jl")
include("BL_matrices.jl")

# BLOCK ELIMINATION ALGORITHM - main loop
global k = 1 # iteration counter
global error = 1000.0
while (error>=tol)
    for j in 1:N_max
        b[j] = mu(g[j])/g[j]
        d[j] = (gamma-1.0)*Mach*Mach*b[j]
        e[j] = 1.0/Prandtl*kcond(g[j])/g[j]
    end

    invDeltaMatrix[:,:,1] = inv(buildA(1)) # builds inverse of matrix A at bottom boundary j = 0
    w[:,1] = buildr(1) # builds vector r at bottom boundary j = 0
    
    # upward sweep
    for j in 2:N_max 
        GammaMatrix[:,:] = buildB(j)*invDeltaMatrix[:,:,j-1] 
        DeltaMatrix[:,:] = buildA(j) - GammaMatrix[:,:]*buildC() 
        invDeltaMatrix[:,:,j] = inv(DeltaMatrix[:,:])
        w[:,j] = buildr(j) - GammaMatrix*w[:,j-1] 
    end

    delta[:,N_max] = invDeltaMatrix[:,:,N_max]*w[:,N_max]
    for j in (N_max-1):-1:1
        delta[:,j] = invDeltaMatrix[:,:,j]*(w[:,j] - buildC()*delta[:,j+1])    
    end

    global error = abs(delta[1,1])+abs(delta[2,1])+abs(delta[3,1])+abs(delta[4,1])+abs(delta[5,1])
    println("iteration $(@sprintf("%d",k)): error = $(@sprintf("%.20f",error))")
    global k+=1

    for j in 1:N_max
        f[j] = f[j] + delta[1,j]
        u[j] = u[j] + delta[2,j]
        v[j] = v[j] + delta[3,j]
        g[j] = g[j] + delta[4,j]
        p[j] = p[j] + delta[5,j]
    end
end

# PRINT THE RESULTS 
fname = "results_eta.dat"
f2 = open(fname, "w")
    println(f2,"eta     f       u       v       g       p")
    for j in 1:N_max
        println(f2,"$(@sprintf("%.5f", (j-1)*h)) $(@sprintf("%.5f", f[j])) $(@sprintf("%.5f", u[j])) $(@sprintf("%.5f", v[j])) $(@sprintf("%.5f", g[j])) $(@sprintf("%.5f", p[j]))")
    end
close(f2)
println("File $fname written successfully.")

fname = "results_y.dat"
y = zeros(Float64,N_max)
for j in 2:N_max
    y[j] = y[j-1] + sqrt(2)*g[j]*h
end   
f2 = open(fname, "w")
    println(f2,"y      f       u       v       g       p")
    for j in 1:N_max
        println(f2,"$(@sprintf("%.5f", y[j])) $(@sprintf("%.5f", f[j])) $(@sprintf("%.5f", u[j])) $(@sprintf("%.5f", v[j])) $(@sprintf("%.5f", g[j])) $(@sprintf("%.5f", p[j]))")
    end
close(f2)
println("File $fname written successfully.")

fname = "mach_y.dat"
mach = zeros(Float64,N_max)
for j in 1:N_max
    mach[j] = u[j]/sqrt(g[j])
end   
f2 = open(fname, "w")
    println(f2,"y      mach")
    for j in 1:N_max
        println(f2,"$(@sprintf("%.5f", y[j])) $(@sprintf("%.5f", mach[j]))")
    end
close(f2)
println("File $fname written successfully.")

fname = "rhoU_y.dat"
rhoU = zeros(Float64,N_max)
for j in 1:N_max
    rhoU[j] = u[j]/g[j]
end   
f2 = open(fname, "w")
    println(f2,"y      rhoU")
    for j in 1:N_max
        println(f2,"$(@sprintf("%.5f", y[j])) $(@sprintf("%.5f", rhoU[j]))")
    end
close(f2)
println("File $fname written successfully.")

################################ END ####################################