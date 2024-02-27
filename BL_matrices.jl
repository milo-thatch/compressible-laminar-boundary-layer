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
##########################################################################
# this module generate the submatrices A, B, and C (Cebeci, 2002)
##########################################################################
module BL_matrices
    using Printf
    export buildA, buildB, buildC, buildr
    import .. f, ..u, ..v, ..g, ..p, ..b, ..d, ..e
    #----------------------------------------------------------------------------------------------------------
    function buildA(n)
        A = zeros(Float64,5,5)
        # allocate elements in matrix A - bottom boundary
        if(n==1)
            A[1,1] = 1.0
            A[2,2] = 1.0
            A[3,4] = Main.alpha0
            A[3,5] = Main.alpha1
            A[4,2] = -1.0
            A[4,3] = -Main.h/2.0
            A[5,4] = -1.0
            A[5,5] = -Main.h/2.0
            return A
        end
        # compute coefficients s and beta
        s1 = b[n]/Main.h + f[n]/2.0
        s3 = v[n]/2.0
        s5 = 0.0
        s7 = 0.0
        beta1 = e[n]/Main.h + f[n]/2.0
        beta3 = p[n]/2.0
        beta5 = 0.0
        beta7 = 0.0
        beta9 = d[n]*v[n]
        # assign the coefficients to the elements in matrix A
        A[1,1] = 1.0
        A[1,2] = -Main.h/2.0
        A[2,1] = s3
        A[2,2] = s5
        A[2,3] = s1
        A[3,1] = beta3
        A[3,2] = beta5
        A[3,3] = beta9
        A[3,4] = beta7
        A[3,5] = beta1
        A[4,2] = -1.0
        A[4,3] = -Main.h/2.0
        A[5,4] = -1.0
        A[5,5] = -Main.h/2.0
        if(n==Main.N_max) # top boundary
            A[4,2] = 1.0
            A[4,3] = 0.0
            A[5,4] = 1.0
            A[5,5] = 0.0
        end
        return A
    end
    #----------------------------------------------------------------------------------------------------------
    function buildB(n)
        B = zeros(Float64,5,5)
        # compute coefficients s and beta
        s2 = - b[n-1]/Main.h + f[n-1]/2.0
        s4 = v[n-1]/2.0
        s6 = 0.0
        s8 = 0.0
        beta2 = - e[n-1]/Main.h + f[n-1]/2.0
        beta4 = p[n-1]/2.0
        beta6 = 0.0
        beta8 = 0.0
        beta10 = d[n-1]*v[n-1]
        # assign coefficients to the elements in matrix B
        B[1,1] = - 1.0
        B[1,2] = - Main.h/2.0
        B[2,1] = s4
        B[2,2] = s6
        B[2,3] = s2
        B[3,1] = beta4
        B[3,2] = beta6
        B[3,3] = beta10
        B[3,4] = beta8
        B[3,5] = beta2
        return B
    end
    function buildC()
        C = zeros(Float64,5,5)
        C[4,2] = 1.0
        C[4,3] = - Main.h/2.0
        C[5,4] = 1.0
        C[5,5] = - Main.h/2.0
        return C
    end
    function buildr(n)
        # compute residual vector
        r = zeros(Float64,5)
        if(n==1) # bottom boundary
            r[4] = u[1] - u[2] + Main.h/2.0*(v[2] + v[1])
            r[5] = g[1] - g[2] + Main.h/2.0*(p[2] + p[1])
            return r
        elseif(n==Main.N_max) # top boundary
            r[1] = f[n-1] - f[n] + Main.h/2.0*(u[n] + u[n-1])
            r[2] = - 1.0/Main.h*(b[n]*v[n] - b[n-1]*v[n-1]) - 0.5*(f[n]*v[n] + f[n-1]*v[n-1])
            r[3] = - 1.0/Main.h*(e[n]*p[n] - e[n-1]*p[n-1]) - 0.5*(f[n]*p[n] + f[n-1]*p[n-1]) - 0.5*(d[n]*v[n]*v[n] + d[n-1]*v[n-1]*v[n-1])
            return r
        end
        r[1] = f[n-1] - f[n] + Main.h/2.0*(u[n] + u[n-1])
        r[2] = - 1.0/Main.h*(b[n]*v[n] - b[n-1]*v[n-1]) - 0.5*(f[n]*v[n] + f[n-1]*v[n-1])
        r[3] = - 1.0/Main.h*(e[n]*p[n] - e[n-1]*p[n-1]) - 0.5*(f[n]*p[n] + f[n-1]*p[n-1]) - 0.5*(d[n]*v[n]*v[n] + d[n-1]*v[n-1]*v[n-1])
        r[4] = u[n] - u[n+1] + Main.h/2.0*(v[n+1] + v[n])
        r[5] = g[n] - g[n+1] + Main.h/2.0*(p[n+1] + p[n])
        return r
    end
    #----------------------------------------------------------------------------------------------------------
end # end module
using .BL_matrices