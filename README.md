# DynamicMovementPrimitives

[![Build Status](https://travis-ci.org/baggepinnen/DynamicMovementPrimitives.jl.svg?branch=master)](https://travis-ci.org/baggepinnen/DynamicMovementPrimitives.jl)

Provides an implementation of Ijspeert et al. 2013

## Installation

```julia
Pkg.add("DynamicMovementPrimitives")
Pkg.checkout("DynamicMovementPrimitives") # Recommended to get the latest version
using DynamicMovementPrimitives
```

## Usage
### Standard DMP
```julia
using DynamicMovementPrimitives
Nbasis  = 15
αz      = 25.
αx      = 1.
opts    = DMPopts(Nbasis,αx,αz)

y       = [zeros(10);linspace(0,2,1000); 2ones(10)]
T       = length(y)
t       = linspace(0,10,T)
h       = t[2]-t[1] # Sample interval
y       = [y 0.5y]
ẏ       = centraldiff(y) /h # Differentiate position to get velocity
ÿ       = centraldiff(ẏ) /h
dmp     = fit(y,ẏ,ÿ,t,opts)

tout,yout,ẏout,xout = solve(dmp,t) # Generate trajectory, see ?solve for options
plot(dmp) # Requires Plots.jl, plots the trajectory from solve with default options
plot(dmp,true)
```

### DMP with two degrees of freedom (Karlsson, Bagge Carlsson et al. 2017)
We start by upgrading the DMP object to incorporate also the controller parameters for the 2DOF controller
```julia
dmp2opts = DMP2dofopts(kp = 25,kv = 10,kc = 10_000,αe = 5)
dmp2 = DMP2dof(dmp, dmp2opts) # Upgrade dmp to 2DOF version

t,yc,ẏc,x,ya,ẏa,e = solve(dmp2,t)
plot(dmp2) # Requires Plots.jl, plots the trajectory from solve with default options
plot(dmp2,true)
```

We test the performance of the 2DOF controller by implementing a custom solver. The solver `euler_disturbance` is the same as the built in solver `euler`, but between `i=500` and `i=700`, we stop the evolution of the physical system by setting `ẏa = 0`.
```julia
function euler_disturbance(time_derivative, state0, t, args...; kwargs...)
    T = length(t)
    n = length(state0)
    res = Matrix{Float64}(T,n)
    res[1,:] = state0
    for i in 2:T
        td = time_derivative(t[i-1],res[i-1,:])*(t[i]-t[i-1])
        if 500 <= i < 700 # Hold ẏa still
            td[4] = 0
        end
        res[i,:] = res[i-1,:] + td
    end
    t,res
end
```

We can call the solve method with our custom solver and plot the result. It should be clear from the figures that this time, the coupled signal `yc` slows down when there is a nonzero error.
```julia
t,yc,ẏc,x,ya,ẏa,e = solve(dmp2,t, solver=euler_disturbance)
plot(t,ẏc, lab="\$ẏ_c\$", c=:red, l=(:dash, 3), layout=(2,2), subplot=1)
plot!(t,yc, lab="\$y_c\$", c=:red, l=(:dash, 3), subplot=2)
plot!(t,ẏa, lab="\$ẏ_a\$", c=:blue, subplot=1)
plot!(t,ya, lab="\$y_a\$", c=:blue, subplot=2)
plot!(t,e, lab="\$e\$", c=:green, subplot=3)
plot!(t,500 .<= 1:T .< 700, lab="Disturbance", c=:green, subplot=4, fillrange=0)
t,yc,ẏc,x,ya,ẏa,e = solve(dmp2,t, solver=euler)
plot!(t,ẏc, lab="\$ẏ_u\$", c=:black, l=(:dashdot, 3), subplot=1)
plot!(t,yc, lab="\$y_u\$", c=:black, l=(:dashdot, 3), subplot=2)
```
In the figure below, the black line represents the evolution with no disturbance, in the paper referred to as the unperturbed evolution. The blue evolution is the actual system evolution whereas the red curve displays the coupled system evolution.
![DMP2dof plot](figs/dmp2.png)
