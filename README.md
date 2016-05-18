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
```julia
using DynamicMovementPrimitives
using Plots
Nbasis  = 15
αz      = 25.
αx      = 1.
opts    = DMPopts(Nbasis,αx,αz)

y       = [zeros(10);linspace(0,1,100); ones(10)]
T       = length(y)
t       = linspace(0,T,T)
h       = t[2]-t[1] # Sample interval
y       = [y 0.5y]
ẏ       = centraldiff(y) /h # Differentiate position to get velocity
ÿ       = centraldiff(ẏ) /h
dmp     = fit(y,ẏ,ÿ,t,opts)

tout,yout,ẏout,xout = solve(dmp,t) # Generate trajectory, see ?solve for options
plotdmp(dmp) # Requires Plots.jl, plots the trajectory from solve with default options
```
