# DynamicMovementPrimitives

[![Build Status](https://travis-ci.org/baggepinnen/DynamicMovementPrimitives.jl.svg?branch=master)](https://travis-ci.org/baggepinnen/DynamicMovementPrimitives.jl)

## Installation

`Pkg.clone("git@github.com:baggepinnen/DynamicMovementPrimitives.jl.git")`
`using DynamicMovementPrimitives`

## Usage
```julia
using DynamicMovementPrimitives
using Plots
math(sl) = map(s->string("\$",s,"\$") ,sl)
Nbasis  = 15
az      = 25.
ax      = 1.
opts    = DMPopts(Nbasis,ax,az)

y       = [zeros(10);linspace(0,1,100); ones(10)]
T       = length(y)
t       = linspace(0,T,T)
h       = t[2]-t[1]
y       = [y 0.5y]
ẏ       = centraldiff(y) / h
ÿ       = centraldiff(ẏ) / h
dmp     = fit(y,ẏ,ÿ,opts)

tout,yout,ẏout,xout = solve(dmp,t)

plot(tout,[yout ẏout xout],lab = ["y_{out}" "ẏ_{out}" "x_{out}"] |> math)
plot!(t,[y ẏ],l=:dash,lab = ["y" "ẏ"] |> math)
```
