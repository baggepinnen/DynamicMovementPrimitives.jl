using DynamicMovementPrimitives
using Base.Test
using Plots, DSP
math(sl) = map(s->string("\$",s,"\$") ,sl)
Nbasis = 15
az = 25.
ax = 1.#az/3
# bz = 1
# close("all")

opts = DMPopts(Nbasis,ax,az)
@test opts.βz == az/4


y = [zeros(10);linspace(0,1,100); ones(10)]
T = length(y)
t = linspace(0,T,T)
h = t[2]-t[1]
y = [y 0.5y]
ẏ = centraldiff(y) / h
ÿ = centraldiff(ẏ) / h

dmp = fit(y,ẏ,ÿ,opts)


tout,yout,ẏout,xout = solve(dmp,t)

plot(tout,[yout ẏout xout],lab = ["y_{out}" "ẏ_{out}" "x_{out}"] |> math)
plot!(t,[y ẏ],l=:dash,lab = ["y" "ẏ"] |> math)
gui()
# plot(force(dmp,solve_canonical(dmp,t)))
# gui()
