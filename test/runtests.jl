using DynamicMovementPrimitives
using Base.Test
using Plots
math(sl) = map(s->string("\$",s,"\$") ,sl)
Nbasis = 10
ax = 1
az = 1
# bz = 1

opts = DMPopts(Nbasis,ax,az)
@test opts.bz == 1/4

T = 100
t = linspace(0,T,T)
y = linspace(0,1,T)''
yd = centraldiff(y)
ydd = centraldiff(yd)

dmp = fit(y,yd,ydd,opts)


tout,yout,ydout,xout = solve(dmp)

plot(tout,[yout ydout xout],lab = ["y_{out}" "yd_{out}" "x_{out}"] |> math)
plot!(0:100,[y yd],l=:dash,lab = ["y" "y_d"] |> math)
gui()
