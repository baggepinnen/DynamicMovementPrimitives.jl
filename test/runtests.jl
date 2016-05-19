using DynamicMovementPrimitives
using Base.Test
import DynamicMovementPrimitives: get_sched_sig, get_centers_linear, kernel_matrix, _1

# Setup ============================================
Nbasis = 15
αz = 25.
αx = 1.#az/3
# bz = 1
# close("all")

optsC = DMPopts(Nbasis,αx,αz, :canonical)
optsT = DMPopts(Nbasis,αx,αz, :time)
optsP = DMPopts(Nbasis,αx,αz, :position)
@test optsC.βz == αz/4

y = [zeros(10);linspace(0,1,100); ones(10)]
T = length(y)
t = linspace(0,T,T)
τ = t[end]/3
h = t[2]-t[1]
y = [y 0.5y]
ẏ = centraldiff(y) / h
ÿ = centraldiff(ẏ) / h
g=y[end,:][:]

# Tests canonical ============================================

x = get_sched_sig(:canonical,αx,τ,t,y)
c,σ2    = get_centers_linear(Nbasis,x)
Ψ = kernel_matrix(x,c,σ2)
@test size(x) == (T,)
@test size(c) == (Nbasis,1)
@test size(σ2) == (Nbasis,1)
@test size(Ψ) == (T,Nbasis)

dmp = fit(y,ẏ,ÿ,t,optsC)
f = force(dmp,x)
@test size(f) == (T,2)
@test all(isfinite(f))
f = force(dmp,x[1],1)
@test isa(f,Number)
@test isfinite(f)

a = acceleration(dmp,y[1],ẏ[1],x[1],g[1],1)
@test isa(a,Number)
@test isfinite(a)
a = acceleration(dmp,_1(y),_1(ẏ),x[1],g)
@test size(a) == (2,)
a = acceleration(dmp,y,ẏ,x,g)
@test size(a) == (T,2)
@test all(isfinite(a))


# Tests time     ============================================
x = get_sched_sig(:time,αx,τ,t,y)
c,σ2    = get_centers_linear(Nbasis,x)
Ψ = kernel_matrix(x,c,σ2)
@test size(x) == (T,)
@test size(c) == (Nbasis,1)
@test size(σ2) == (Nbasis,1)
@test size(Ψ) == (T,Nbasis)

dmp = fit(y,ẏ,ÿ,t,optsT)
f = force(dmp,x)
@test size(f) == (T,2)
@test all(isfinite(f))
f = force(dmp,x[1],1)
@test isa(f,Number)
@test isfinite(f)

# Tests position ============================================
x = get_sched_sig(:position,αx,τ,t,y)
c,σ2    = get_centers_linear(Nbasis,x)
Ψ = kernel_matrix(x[:,1],c[:,1],σ2[:,1])
@test size(x) == (T,2)
@test size(c) == (Nbasis,2)
@test size(σ2) == (Nbasis,2)
@test size(Ψ) == (T,Nbasis)

dmp = fit(y,ẏ,ÿ,t,optsP)
f = force(dmp,x)
@test size(f) == (T,2)
@test all(isfinite(f))
f = force(dmp,x[1],1)
@test isa(f,Number)
@test isfinite(f)

a = acceleration(dmp,y[1],ẏ[1],x[1],g[1],1)
@test isa(a,Number)
@test isfinite(a)
a = acceleration(dmp,_1(y),_1(ẏ),x[1,:][:],g)
@test size(a) == (2,)
a = acceleration(dmp,y,ẏ,x,g)
@test size(a) == (T,2)
@test all(isfinite(a))





dmp = fit(y,ẏ,ÿ,t,optsC)
tout,yout,ẏout,xout = solve(dmp,t)
@test tout == t
@test abs(yout .- y) |> sum < 3
@test abs(ẏout .- ẏ) |> sum < 0.3
plotdmp(dmp)

dmp = fit(y,ẏ,ÿ,t,optsT)
tout,yout,ẏout,xout = solve(dmp,t)
@test abs(yout .- y) |> sum < 3
@test abs(ẏout .- ẏ) |> sum < 0.3
plotdmp(dmp)

dmp = fit(y,ẏ,ÿ,t,optsP)
tout,yout,ẏout,xout = solve(dmp,t)
@test abs(yout .- y) |> sum < 3
@test abs(ẏout .- ẏ) |> sum < 0.3
plotdmp(dmp)
