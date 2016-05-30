using DynamicMovementPrimitives
using Base.Test
import DynamicMovementPrimitives: get_sched_sig, get_centers_linear, kernel_matrix, _1

# Setup ============================================
Nbasis = 20
αzC = 25.
αzPT = 10.
αx = 1.#αz/3
# bz = 1
# close("all")

optsC = DMPopts(Nbasis,αx,αzC, :canonical)
optsT = DMPopts(Nbasis,αx,αzPT, :time)
optsP = DMPopts(Nbasis,αx,αzPT, :position)
@test optsC.βz == αzC/4

y = [linspace(0,1,100); ones(10)]
T = length(y)
t = linspace(0,T,T)
τ = t[end]/3
h = t[2]-t[1]
y = [y 0.5y]
ẏ = centraldiff(y) / h
ÿ = centraldiff(ẏ) / h
g=y[end,:][:]

# Tests canonical ============================================

x = get_sched_sig(:canonical,αx,τ,t,y,g)
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
x = get_sched_sig(:time,αx,τ,t,y,g)
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
x = get_sched_sig(:position,αx,τ,t,y,g)
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
tout,youtC,ẏout,xout = solve(dmp)
@test tout == t
@test abs(youtC .- y) |> sum < 2
@test abs(ẏout .- ẏ) |> sum < 0.3
# plotdmp(dmp)

dmp = fit(y,ẏ,ÿ,t,optsT)
tout,youtT,ẏout,xout = solve(dmp)
@test abs(youtT .- y) |> sum < 2
@test abs(ẏout .- ẏ) |> sum < 0.3
# plotdmp(dmp)

dmp = fit(y,ẏ,ÿ,t,optsP)
tout,youtP,ẏout,xout = solve(dmp)
@test abs(youtP .- y) |> sum < 6
@test abs(ẏout .- ẏ) |> sum < 0.3
# plotdmp(dmp)
# plotdmpphase(dmp)





# Test twolink
include("../src/two_link.jl")
using TwoLink
p,pd,pdd = traj(0,1,0:100)
@test p[1] == 0
@test p[end] == 1
@test pd[1] == 0
@test pd[end] == 0
@test length(p) == length(pd) == length(pdd) == 101




cpoints = [0.5 -0.5;
        0.5 0;
        1.5 0;
        1.5 -0.5]

ctraj = connect_points(cpoints,20)[1]
jtraj = inverse_kin(ctraj,:up)
ctraj2 = forward_kin(jtraj)

@test ctraj ≈ ctraj2
