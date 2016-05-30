include("two_link.jl")
using TwoLink
using DynamicMovementPrimitives
using Plots

cpoints = [0.5 -0.5;
        0.5 0;
        1.5 0;
        1.5 -0.5]

y,yd,ydd = connect_points(cpoints,20)
q = inverse_kin(y,:up)


# Setup ============================================
Nbasis = 20
αzC = 25.
αzPT = 10.
αx = 1.

opts = DMPopts(Nbasis,αx,αzC, :canonical)

T = size(y,1)
t = linspace(0,T,T)
τ = t[end]/3
h = t[2]-t[1]
qd = centraldiff(q) / h
qdd = centraldiff(qd) / h
g=q[end,:][:]
dmp = fit(q,qd,qdd,t,opts)


to,qo,zo,xo = solve(dmp)
yo = forward_kin(qo)


fig1 = plot(q, title="Joint space reference and output")
plot!(qo, l=:dash)

fig2 = plot(y, title="Cartesian space reference and output")
plot!(yo, l=:dash)

subplot(fig1,fig2); gui()
