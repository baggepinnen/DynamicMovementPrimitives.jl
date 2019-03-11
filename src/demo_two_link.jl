include("two_link.jl")
using TwoLink
using DynamicMovementPrimitives
using Requires

cpoints = [0.5 -0.5;
0.5 0;
1.5 0;
1.5 -0.5]

y,yd,ydd = connect_points(cpoints,20)
q = inverse_kin(y,:up)


# Setup ============================================
Nbasis = 20
αz = 25.
αx = 1.

opts = DMPopts(Nbasis,αx,αz, :canonical)

T = size(y,1)
t = LinRange(0,T,T)
τ = t[end]/3
h = t[2]-t[1]
qd = centraldiff(q) / h
qdd = centraldiff(qd) / h
g=q[end,:][:]
dmp = fit(q,qd,qdd,t,opts)


to,qo,zo,xo = solve(dmp)
yo = forward_kin(qo)


@require Plots begin
fig1 = plot(q, xlabel="Time", title="Joint space reference and output", lab=["Reference_"].*["x" "y"])
plot!(qo, l=:dash, lab=["Fit_"].*["x" "y"])

fig2 = plot(y, xlabel="Time", title="Cartesian space reference and output", lab=["Reference_"].*["x" "y"])
plot!(yo, l=:dash, lab=["Fit_"].*["x" "y"])

fig3 = plot(q[:,1],q[:,2], xlabel="x", ylabel="y", title="Joint space reference and output", lab="Reference")
plot!(qo[:,1],qo[:,2], l=:dash, lab="Fit")

fig4 = plot(y[:,1],y[:,2], xlabel="x", ylabel="y", title="Cartesian space reference and output", lab="Reference")
plot!(yo[:,1],yo[:,2], l=:dash, lab="Fit")

subplot(fig1,fig2,fig3,fig4); gui()
end
