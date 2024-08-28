mutable struct DMP2dofopts
    kp::Float64
    kv::Float64
    kc::Float64
    αe::Float64
end

DMP2dofopts(;kp = 25,kv = 10,kc = 10_000,αe = 5) = DMP2dofopts(kp,kv,kc,αe)

"""
Same as DMP but contains an extra struct `opts2` with 2DOF controller parameters
Upgrade a `DMP` to a `DMP2dof` using

 ```
 dmp2opts = DMP2dofopts(kp = 25,kv = 10,kc = 10_000,αe = 5) # Specify parameters here\n
 dmp2 = DMP2dof(dmp, dmp2opts) # Upgrade dmp to 2DOF version
 ```
"""
mutable struct DMP2dof <: AbstractDMP
    opts::DMPopts
    g::Vector{Float64}
    y::Matrix{Float64}
    ẏ::Matrix{Float64}
    ÿ::Matrix{Float64}
    t::AbstractVector
    w::Matrix{Float64}
    τ::Float64
    c::VecOrMat{Float64}
    σ2::VecOrMat{Float64}
    opts2::DMP2dofopts
end

DMP2dof(d::DMP, dmp2opts::DMP2dofopts) = DMP2dof(d.opts,d.g,d.y,d.ẏ,d.ÿ,d.t,d.w,d.τ,d.c,d.σ2,dmp2opts)

function acceleration(d::DMP2dof, yc::Number,ẏc::Number,x::Number,ya,ẏa,e,g::Number,i=1)
    f = force(d,x,i)
    _acceleration(d,f,yc,ẏc,x,ya,ẏa,e,g)
end

function acceleration(d::DMP2dof, yc,ẏc,x,ya,ẏa,e,g=d.g)
    f = force(d,x)
    _acceleration(d,f,yc,ẏc,x,ya,ẏa,e,g)
end

function _acceleration(d,f, yc,ẏc,x,ya,ẏa,e,g)
    αx,αz,βz,αe,τ,kc,kp,kv = d.opts.αx,d.opts.αz,d.opts.βz,d.opts2.αe,d.τ,d.opts2.kc,d.opts2.kp,d.opts2.kv
    τa  = τ*(1+kc*e^2)
    ẋ   = -αx/τa * x
    z   = τa*ẏc
    ż   = (αz*(βz*(g-yc)-z)+f)/τa
    ė   = αe*(ya-yc-e)
    ÿc  = (ż*τa - 2τ*kc*z*e*ė)/τa^2
    ÿa  = kp*(yc-ya) + kv*(ẏc-ẏa) + ÿc
    ẏc,ÿc,ẏa,ÿa,ė,ẋ
end

function solve_canonical(dmp::DMP2dof, t, y0, g, solver; kwargs...)
    T,n = size(dmp.y)
    αx  = dmp.opts.αx
    τ   = dmp.τ
    yc  = zeros(T,n)
    ẏc  = zeros(T,n)
    ya  = zeros(T,n)
    ẏa  = zeros(T,n)
    e   = zeros(T,n)
    x   = zeros(T)
    for i = 1:n
        function time_derivative(dstate, state, p, t)
            local yc  = state[1]
            local ẏc  = state[2]
            local ya  = state[3]
            local ẏa  = state[4]
            local e   = state[5]
            local x   = state[6]
            ẏc,ÿc,ẏa,ÿa,ė,ẋ = acceleration(dmp,yc,ẏc,x,ya,ẏa,e,g[i],i)
            dstate[1] = ẏc
            dstate[2] = ÿc
            dstate[3] = ẏa
            dstate[4] = ÿa
            dstate[5] = ė
            dstate[6] = ẋ
        end
        state0 = [y0[i], dmp.ẏ[1,i], y0[i], dmp.ẏ[1,i], 0, 1.]
        prob = OrdinaryDiffEq.ODEProblem(time_derivative,state0,(t[1],t[end]))
        sol = OrdinaryDiffEq.solve(prob,solver;saveat=t,dt=t[2], kwargs...)
        yc[:,i]  .= sol(t,idxs=1)
        ẏc[:,i]  .= sol(t,idxs=2)
        ya[:,i]  .= sol(t,idxs=3)
        ẏa[:,i]  .= sol(t,idxs=4)
        e[:,i]   .= sol(t,idxs=5)
        x        .= sol(t,idxs=6) # TODO: se till att denna är samma för alla DOF
    end
    t,yc,ẏc,x,ya,ẏa,e
end
