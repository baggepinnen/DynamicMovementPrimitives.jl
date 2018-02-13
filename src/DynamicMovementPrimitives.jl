module DynamicMovementPrimitives
using ODE, Requires, RecipesBase


export DMP, DMPopts, centraldiff,fit, solve, force, acceleration, solve_canonical, kernel_vector, plotdmp, plotdmpphase, euler

export DMP2dofopts, DMP2dof

"""
`DMPopts(Nbasis,αx,αz) = DMPopts(Nbasis, αx, αz, βz = αz/4)`\n
Holds parameters for fitting a DMP
# Fields
`Nbasis,αx,αz,βz,sched_sig,fitmethod`\n
`sched_sig` can be chosen as `:canonical` (default), `:time` or `position`\n
`fitmethod` can be chosen as `:lwr` or `:leastsquares`

See example file or the paper by Ijspeert et al. 2013
"""
struct DMPopts
    Nbasis::Int
    αx::Float64
    αz::Float64
    βz::Float64
    sched_sig::Symbol
    fitmethod::Symbol
end

DMPopts(Nbasis,αx,αz,sched_sig::Symbol=:canonical,fitmethod::Symbol = :lwr) = DMPopts(Nbasis,αx,αz,αz/4,sched_sig,fitmethod)

abstract type AbstractDMP end

"""
The result of fitting a DMP
#Fields
`opts,g,y,ẏ,ÿ,w,τ,c,σ2`\n
See example file or the paper by Ijspeert et al. 2013
"""
mutable struct DMP <: AbstractDMP
    opts::DMPopts
    g::Vector{Float64}
    y::Matrix{Float64}
    ẏ::Matrix{Float64}
    ÿ::Matrix{Float64}
    t::AbstractVector
    w::Matrix
    τ::Float64
    c::VecOrMat{Float64}
    σ2::VecOrMat{Float64}
end

include("utilities.jl")

"""
`solve_canonical(αx,τ,T::AbstractVector)`\n
`solve_canonical(αx,τ,T::Real)`\n
`solve_canonical(dmp::DMP [,t])`
"""
solve_canonical(αx,τ,T::AbstractVector) = exp.(-αx/τ.*T)
solve_canonical(αx,τ,T::Real)           = solve_canonical(αx,τ,(0:T-1))
solve_canonical(dmp::DMP,t)             = solve_canonical(dmp.opts.αx, dmp.τ, t)
solve_canonical(dmp::DMP)               = solve_canonical(dmp.opts.αx, dmp.τ, dmp.t)


function get_sched_sig(s,αx,τ,t,y,g)
    if s == :canonical
        return solve_canonical(αx,τ,t)
    elseif s == :position
        return g'.-y
    elseif s == :time
        return t
    end
    warn("Scheduling signal $s unknown")
    return solve_canonical(αx,τ,t)
end
get_sched_sig(dmp::DMP) = get_sched_sig(dmp.opts.sched_sig,dmp.opts.αx,dmp.τ,dmp.t,dmp.g)


"""
`fit(y, ẏ, ÿ, opts, g=y[end])`

Fits a DMP to data\n
`y, ẏ, ÿ` are position, velocity and acceleration respectively, `T×n` matrices where `T` is the number of time steps and `n` is the number of degrees of freedom.

see also `solve`, `plotdmp`
"""
function fit(y,ẏ,ÿ,t,opts,g=y[end,:][:])
    if opts.sched_sig != :canonical
        error("Scheduling signal $(opts.sched_sig) currently not supported")
    end
    T       = t[end]
    N       = size(y,1)
    n       = isa(y,Matrix) ? size(y,2) : 1
    τ       = T/3 # After three time constants we have reached 1-exp(-3) ≈ 0.95
    Nbasis  = opts.Nbasis
    βz      = opts.βz
    αz      = opts.αz
    αx      = opts.αx
    y0      = _1(y)
    x       = get_sched_sig(opts.sched_sig,αx,τ,t,y,g)
    ft      = τ^2*ÿ .- αz*(βz*(g'.-y).-τ*ẏ)
    ξ       = x.*(g-y0)'
    c,σ2    = get_centers_linear(Nbasis,x)
    if opts.sched_sig != :position
        Ψ = kernel_matrix(x,c,σ2)
    end
    w = zeros(Nbasis,n)
    for i = 1:n #joints
        if opts.sched_sig == :position
            Ψ = kernel_matrix(x[:,i],c[:,i],σ2[:,i])
        end
        if opts.fitmethod == :leastsquares
            w = Ψ\(ft./ξ)
        else # LWR
            for j = 1:Nbasis
                sTΓ = ξ[:,i].*Ψ[:,j]
                w[j,i] = vecdot(sTΓ,ft[:,i])/vecdot(sTΓ,ξ[:,i])
            end
        end
    end
    return DMP(opts, g, y, ẏ, ÿ,t, w, τ,c,σ2)
end

function force_single(d::AbstractDMP,x::Number, i)
    # ODE Point case
    y0 = _1(d)
    Ψ  = kernel_vector(x, d.c, d.σ2)
    f = vecdot(Ψ,d.w[:,i]) * x*(d.g[i]-y0[i])
end

function force_single(d::AbstractDMP,x::Number)
    # Point case
    y0 = _1(d)
    Ψ  = kernel_matrix([x], d.c, d.σ2)
    f  = (Ψ*d.w)[:] .* (x*(d.g-y0))
end

function force_single(d::AbstractDMP,x::AbstractVector)
    # Trajectory case
    y0 = _1(d)
    Ψ    = kernel_matrix(x, d.c, d.σ2)
    f = Ψ*d.w .* (x.*(d.g-y0)')
end

function force_multiple(d::AbstractDMP,x::Number, i)
    # ODE Point case
    y0 = _1(d)
    Ψ  = kernel_vector(x, d.c[:,i], d.σ2[:,i])
    f  = vecdot(Ψ,d.w[:,i]) * x*(d.g[i]-y0[i])
end

function force_multiple(d::AbstractDMP,x::AbstractVector)
    # Point case
    n = size(d.c,2) # Number of DOF
    y0 = _1(d)
    f = zeros(n)
    for i = 1:n
        Ψ    = kernel_matrix([x[i]], d.c[:,i], d.σ2[:,i])
        f[i] = vecdot(Ψ,d.w[:,i]) * x[i].*(d.g[i]-y0[i])
    end
    f
end

function force_multiple(d::AbstractDMP,x::AbstractMatrix)
    # Trajectory case
    T,n = size(x)
    y0 = _1(d)
    f = zeros(T,n)
    for i = 1:n
        Ψ    = kernel_matrix(x[:,i], d.c[:,i], d.σ2[:,i])
        f[:,i] = Ψ*d.w[:,i] .* (x[:,i].*(d.g[i]-y0[i])')
    end
    f
end

"""
`force(dmp,x)`

Calculate the forcing term for `dmp` when the phase variable is `x`\n
The return value will be an `n` Vector or `T×n` Matrix depending on `typeof(x)`
"""
function force(d::AbstractDMP, x)
    if d.opts.sched_sig == :position
        force_multiple(d,x)
    else
        force_single(d,x)
    end
end

function force(d::AbstractDMP, x,i)
    if d.opts.sched_sig == :position
        force_multiple(d,x,i)
    else
        force_single(d,x,i)
    end
end

function acceleration(d::DMP, y::Number,ẏ::Number,x::Number,g::Number,i=1)
    f = force(d,x,i)
    (d.opts.αz*(d.opts.βz*(g-y)-d.τ*ẏ)+f)/d.τ^2
end

function acceleration(d::DMP, y::AbstractVector,ẏ::AbstractVector,x ,g::AbstractVector = d.g)
    f = force(d,x)
    (d.opts.αz*(d.opts.βz*(g-y)-d.τ*ẏ)+f)/d.τ^2
end

function acceleration(d::DMP, y::AbstractMatrix,ẏ::AbstractMatrix,x::AbstractVecOrMat,g::AbstractVector)
    f = force(d,x)
    (d.opts.αz*(d.opts.βz*(g'.-y)-d.τ*ẏ)+f)/d.τ^2
end

function solve_canonical(dmp::DMP, t, y0, g, solver)
    T,n = size(dmp.y)
    αx  = dmp.opts.αx
    τ   = dmp.τ
    z   = zeros(T,n)
    y   = zeros(T,n)
    x   = zeros(T)
    for i = 1:n
        function time_derivative(t,state)
            local z   = state[1]
            local y   = state[2]
            local x   = state[3]
            zp  = acceleration(dmp, y, z, x,g[i],i)
            yp  = z
            xp  = -αx/τ * x
            [zp;yp;xp]
        end
        state0 = [dmp.ẏ[1,i]; y0[i]; 1.]
        tout,state_history = solver(time_derivative, state0, t,points=:specified)
        res    = vv2m(state_history)
        z[:,i] = res[:,1]
        y[:,i] = res[:,2]
        x      = res[:,3] # TODO: se till att denna är samma för alla DOF
    end
    t,y,z,x
end

function solve_position(dmp, t, y0, g, solver)
    T,n = size(dmp.y)
    αx  = dmp.opts.αx
    τ   = dmp.τ
    z   = zeros(T,n)
    y   = zeros(T,n)
    for i = 1:n
        function time_derivative(t,state)
            local z   = state[1]
            local y   = state[2]
            zp  = acceleration(dmp, y, z, g[i]-y,g[i],i)
            yp  = z
            [zp;yp]
        end
        state0  = [0; y0[i]]
        tout,state_history = solver(time_derivative, state0, t,points=:specified)
        res = vv2m(state_history)
        z[:,i] = res[:,1]
        y[:,i] = res[:,2]
    end
    t,y,z,y
end

function solve_time(dmp, t, y0, g, solver)
    T,n     = size(dmp.y)
    αx      = dmp.opts.αx
    τ       = dmp.τ
    z       = zeros(T,n)
    y       = zeros(T,n)
    for i = 1:n
        function time_derivative(t,state)
            local z   = state[1]
            local y   = state[2]
            zp  = acceleration(dmp, y, z, t ,g[i],i)
            yp  = z
            [zp;yp]
        end
        state0  = [0; y0[i]]
        tout,state_history = solver(time_derivative, state0, t,points=:specified)
        res = vv2m(state_history)
        z[:,i] = res[:,1]
        y[:,i] = res[:,2]
    end
    t,y,z,t
end

"""
`t,y,z,x = solve(dmp, t = 0:length(dmp.t)-1; y0 = _1(dmp), g = dmp.g, solver=ode45)`

`t` time vector

## Keyword arguments: \n
`y0` start position, defaults to the initial point in training data from `dmp`
`g` goal, defaults to goal from `dmp`\n
`solver` the ode solver to use, see https://github.com/JuliaLang/ODE.jl \n
The default solver is `solver=ode54`, a faster alternative is `solver=euler`

see also `plotdmp`
"""
function solve(dmp::AbstractDMP, t = dmp.t; y0 = _1(dmp), g = dmp.g, solver=ode45)
    if dmp.opts.sched_sig == :position
        return solve_position(dmp, t, y0, g, solver)
    elseif dmp.opts.sched_sig == :time
        return solve_time(dmp, t, y0, g, solver)
    end
    solve_canonical(dmp, t, y0, g, solver)
end


"""
    A = linearize(dmp::AbstractDMP, y, ẏ, x, y0)

Linearize the dynamics of the DMP around the point (y,ẏ,x). The resulting linear system is on the form
`q̈ = A*[y; ẏ; x]`
"""
function linearize(dmp::AbstractDMP, y, ẏ, x, y0)
    error("g not taken care of")
    n = size(y,1)
    Ø = 0I
    ϕ = kernel_vector(dmp,x)
    ∇ϕ = - ϕ * (x-dmp.c)/dmp.σ2
    ∇ϕwx = ∇ϕ*dmp.w*x + ϕ*dmp.w
    A = -dmp.αz/dmp.τ * [Ø I Ø; dmp.βz/dmp.τ*I I Ø; Ø Ø dmp.αx*I] # Linear part
    A[n+1:2n, 2n+1:3n] = ∇ϕwx * (dmp.g-y0) / dmp.τ^2# Linearized force term
end


include("DMP2dof.jl")
end # module
