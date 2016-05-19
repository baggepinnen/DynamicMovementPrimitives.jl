module DynamicMovementPrimitives
using ODE, Requires


export DMPopts, centraldiff,fit, solve, force, acceleration, solve_canonical, plotdmp


function centraldiff(v::AbstractMatrix)
    dv = diff(v)/2
    a1 = [dv[1,:];dv]
    a2 = [dv;dv[end,:]]
    a = a1+a2
end

function centraldiff(v::AbstractVector)
    dv = diff(v)/2
    a1 = [dv[1];dv]
    a2 = [dv;dv[end]]
    a = a1+a2
end

"""Takes an n vector of m vectors and creates a n×m matrix"""
vv2m(x::Vector) = [x[i][j] for i in eachindex(x), j in eachindex(x[1])]

"""
`DMPopts(Nbasis,αx,αz) = DMPopts(Nbasis,αx,αz,αz/4)`
Holds parameters for fitting a DMP
# Fields
`Nbasis,αx,αz,βz,sched_sig`\n
`sched_sig` can be chosen as `:canonical` (default), `:time` or `position`

See example file or the paper by Ijspeert et al. 2013
"""
immutable DMPopts
    Nbasis::Int
    αx::Real
    αz::Real
    βz::Real
    sched_sig::Symbol
end

DMPopts(Nbasis,αx,αz) = DMPopts(Nbasis,αx,αz,αz/4,:position)
DMPopts(Nbasis,αx,αz,sched_sig::Symbol) = DMPopts(Nbasis,αx,αz,αz/4,sched_sig)

"""
The result of fitting a DMP
#Fields
`opts,g,y,ẏ,ÿ,w,τ,c,σ2`\n
See example file or the paper by Ijspeert et al. 2013
"""
immutable DMP
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
end

function get_centers_linear(Nbasis,x)
    ma = maximum(x,1)
    mi = minimum(x,1)
    n = size(x,2)
    d = ma-mi
    Ni = d./(Nbasis+1)


    σ2 = zeros(Nbasis,n)
    c = zeros(Nbasis,n)
    for i = 1:n
        σ2[:,i] = d[i]*(0.5/Nbasis)^2 * ones(Nbasis)
        c[:,i]  = linspace(mi[i]+Ni[i],ma[i]-Ni[i],Nbasis)
    end
    return c, σ2
end

# function get_centers_log(Nbasis)
#     Ni = 1/(Nbasis+1)
#     return (logspace(log10(Ni),log10(1-Ni),Nbasis))[end:-1:1]
# end

comp(x) = (x)

function kernel_matrix(x::AbstractVecOrMat,c,σ2)
    Ψ = Float64[exp(-1/(2σ2[j])*(comp(x)-(c[j]))^2) for x in x, j in eachindex(c)]
    Ψ ./= sum(Ψ,2)
    Ψ
end

function kernel_vector(x::Real,c,σ2)
    Ψ = Float64[exp(-1/(2σ2[j])*(comp(x)-(c[j]))^2) for j in eachindex(c)]
    Ψ ./= sum(Ψ)
    Ψ
end

solve_canonical(αx,τ,T::AbstractVector) = exp(-αx/τ.*T)
solve_canonical(αx,τ,T::Real) = solve_canonical(αx,τ,(0:T-1))
solve_canonical(dmp::DMP,t) = solve_canonical(dmp.opts.αx, dmp.τ, t)
_1(y::VecOrMat) = y[1,:][:]
_1(dmp::DMP) = _1(dmp.y)
_T(dmp::DMP) = size(dmp.y,1)

function get_sched_sig(s,αx,τ,t,y)

    if s == :canonical
        return solve_canonical(αx,τ,t)
    elseif s == :position
        return y
    elseif s == :time
        return t
    end
    warn("Scheduling signal $s unknown")
    return solve_canonical(αx,τ,t)
end
get_sched_sig(dmp::DMP) = get_sched_sig(dmp.opts.sched_sig,dmp.opts.αx,dmp.τ,dmp.t)


"""
`fit(y, ẏ, ÿ, opts, g=y[end])`

Fits a DMP to data\n
`y, ẏ, ÿ` are position, velocity and acceleration respectively, `T×n` matrices where `T` is the number of time steps and `n` is the number of degrees of freedom.

see also `solve`, `plotdmp`
"""
function fit(y,ẏ,ÿ,t,opts,g=y[end,:][:])
    T       = t[end]
    N       = size(y,1)
    n       = isa(y,Matrix) ? size(y,2) : 1
    τ       = T/3 # After three time constants we have reached 1-exp(-3) ≈ 0.95
    Nbasis  = opts.Nbasis
    βz      = opts.βz
    αz      = opts.αz
    αx      = opts.αx
    y0      = _1(y)
    x       = get_sched_sig(opts.sched_sig,αx,τ,t,y)
    ft      = τ^2*ÿ - αz*(βz*(g'.-y)-τ*ẏ)
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
        for j = 1:Nbasis
            sTΓ = ξ[:,i].*Ψ[:,j]
            w[j,i] = vecdot(sTΓ,ft[:,i])/vecdot(sTΓ,ξ[:,i])
        end
    end
    return DMP(opts, g, y, ẏ, ÿ,t, w, τ,c,σ2)
end

"""
`force(dmp,x)`

Calculate the forcing term for `dmp` when the phase variable is `x`\n
The return value will be an `n` Vector or `T×n` Matrix depending on `typeof(x)`
"""
function force(d::DMP, x)
    if d.opts.sched_sig == :position
        force_multiple(d,x)
    else
        force_single(d,x)
    end
end

function force_single(d::DMP,x::Number, i)
    # ODE Point case
    y0 = _1(d)
    Ψ    = kernel_vector(x, d.c, d.σ2)
    f = vecdot(Ψ,d.w[:,i]) * x*(d.g[i]-y0[i])
end

function force_single(d::DMP,x::Number)
    # Point case
    y0 = _1(d)
    Ψ  = kernel_matrix([x], d.c, d.σ2)
    f  = (Ψ*d.w)[:] .* (x*(d.g-y0))
end

function force_single(d::DMP,x::AbstractVector)
    # Trajectory case
    y0 = _1(d)
    Ψ    = kernel_matrix(x, d.c, d.σ2)
    f = Ψ*d.w .* (x.*(d.g-y0)')
end

function force_multiple(d::DMP,x::Number, i)
    # ODE Point case
    y0 = _1(d)
    Ψ  = kernel_vector(x[i], d.c[:,i], d.σ2[:,i])
    f  = vecdot(Ψ,d.w[:,i]) * x*(d.g[i]-y0[i])
end

function force_multiple(d::DMP,x::AbstractVector)
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

function force_multiple(d::DMP,x::AbstractMatrix)
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

function force(d::DMP,x::Number, i)
    p = size(d.c,2) # Number of sched signals, can be either 1 or n
    i = p > 1 ? i : 1
    y0 = _1(d)
    Ψ  = kernel_vector(x, d.c[:,i], d.σ2[:,i])
    f = vecdot(Ψ,d.w[:,i]) * x*(d.g[i]-y0[i])
end


function acceleration(d::DMP, y::Number,ẏ::Number,x::Number,g::Number,i=1)
    f = force(d,x,i)
    d.opts.αz*(d.opts.βz*(g-y)-ẏ)+f
end

function acceleration(d::DMP, y::AbstractVector,ẏ::AbstractVector,x ,g::AbstractVector)
    f = force(d,x)
    d.opts.αz*(d.opts.βz*(g-y)-ẏ)+f
end

function acceleration(d::DMP, y::AbstractMatrix,ẏ::AbstractMatrix,x::AbstractVecOrMat,g::AbstractVector)
    f = force(d,x)
    d.opts.αz*(d.opts.βz*(g'.-y)-ẏ)+f
end

"""
`solve(dmp::DMP, t = 0:_T(dmp)-1; y0 = _1(dmp), g = dmp.g, solver=ode45)`

`t` time vector

## Keyword arguments: \n
`y0` start position, defaults to the initial point in training data from `dmp`
`g` goal, defaults to goal from `dmp`\n
`solver` the ode solver to use, see https://github.com/JuliaLang/ODE.jl \n

see also `plotdmp`
"""
function solve(dmp::DMP, t = dmp.t; y0 = _1(dmp), g = dmp.g, solver=ode45)
    if dmp.opts.sched_sig == :position
        return solve_position(dmp, t, y0, g, solver)
    elseif dmp.opts.sched_sig == :time
        return solve_time(dmp, t, y0, g, solver)
    end
    solve_canonical(dmp, t, y0, g, solver)
end

function solve_canonical(dmp, t, y0, g, solver)
    T,n     = size(dmp.y)
    αx      = dmp.opts.αx
    τ       = dmp.τ
    z       = zeros(T,n)
    y       = zeros(T,n)
    x       = zeros(T)
    for i = 1:n
        function time_derivative(t,state)
            local z   = state[1]
            local y   = state[2]
            local x   = state[3]
            zp  = acceleration(dmp, y, z, x,g[i],i)
            yp  = z
            xp  = -αx * x
            [zp;yp;xp] / τ
        end
        state0  = [0; y0[i]; 1.]
        tout,state_history = solver(time_derivative, state0, t,points=:specified)
        res = vv2m(state_history)
        z[:,i] = res[:,1]/τ
        y[:,i] = res[:,2]
        x = res[:,3] # TODO: se till att denna är samma för alla DOF
    end
    t,y,z,x
end

function solve_position(dmp, t, y0, g, solver)
    T,n     = size(dmp.y)
    αx      = dmp.opts.αx
    τ       = dmp.τ
    z       = zeros(T,n)
    y       = zeros(T,n)
    for i = 1:n
        function time_derivative(t,state)
            local z   = state[1]
            local y   = state[2]
            zp  = acceleration(dmp, y, z, y,g[i],i)
            yp  = z
            [zp;yp] / τ
        end
        state0  = [0; y0[i]]
        tout,state_history = solver(time_derivative, state0, t,points=:specified)
        res = vv2m(state_history)
        z[:,i] = res[:,1]/τ
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
            [zp;yp] / τ
        end
        state0  = [0; y0[i]]
        tout,state_history = solver(time_derivative, state0, t,points=:specified)
        res = vv2m(state_history)
        z[:,i] = res[:,1]/τ
        y[:,i] = res[:,2]
    end
    t,y,z,t
end


plotdmp(dmp::DMP, args...) = println("To plot a DMP, install package Plots.jl or call tout,yout,ẏout,xout = solve(dmp) to produce your own plot.")


@require Plots begin
math(sl) = map(s->string("\$",s,"\$") ,sl)
function plotdmp(dmp::DMP; kwargs...)
    tout,yout,ẏout,xout = solve(dmp)
    n = size(dmp.y,2)
    fig = Plots.subplot(n = n, nc = 1)
    for i = 1:n
        Plots.plot!(fig[i,1],tout,[yout[:,i] ẏout[:,i]],lab = ["y_{out}" "ẏ_{out}"] |> math; kwargs...)
        Plots.plot!(fig[i,1],tout,[dmp.y[:,i] dmp.ẏ[:,i]],l=:dash,lab = ["y" "ẏ"] |> math; kwargs...)
    end
end
end

end # module
