module DynamicMovementPrimitives
using Plots, ODE, Debug


export DMPopts, centraldiff,fit, solve, force, acceleration, solve_canonical


function centraldiff(v)
    c = size(v,2)
    dv = diff(v)/2
    a1 = [zeros(1,c);dv]
    a2 = [dv;zeros(1,c)]
    a = a1+a2
end

"""Takes a n vector of m vectors and creates a n×m matrix"""
vv2m(x::Vector) = [x[i][j] for i in eachindex(x), j in eachindex(x[1])]

"""
`DMPopts(Nbasis,αx,αz) = DMPopts(Nbasis,αx,αz,αz/4)`
Holds parameters for fitting a DMP
# Fields
`Nbasis,αx,αz,βz`
"""
type DMPopts
    Nbasis::Int
    αx::Real
    αz::Real
    βz::Real
end

DMPopts(Nbasis,αx,αz) = DMPopts(Nbasis,αx,αz,αz/4)

"""
The result of fitting a DMP
#Fields
`opts,g,y,ẏ,ÿ,w,τ,c,σ2`
"""
type DMP
    opts::DMPopts
    g::Vector{Float64}
    y::Matrix{Float64}
    ẏ::Matrix{Float64}
    ÿ::Matrix{Float64}
    w::Matrix{Float64}
    τ::Float64
    c::Vector{Float64}
    σ2::Vector{Float64}
end



function get_centers_linear(Nbasis)
    Ni = 1/(Nbasis+1)
    return Ni:Ni:1-Ni
end

function get_centers_log(Nbasis)
    Ni = 1/(Nbasis+1)
    return (logspace(log10(Ni),log10(1-Ni),Nbasis))[end:-1:1]
end

comp(x) = (x)

function kernel_matrix(x::AbstractVector,c,σ2)
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
_y0(y::VecOrMat) = y[1,:][:]
_y0(dmp::DMP) = _y0(dmp.y)
_T(dmp::DMP) = size(dmp.y,1)


"""
`fit(y, ẏ, ÿ, opts, g=y[end])`

Fits a DMP to data\n
`y, ẏ, ÿ` are position, velocity and acceleration respectively, `T×n` matrices where `T` is the number of time steps and `n` is the number of degrees of freedom.

"""
function fit(y,ẏ,ÿ,opts,g=y[end,:][:])

    T,n = size(y)
    τ = T/3 # After three time constants we have reached 1-exp(-3) ≈ 0.95
    Nbasis = opts.Nbasis
    βz = opts.βz
    αz = opts.αz
    αx = opts.αx
    y0 = _y0(y)
    x = solve_canonical(αx,τ,T)
    σ2 = (0.5/Nbasis)^2 * ones(Nbasis)

    ft = τ^2*ÿ - αz*(βz*(g'.-y)-τ*ẏ)
    ξ = x.*(g-y0)'
    c = get_centers_linear(Nbasis)
    Ψ = kernel_matrix(x,c,σ2)
    plot(Ψ);gui()

    w = zeros(Nbasis,n)
    for i = 1:n #joints
        for j = 1:Nbasis
            sTΓ = (ξ[:,i].*Ψ[:,j])'
            w[j,i] = (sTΓ*ft[:,i]/(sTΓ*ξ[:,i]))[1]
        end
    end
    return DMP(opts, g, y, ẏ, ÿ, w, τ,c,σ2)
end

"""
`force(dmp,x)`

Calculate the forcing term for `dmp` when the phase variable is `x`\n
`x` can be a scalar or a vector\n
The return value will be a `n` Vector or `T×n` Matrix depending on `typeof(x)`
"""
function force(d::DMP,x::AbstractVector)
    Ψ   = kernel_matrix(x,d.c,d.σ2)
    f   = Ψ*d.w .* (x.*(d.g-_y0(d))')
end

function force(d::DMP,x::Number)
    Ψ   = kernel_vector(x,d.c,d.σ2)
    f   = vec(Ψ'd.w) .* (x*(d.g-_y0(d)))
end

function acceleration(d::DMP,y,ẏ,x,g,n)
    αz = d.opts.αz
    βz = d.opts.βz
    f = force(d,x)[n]
    αz*(βz*(g'-y)-ẏ)+f
end

"""
    `solve(dmp::DMP, t = 0:_T(dmp)-1; y0 = _y0(dmp), g = dmp.g, solver=ode45)`

`t` time vector

## Keyword arguments: \n
`y0` start position, defaults to the initial point in training data from `dmp`
`g` goal, defaults to goal from `dmp`\n
`solver` the ode solver to use, see https://github.com/JuliaLang/ODE.jl \n
"""
function solve(dmp::DMP, t = 0:_T(dmp)-1; y0 = _y0(dmp), g = dmp.g, solver=ode45)
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



end # module
