module DynamicMovementPrimitives
using Plots, ODE, Debug


export DMPopts, centraldiff,fit, solve, force, acceleration




function centraldiff(v)
    c = size(v,2)
    dv = diff(v)/2
    a1 = [zeros(1,c);dv]
    a2 = [dv;zeros(1,c)]
    a = a1+a2
end

type DMPopts
    Nbasis
    αx
    αz
    βz
end

type DMP
    opts::DMPopts
    g
    y
    yd
    ydd
    w
    τ
    c
    σ2
end

DMPopts(Nbasis,αx,αz) = DMPopts(Nbasis,αx,αz,αz/4)

function get_centers_linear(Nbasis)
    Ni = 1/(Nbasis+1)
    return Ni:Ni:1-Ni
end

# function get_centers_log(Nbasis)
#     Ni = 1/(Nbasis+1)
#     return (1-logspace(log10(Ni),log10(1-Ni),Nbasis))[end:-1:1]
# end

function kernel_matrix(x::AbstractVector,c,σ2)
    Ψ = Float64[exp(-1/(2σ2[j])*(x[t]-c[j])^2) for t in eachindex(x), j in eachindex(c)]
    Ψ ./= sum(Ψ,2)
    Ψ
end

function kernel_vector(x::Real,c,σ2)
    Ψ = Float64[exp(-1/(2σ2[j])*(x-c[j])^2) for j in eachindex(c)]
    Ψ ./= sum(Ψ)
    Ψ
end

solve_canonical(αx,τ,T::AbstractVector) = exp(-αx/τ.*T)
solve_canonical(αx,τ,T::Real) = solve_canonical(αx,τ,(0:T-1))
_y0(y::VecOrMat) = y[1,:][:]
_y0(dmp::DMP) = _y0(dmp.y)


function fit(y,yd,ydd,opts,g=y[end,:][:])

    T,n = size(y)
    τ = T/3 # After three time constants we have reached 1-exp(-3) ≈ 0.95
    Nbasis = opts.Nbasis
    βz = opts.βz
    αz = opts.αz
    αx = opts.αx
    y0 = _y0(y)
    x = solve_canonical(αx,τ,T)
    σ2 = (1/Nbasis)^2 * ones(Nbasis)

    ft = τ^2*ydd - αz*(βz*(g.-y)-τ*yd)
    ξ = x.*(g-y0)'
    c = get_centers_linear(Nbasis)
    Ψ = kernel_matrix(x,c,σ2)


    w = zeros(Nbasis,n)
    for i = 1:n
        for j = 1:Nbasis
            sTΓ = (ξ[:,i].*Ψ[:,j])'
            w[j,i] = (sTΓ*ft[:,i]/(sTΓ*ξ[:,i]))[1]
        end
    end
    return DMP(opts, g, y, yd, ydd, w, τ,c,σ2)
end

function force(dmp::DMP,x::AbstractVector)
    Ψ   = kernel_matrix(x,dmp.c,dmp.σ2)
    f   = Ψ*dmp.w .* (x.*(dmp.g-_y0(dmp))')
end

function force(dmp::DMP,x::Number)
    Ψ   = kernel_vector(x,dmp.c,dmp.σ2)
    f   = vec(Ψ'dmp.w) .* (x*(dmp.g-_y0(dmp)))
end

function acceleration(dmp,y,yd,x,g)
    αz = dmp.opts.αz
    βz = dmp.opts.βz
    f = force(dmp,x)
    (αz*(βz*(g-y)-yd)+f)
end

function solve(dmp::DMP, t = 0:size(dmp.y,1)-1; y0 = _y0(dmp), g = dmp.g, solver=ode45)
    n = size(dmp.y,2)

    αx = dmp.opts.αx
    τ = dmp.τ
    # x = solve_canonical(αx,τ,tgrid)
    state0 = [zeros(y0); y0; 1.]
    function time_derivative(t,state)
        z = state[1:n]
        y = state[n+1:2n]
        x = state[end]
        zp = acceleration(dmp, y, z, x,g)
        yp = z
        xp = -αx * x
        [zp;yp;xp] ./ τ
    end
    # @bp
    tout,state_history = solver(time_derivative, state0, t)
    res = cat(2,state_history...)'

    z = res[:,1:n]
    y = res[:,n+1:2n]
    x = res[:,end]
    tout,y,z,x

end



end # module
