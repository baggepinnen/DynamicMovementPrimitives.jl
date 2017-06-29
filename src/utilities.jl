
"""
`c, σ2 = get_centers_linear(Nbasis,x)`\n
`x` is the (multidimensional) scheduling signal
"""
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
        c[:,i]  = linspace(mi[i]+0Ni[i],ma[i]-0Ni[i],Nbasis)
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

function kernel_vector(dmp::AbstractDMP, x::Real)
    kernel_vector(x, dmp.c, dmp.σ2)
end


function centraldiff(v::AbstractMatrix)
    dv = diff(v)/2
    a1 = [dv[[1],:];dv]
    a2 = [dv;dv[[size(dv,1)],:]]
    a = a1+a2
end

function centraldiff(v::AbstractVector)
    dv = diff(v)/2
    a1 = [dv[1];dv]
    a2 = [dv;dv[end]]
    a = a1+a2
end

_1(y::VecOrMat) = y[1,:][:]
_1(dmp::AbstractDMP) = _1(dmp.y)
_T(dmp::AbstractDMP) = size(dmp.y,1)

math(sl) = map(s->string("\$",s,"\$") ,sl)

function euler(time_derivative, state0, t, args...; kwargs...)
    T = length(t)
    n = length(state0)
    res = Matrix{Float64}(T,n)
    res[1,:] = state0
    for i in 2:T
        res[i,:] = res[i-1,:] + time_derivative(t[i-1],res[i-1,:])*(t[2]-t[1])
    end
    t,res
end

"""
plot(dmp::DMP, phase::Bool=false; [y0])
"""
@recipe function plotdmp(dmp::AbstractDMP, phase::Bool=false; y0 = _1(dmp))
    n = size(dmp.y,2)
    if n < 2 && phase
        warn("Phase plotting only supported for DMP of dimension greater than 1")
        phase = false
    end
    layout := (n,1)
    tout,yout,ẏout,xout = solve(dmp, y0=y0)[1:4]
    delete!(d,:y0)
    if phase
        xguide := "Position 1"
        yguide := "Position 2"
    else
        xguide := "Time"
        yguide := "Output"
    end
    for i = 1:n
        @series begin
            title       := "Dynamic Movement Primitive $i"
            subplot     := i
            seriestype  := :path
            label       := ["y_{out}" "ẏ_{out}"] |> math
            if phase
                xguide := "Position 1"
                yguide := "Position 2"
                (yout[:,1], yout[:,2])
            else
                xguide := "Time"
                yguide := "Output"
                (tout, [yout[:,i] ẏout[:,i]])
            end

        end
        @series begin
            title       := "Dynamic Movement Primitive $i"
            subplot     := i
            seriestype  := :path
            linestyle   := :dash
            label       := ["y" "ẏ"] |> math
            phase ? (dmp.y[:,1], dmp.y[:,2]) : (tout, [dmp.y[:,i] dmp.ẏ[:,i]])
        end
    end

end





plotdmp(dmp::AbstractDMP, args...) = println("To plot a DMP, install package Plots.jl or call tout,yout,ẏout,xout = solve(dmp) to produce your own plot.")

@require Plots begin
math(sl) = map(s->string("\$",s,"\$") ,sl)
function plotdmp(dmp::AbstractDMP; kwargs...)
    tout,yout,ẏout,xout = solve(dmp)[1:4]
    n = size(dmp.y,2)
    fig = Plots.plot(layout = (n,1))
    for i = 1:n
        Plots.plot!(tout,[yout[:,i] ẏout[:,i]],lab = ["y_{out}" "ẏ_{out}"] |> math; kwargs...)
        Plots.plot!(tout,[dmp.y[:,i] dmp.ẏ[:,i]],l=:dash,lab = ["y" "ẏ"] |> math; kwargs...)
    end
    Plots.gui()
end

function plotdmpphase(dmp::AbstractDMP; kwargs...)
    tout,yout,ẏout,xout = solve(dmp)[1:4]
    Plots.plot(yout[:,1],yout[:,2],lab = ["y_{out}" "ẏ_{out}"] |> math; kwargs...)
    Plots.plot!(dmp.y[:,1],dmp.y[:,2],l=:dash,lab = ["y" "ẏ"] |> math; kwargs...)
    Plots.gui()
end

end
