
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


_1(y::VecOrMat) = y[1,:][:]
_1(dmp::DMP) = _1(dmp.y)
_T(dmp::DMP) = size(dmp.y,1)


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

function plotdmpphase(dmp::DMP; kwargs...)
    tout,yout,ẏout,xout = solve(dmp)
    Plots.plot(yout[:,1],yout[:,2],lab = ["y_{out}" "ẏ_{out}"] |> math; kwargs...)
    Plots.plot!(dmp.y[:,1],dmp.y[:,2],l=:dash,lab = ["y" "ẏ"] |> math; kwargs...)
end
end
