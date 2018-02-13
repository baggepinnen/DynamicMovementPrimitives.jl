module TwoLink

export torque, forward_kin, inverse_kin, inverse_kin_up, inverse_kin_down, traj, connect_points, acceleration, time_derivative, time_derivative!, inertia, ⊕
⊕(a::Tuple{Number,Number}, b::Tuple{Number,Number}) = (a[1]+b[1], a[2]+b[2])

const v1 = const v2 = 2
const m1 = const m2 = 0.2
const l1 = const l2 = 1
const k1 = const k2 = 0.1
const g = 1.82

function inertia(q2)
    x = cos(q2)*l1*l2*m2 + l2^2*m2
    [2*cos(q2)*l1*l2*m2+l2^2*m2+1    x;
    x                    l2^2*m2]
end
inertia(q1,q2) = inertia(q2)

signfunc(x) = tanh(100x) # A friendlier (smoother) version of sign(x)

"""
`τ1,τ2 = torque(q,qd,qdd)`\n
Inverse model
"""
function torque(q,qd,qdd)
    q1,q2     = q
    qd1,qd2   = qd
    qdd1,qdd2 = qdd

    c2  = cos(q2)
    s1  = sin(q1)
    s2  = sin(q2)
    s12 = sin(q1+q2)

    τ1 = m2*l2^2*(qdd1+qdd2) + m2*l1*l2*c2*(2qdd1+qdd2) +
    (m1+m2)*l1^2+qdd1 - m2*l1*l2*s2*qd2^2 - 2m2*l1*l2*s2*qd1*qd2 +
    g*l2*m2*s12 + (m1+m2)*l1*g*s1 + v1*qd1 +
    k1*signfunc(qd1)

    τ2 = m2*l1*l2*c2*qdd1 + m2*l1*l2*s2*qd1^2 + g*l2*m2*s12 +
    m2*l2^2*(qdd1 + qdd2) + v2*qd2 + k2*signfunc(qd2)

    (τ1,τ2)
end

function torque(q::AbstractMatrix,qd::AbstractMatrix,qdd::AbstractMatrix)
    tau = similar(q)
    for i = 1:size(q,2)
        tau1,tau2 = torque(q[:,i], qd[:,i],qdd[:,i])
        tau[1,i] = tau1
        tau[2,i] = tau2
    end
    tau
end

"""
`q̈ = acceleration(τ::Tuple, q::Tuple, qd::Tuple)`\n
`q̈1,q̈2 = acceleration(q1,q2,qd1,qd2,τ1,τ2)`\n
Model
"""
function acceleration(q1,q2,qd1,qd2,τ1,τ2)
    c2 = cos(q2)
    s1 = sin(q1)
    s2 = sin(q2)
    s12 = sin(q1+q2)

    qdd1 = (-l2*(g*l1*m1*s1 + g*l1*m2*s1 + g*l2*m2*s12 + k1*signfunc(qd1) + l1^2*m1 + l1^2*m2 - 2*l1*l2*m2*qd1*qd2*s2 - l1*l2*m2*qd2^2*s2 + qd1*v1 - τ1) + (c2*l1 + l2)*(g*l2*m2*s12 + k2*signfunc(qd2) + l1*l2*m2*qd1^2*s2 + qd2*v2 - τ2))/(l2*(2*c2*l1*l2*m2 + l2^2*m2 - m2*(c2*l1 + l2)^2 + 1))

    qdd2 = (l2*m2*(c2*l1 + l2)*(g*l1*m1*s1 + g*l1*m2*s1 + g*l2*m2*s12 + k1*signfunc(qd1) + l1^2*m1 + l1^2*m2 - 2*l1*l2*m2*qd1*qd2*s2 - l1*l2*m2*qd2^2*s2 + qd1*v1 - τ1) - (2*c2*l1*l2*m2 + l2^2*m2 + 1)*(g*l2*m2*s12 + k2*signfunc(qd2) + l1*l2*m2*qd1^2*s2 + qd2*v2 - τ2))/(l2^2*m2*(2*c2*l1*l2*m2 + l2^2*m2 - m2*(c2*l1 + l2)^2 + 1))

    qdd1,qdd2
end

function acceleration(τ::Tuple, q::Tuple, qd::Tuple)
    q1,q2   = q
    qd1,qd2 = qd
    τ1,τ2   = τ
    acceleration(q1,q2,qd1,qd2,τ1,τ2)
end

function acceleration(τ, q, q̇)
    q̈ = similar(q, dims=(2,size(q,2)))
    for i = 1:size(q,2)
        q̈1,q̈2 = acceleration(q[:,i]..., q̇[:,i]...,τ[:,i]...)
        q̈[1,i] = q̈1
        q̈[2,i] = q̈2
    end
    q̈
end

function time_derivative(state, τ)
    qdd1,qdd2 = acceleration(state[1],state[2],state[3],state[4],τ[1],τ[2])
    return [state[3],state[4], qdd1,qdd2]
end

function time_derivative!(state, τ, deriv)
    qdd1,qdd2 = acceleration(state[1],state[2],state[3],state[4],τ[1],τ[2])
    deriv[1] = state[3]
    deriv[2] = state[4]
    deriv[3] = qdd1
    deriv[4] = qdd2
end

# @syms v1 v2 k1 k2 m1 m2 l1 l2 g q1 q2 qd1 qd2 qdd1 qdd2 t1 t2 c1 c2 s1 s2 s12
# tau1 = m2*l2^2*(qdd1+qdd2) + m2*l1*l2*c2*(2qdd1+qdd2) +
# (m1+m2)*l1^2+qdd1 - m2*l1*l2*s2*qd2^2 - 2m2*l1*l2*s2*qd1*qd2 +
# g*l2*m2*s12 + (m1+m2)*l1*g*s1 + v1*qd1 +
# k1*sign(qd1)
# tau2 = m2*l1*l2*c2*qdd1 + m2*l1*l2*s2*qd1^2 + g*l2*m2*s12 +
# m2*l2^2*(qdd1 + qdd2) + v2*qd2 + k2*sign(qd2)
# sol = solve([t1-tau1,t2-tau2],[qdd1,qdd2])

function jacobian(q)
    q1,q2 = q[1],(q[1]+q[2])
    J = [l1*cos(q1)+l2*cos(q2) l2*cos(q2); l1*sin(q1)+l2*sin(q2) l2*sin(q2)]
end

function forward_kin(q)
    q1,q2 = q[1],(q[1]+q[2])
    p1 = (l1*sin(q1), -l1*cos(q1))
    p = p1 ⊕ (l2*sin(q2), -l2*cos(q2))
    p1, p
end

function forward_kin(jtraj::Matrix)
    N = size(jtraj,2)
    p = zeros(eltype(jtraj), 2, N)
    for i = 1:N
        q = jtraj[1:2,i]
        p[:,i] = forward_kin(q)[2] |> collect
    end
    p
end

"""
    (q11,q21), (q12,q22) = inverse_kin(p)
"""
function inverse_kin(p)
    x,y = p[1], p[2]
    a2 = x^2 + y^2
    a = sqrt(a2)
    α = acos((l1^2 + l2^2 - a2)/(2l1*l2))
    α2 = 2π-α
    q21 = π - α
    q22 = π - α2

    β = asin(sin(-α)/a*l2)
    β2 = -β
    q11 = asin(x/a)-β2
    q12 = asin(x/a)-β
    (q11,q21), (q12,q22)
end

function inverse_kin_up(p)
    q1,q2 = inverse_kin(p)
    return q1[2] > q2[2] ? q1 : q2
end

function inverse_kin_down(p)
    q1,q2 = inverse_kin(p)
    return q1[2] < q2[2] ? q1 : q2
end

"""
    p = inverse_kin(ctraj::Matrix, dir = :up)

`p` is of size N × 2
"""
function inverse_kin(ctraj::Matrix, dir = :up)
    if dir == :up
        kin_fun = inverse_kin_up
    elseif dir == :down
        kin_fun = inverse_kin_down
    else
        error("Symbol direction unknown")
    end
    N = size(ctraj,2)
    p = zeros(eltype(ctraj), 2, N)
    for i = 1:N
        q = ctraj[:,i]
        p[:,i] = kin_fun(q) |> collect
    end
    p
end

"""
    p, pd, pdd = traj(q0,q1,t)

Move from `q0` to `q1` during total time `t` using quadratic blends.
"""
function traj(q0,q1,t)
    tf = maximum(t)
    V = (q1-q0)/tf * 1.5
    traj(q0,q1,t, V)
end

"""
    p, pd, pdd = traj(q0,q1,t, V)

Move from `q0` to `q1` with velocity `V` and total time `t` using quadratic blends.
"""
function traj(q0,q1,t, V)
    tf = maximum(t)
    V  = abs(V) * sign(q1-q0)
    if abs(V) < abs(q1-q0)/tf
        error("V too small")
    elseif abs(V) > 2*abs(q1-q0)/tf
        error("V too big")
    end

    if q0 == q1
        s = ones(eltype(q0),size(t)) * q0
        sd = zeros(eltype(q0),size(t))
        sdd = zeros(eltype(q0),size(t))
        return s, sd, sdd
    end

    tb = (q0 - q1 + V*tf)/V
    a = V/tb

    p = Array{eltype(q0)}(size(t))
    pd = Array{eltype(q0)}(size(t))
    pdd = Array{eltype(q0)}(size(t))

    for (i,t) = enumerate(t)
        if t <= tb
            # initial blend
            p[i]   = q0 + a/2*t^2
            pd[i]  = a*t
            pdd[i] = a
        elseif t <= tf-tb
            # linear motion
            p[i]   = (q1+q0-V*tf)/2 + V*t
            pd[i]  = V
            pdd[i] = 0
        else
            # final blend
            p[i]   = q1 - a/2*tf^2 + a*tf*t - a/2*t^2
            pd[i]  = a*tf - a*t;
            pdd[i] = -a;
        end
    end
    return p, pd, pdd

end

"""
    p, pd, pdd = connect_points(points,ni)
Connect the points in matrix `points` ∈ R(ni+1 × 2) using third degree polinomials

# Example
```
traj = connect_points([0 1; 1 -0.5],ni)
q    = inverse_kin(traj[1],:up)'
qd   = centraldiff(q')'/sys.h
qdd  = centraldiff(qd')'/sys.h
u    = torque.(q,qd,qdd)
```
"""
function connect_points(points,ni)
    fx(i)    = traj(points[i,1],points[i+1,1],1:ni)
    fy(i)    = traj(points[i,2],points[i+1,2],1:ni)
    n        = size(points,2)
    x,xd,xdd = fx(1)
    y,yd,ydd = fy(1)
    p        = [x y]
    pd       = [xd yd]
    pdd      = [xdd ydd]
    for i = 2:n-1
        x,xd,xdd = fx(i)
        y,yd,ydd = fy(i)
        p        = vcat(p,[x y])
        pd       = vcat(pd,[xd yd])
        pdd      = vcat(pdd,[xdd ydd])
    end
    p', pd', pdd'
end


function plotworkspace(N=1000)
    q = [(2π*rand(),2π*rand()) for n in 1:N]
    p = getindex.(forward_kin.(q), 2)
    p1,p2 = getindex.(p,1), getindex.(p,2)
    scatter(p1,p2)
    p1,p2
end

end
