module TwoLink


export torque, forward_kin, inverse_kin, inverse_kin_up, inverse_kin_down, traj, connect_points, acceleration, time_derivative, time_derivative!, inertia

import Base: +
TwoTuple = Tuple{Number,Number}
+(a::TwoTuple,b::TwoTuple) = (a[1]+b[1],a[2]+b[2])

const v1 = const v2 = 2
const m1 = const m2 = 0.2
const l1 = const l2 = 1
const k1 = const k2 = 0.00
const g = 1.82

function inertia(q2)
    x = cos(q2)*l1*l2*m2 + l2^2*m2
    [2*cos(q2)*l1*l2*m2+l2^2*m2+1    x;
              x                    l2^2*m2]
end
inertia(q1,q2) = inertia(q2)

signfunc(x) = tanh(0.01x) # A friendlier (smoother) version of sign(x)

"""
`[τ1,τ2] = torque(q,qd,qdd)`\n
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
    m2*l2*g*s12 + (m1+m2)*l1*g*s1 + v1*qd1 +
    k1*signfunc(qd1)

    τ2 = m2*l1*l2*c2*qdd1 + m2*l1*l2*s2*qd1^2 + m2*l2*g*s12 +
    m2*l2^2*(qdd1 + qdd2) + v2*qd2 + k2*signfunc(qd2)

    [τ1,τ2]
end

"""
`q̈1,q̈2 = acceleration(τ::VecOrTuple, q::VecOrTuple, qd::VecOrTuple)`\n
`q̈1,q̈2 = acceleration(q1,q2,qd1,qd2,τ1,τ2)`\n
Model
"""
function acceleration(q1,q2,qd1,qd2,τ1,τ2)
    c2 = cos(q2)
    s1 = sin(q1)
    s2 = sin(q2)
    s12 = sin(q1+q2)

    qdd1 = (-l2*(g*l1*m1*s1 + g*l1*m2*s1 + g*l2*m2*s12 + k1*signfunc(qd1) + l1^2*m1 + l1^2*m2 - 2*l1*l2*m2*qd1*qd2*s2 - l1*l2*m2*qd2^2*s2 + qd1*v1 - τ1) + (c2*l1 + l2)*(g*l2*m2*s12 + k2*signfunc(qd2) + l1*l2*m2*qd1^2*s2 + qd2*v2 - τ2))/(l2*(2*c2*l1*l2*m2 + l2^2*m2 - m2*(c2*l1 + l2)^2 + 1))

    qdd2 = (l2*m2*(c2*l1 + l2)*(g*l1*m1*s1 + g*l1*m2*s1 + g*l2*m2*s12 + k1*signfunc(qd1) + l1^2*m1 + l1^2*m2 - 2*l1*l2*m2*qd1*qd2*s2 - l1*l2*m2*qd2^2*s2 + qd1*v1 - τ1) - (2*c2*l1*l2*m2 + l2^2*m2 + 1)*(g*l2*m2*s12 + k2*sign(qd2) + l1*l2*m2*qd1^2*s2 + qd2*v2 - τ2))/(l2^2*m2*(2*c2*l1*l2*m2 + l2^2*m2 - m2*(c2*l1 + l2)^2 + 1))

    qdd1,qdd2
end

function acceleration(τ, q, qd)
    q1,q2   = q
    qd1,qd2 = qd
    τ1,τ2   = τ
    acceleration(q1,q2,qd1,qd2,τ1,τ2)
end

function time_derivative(state, τ)
    qdd1,qdd2 = acceleration(state[1],state[2],state[3],state[4],τ[1],τ[2])
    return [state[3],state[4], qdd1,qdd2]
end

function time_derivative!(state, τ, deriv)
    qdd1,qdd2 = acceleration(state[1],state[2],state[3],state[4],τ[1],τ[2])
    deriv[:] = [state[3],state[4], qdd1,qdd2]
end

# @syms v1 v2 k1 k2 m1 m2 l1 l2 g q1 q2 qd1 qd2 qdd1 qdd2 t1 t2 c1 c2 s1 s2 s12
# tau1 = m2*l2^2*(qdd1+qdd2) + m2*l1*l2*c2*(2qdd1+qdd2) +
# (m1+m2)*l1^2+qdd1 - m2*l1*l2*s2*qd2^2 - 2m2*l1*l2*s2*qd1*qd2 +
# m2*l2*g*s12 + (m1+m2)*l1*g*s1 + v1*qd1 +
# k1*sign(qd1)
# tau2 = m2*l1*l2*c2*qdd1 + m2*l1*l2*s2*qd1^2 + m2*l2*g*s12 +
# m2*l2^2*(qdd1 + qdd2) + v2*qd2 + k2*sign(qd2)
# sol = solve([t1-tau1,t2-tau2],[qdd1,qdd2])

function forward_kin(q)
    q1,q2 = q[1],(q[1]+q[2])
    p1 = (l1*sin(q1), -l1*cos(q1))
    p = p1 + (l2*sin(q2), -l2*cos(q2))
    p1, p
end

function forward_kin(jtraj::Matrix)
    p = zeros(jtraj)
    for i = 1:size(jtraj,1)
        q = jtraj[i,:][:]
        p[i,:] = forward_kin(q)[2] |> collect
    end
    p
end


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

function inverse_kin(ctraj::Matrix, dir = :up)
    if dir == :up
        kin_fun = inverse_kin_up
    elseif dir == :down
        kin_fun = inverse_kin_down
    else
        error("Symbol direction unknown")
    end

    p = zeros(ctraj)
    for i = 1:size(ctraj,1)
        q = ctraj[i,:][:]
        p[i,:] = kin_fun(q) |> collect
    end
    p
end

function traj(q0,q1,t)
    tf = maximum(t)
    V = (q1-q0)/tf * 1.5
    traj(q0,q1,t, V)
end

function traj(q0,q1,t, V)
    tf = maximum(t)
    V  = abs(V) * sign(q1-q0)
    if abs(V) < abs(q1-q0)/tf
        error("V too small")
    elseif abs(V) > 2*abs(q1-q0)/tf
        error("V too big")
    end

    if q0 == q1
        s = ones(Float64,size(t)) * q0
        sd = zeros(Float64,size(t))
        sdd = zeros(Float64,size(t))
        return s, sd, sdd
    end

    tb = (q0 - q1 + V*tf)/V
    a = V/tb

    p = Array(Float64,size(t))
    pd = Array(Float64,size(t))
    pdd = Array(Float64,size(t))

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

function connect_points(points,ni)
    fx(i)    = traj(points[i,1],points[i+1,1],1:ni)
    fy(i)    = traj(points[i,2],points[i+1,2],1:ni)
    n        = size(points,1)
    x,xd,xdd = fx(1)
    y,yd,ydd = fy(1)
    p        = [x y]
    pd       = [xd yd]
    pdd      = [xdd ydd]
    for i = 2:n-1
        x,xd,xdd = fx(i)
        y,yd,ydd = fy(i)
        p        = cat(1,p,[x y])
        pd       = cat(1,pd,[xd yd])
        pdd      = cat(1,pdd,[xdd ydd])
    end
    p, pd, pdd
end

end
