using Distributed
@everywhere using LinearAlgebra
@everywhere using Printf

@everywhere function integration( par, init_cond, t0, tf, dt)

    function derivs_rossler(y, t, par)
        a,b,c = par
        return [-y[2]-y[3],
                y[1]+a*y[2],
                b+y[3]*(y[1]-c)]
    end

    function derivs_duff_simple(y, t, par)
        return [ y[2]
                -0.1*y[2]+y[1]-y[1]^3]
    end

    function derivs_Li_MSSP_2023_nondim(y, t, par)
        omega, P, lambd, ni1, ni2, phi1, phi2, mi1, mi2, theta, epsilon, ro = par
        
        dydt = [
            y[2],
            -y[1] * (1 - 1 / sqrt(y[1]^2 + ni1^2 * (1 - y[3])^2)) - phi1 * y[2] + P * cos(omega * t) - ro * y[5],
            y[4],
            -lambd / mi1 * y[3] - 2 * lambd / mi2 * y[3] * (1 - 1 / sqrt(1 + ni2^2 * y[3]^2)) + lambd * (1 - y[3]) * (1 - 1 / sqrt(y[1]^2 + ni1^2 * (1 - y[3])^2)) - phi2 * y[4],
            -theta * y[5] + epsilon * y[2]
        ]
        
        return dydt
    end
    
    function rk4_np!(fun, x_k, t_k, dt, par, x_k1)
        f1 = fun(x_k, t_k, par)
        f2 = fun(x_k .+ (dt / 2.0) .* f1, t_k + dt / 2.0, par)
        f3 = fun(x_k .+ (dt / 2.0) .* f2, t_k + dt / 2.0, par)
        f4 = fun(x_k .+ dt .* f3, t_k + dt, par)
        x_k1 .= x_k .+ (dt / 6.0) .* (f1 .+ 2.0 .* f2 .+ 2.0 .* f3 .+ f4)
    
        #return  x_k .+ (dt / 6.0) .* (f1 .+ 2.0 .* f2 .+ 2.0 .* f3 .+ f4)
    
    end

    function derivs(y, t, par)
        println(y)
        eval(Meta.parseall(raw_eq))
        return dydt
    end

    #raw_eq = """dydt = [y[2], -0.3*y[2]+y[1]-y[1]^3]"""
    #raw_derivs = @eval function ()
    #    $raw_eq
    #end
    time = t0:dt:tf
    num_steps = length(time)
    y_arr = zeros(num_steps, length(init_cond))
    y = copy(init_cond)
    y_out = similar(y)

    for i in 1:num_steps
        y_arr[i, :] .= y
        rk4_np!(derivs_duff_simple, y, time[i], dt, par, y_out)
        y .= y_out
    end
    return hcat(time, y_arr)

end

