using Distributed
@everywhere using LinearAlgebra
@everywhere using Printf

function rk4_np(fun, x_k, t_k, dt, par)
    f1 = fun(x_k, t_k, par)
    f2 = fun(x_k .+ (dt / 2.0) .* f1, t_k + dt / 2.0, par)
    f3 = fun(x_k .+ (dt / 2.0) .* f2, t_k + dt / 2.0, par)
    f4 = fun(x_k .+ dt .* f3, t_k + dt, par)
    #x_k1 .= x_k .+ (dt / 6.0) .* (f1 .+ 2.0 .* f2 .+ 2.0 .* f3 .+ f4)

    return  x_k .+ (dt / 6.0) .* (f1 .+ 2.0 .* f2 .+ 2.0 .* f3 .+ f4)

end


function derivs_duff_simple(y, t, par)
    return [ y[2]
            -0.1*y[2]+y[1]-y[1]^3]
end


@everywhere function integration( par, init_cond, t0, tf, dt)



    function derivs_duff_simple(y, t, par)
        return [ y[2]
                -0.1*y[2]+y[1]-y[1]^3]
    end

    function derivs_Li_MSSP_2023(y, t, par)
        omega, m1, m2, k1, k2, k3, l1, l2, l3, h, c1, c2, induc, alpha, res, amp = par
        dydt = [
            y[2]
            -k1/m1 * y[1] * (1 - l1 / sqrt(y[1]^2 + (h - y[3])^2)) - c1/m1 * y[2] - (-amp * cos(omega * t)) - alpha/m1 * y[6]
            y[4]
            -k2/m2 * y[3] - 2 * k3/m2 * y[3] * (1 - l3 / sqrt(l3^2 + y[3]^2)) + k1/m2 * (h - y[3]) * (1 - l1 / sqrt(y[1]^2 + (h - y[3])^2)) - c2/m2 * y[4]
            y[6]
            -res/induc * y[6] + alpha/induc * y[2]
        ]
        return dydt
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
        eval(Meta.parse(raw_eq))
        return dydt
    end

    raw_eq = "p1, = par\ndydt = [y[1], -p1*y[1]+y[0]-pow(y[0],3)]"  
    time = t0:dt:tf
    num_steps = length(time)
    y_arr = zeros(num_steps, length(init_cond))
    y = copy(init_cond)
    y_out = similar(y)
    #poincare_matrix = Matrix{Float64}(undef, 0, sys_dim+1)

    for i in 1:num_steps
        y_arr[i, :] .= y
        rk4_np!(derivs_Li_MSSP_2023_nondim, y, time[i], dt, par, y_out)
        y .= y_out
    end
    return hcat(time, y_arr)

end

