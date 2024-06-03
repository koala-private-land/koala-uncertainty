function declining_discount_annuity(t::Integer = 1, schedule::String = "green_book")
    # t: Annuity starting at year t
    # schedule: declining discount rate schedule (either "green_book" following UK Green Book guidance, or "french", for 4% for the first 30 years, then 2% afterwards)

    if cmp(schedule, "green_book")==0
        T = 1:300
        if (t > T[end])
            throw("t > 300 not supported yet in the Green Book schedule discount rate")
        end
        rho = ones(301,1)
        r = vcat(ones(30,1) * 0.035, ones(45,1) * 0.03, ones(50, 1) * 0.025, ones(75, 1) * 0.02, ones(100, 1) * 0.015)
        for tt in T
            rho[tt+1] = rho[tt] * (1/(1+r[tt]))
        end
        rho_inf = rho[T[end]] / 0.01 # Annuity value to infinity
        rho_npv = sum(rho[t:end]) + rho_inf
    elseif cmp(schedule, "french")==0
        # Note: t > 30 is not currently supported
        T = 1:30
        if (t > T[end])
            throw("t > 30 not supported yet in the French schedule discount rate")
        end
        rho = ones(31,1)
        r = vcat(ones(30,1) * 0.04)
        for tt in T
            rho[tt+1] = rho[tt] * (1/(1+r[tt]))
        end
        rho_inf = rho[T[end]] / 0.02 # Annuity value to infinity
        rho_npv = sum(rho[t:end]) + rho_inf
    else
        throw("schedule argument invalid")
    end
    return(rho_npv)
end