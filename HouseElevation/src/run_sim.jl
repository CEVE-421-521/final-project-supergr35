using Distributions

"""Helper function for trapezoidal rule"""
function trapz(x, y)
    return sum((x[2:end] - x[1:(end - 1)]) .* (y[2:end] + y[1:(end - 1)])) * 0.5
end

"""
Run the model for a given action and SOW

Expected Annual Damages are computed using the trapezoidal rule
"""
function run_sim(a::Action, sow::SOW, p::ModelParams)

    # first, we calculate the cost of elevating the house
    construction_cost = elevation_cost(p.house, a.Δh_ft)

    # we don't need to recalculate the steps of the trapezoidal integral for each year
    storm_surges_ft = range(
        quantile(sow.surge_dist, 0.0005); stop=quantile(sow.surge_dist, 0.9995), length=130
    )

    eads = map(p.years) do year

        # get the sea level for this year
        slr_ft = sow.slr(year)

        # Compute EAD using trapezoidal rule
        pdf_values = pdf.(sow.surge_dist, storm_surges_ft) # probability of each
        depth_ft_gauge = storm_surges_ft .+ slr_ft # flood at gauge
        depth_ft_house = depth_ft_gauge .- (p.house.height_above_gauge_ft + a.Δh_ft) # flood @ house
        damages_frac = p.house.ddf.(depth_ft_house) ./ 100 # damage
        weighted_damages = damages_frac .* pdf_values # weighted damage
        # Trapezoidal integration of weighted damages
        ead = trapz(storm_surges_ft, weighted_damages) * p.house.value_usd
    end

    years_idx = p.years .- minimum(p.years)
    discount_fracs = (1 - sow.discount_rate) .^ years_idx
    ead_npv = sum(eads .* discount_fracs)
    return -(ead_npv + construction_cost)
end

"""
Run the sequential decision model for a given SeqAction and SOW
"""
function run_sim_seq(a::SeqAction, sow::SOW, p::ModelParams)
    # get an expected storm surge value
    storm_surges_ft = range(
        quantile(sow.surge_dist, 0.0005); stop=quantile(sow.surge_dist, 0.9995), length=130
    )
    pdf_values = pdf.(sow.surge_dist, storm_surges_ft) # probability of each
    exp_surge = sum(storm_surges_ft .* pdf_values)/sum(pdf_values) # weighted average to use for heightening heuristic
    house_height = p.house.height_above_gauge_ft

    eads = map(p.years) do year 
        # first, we need to determine the actual elevation height based on our policy
        slr_ft = sow.slr(year)
        exp_depth = slr_ft + exp_surge
        # only heighten if we expect to be below the buffer height
        if a.buff < house_height - exp_depth
            Δh = 0
        else
            Δh = exp_depth - (house_height - a.buff) + a.free
            # note that we can't elevate more than 14 feet
            if Δh > 14.0
                Δh = 14.0
            end
        end
        construction_cost = elevation_cost(p.house, Δh)
        depth_ft_gauge = storm_surges_ft .+ slr_ft # flood at gauge including uncertainty
        # now calculate damages for our heightening
        depth_ft_house = depth_ft_gauge .- (house_height += Δh) # flood @ house, updating height
        damages_frac = p.house.ddf.(depth_ft_house) ./ 100 # damage 
        weighted_damages = damages_frac .* pdf_values # weighted damage
        # Trapezoidal integration of weighted damages
        ead = trapz(storm_surges_ft, weighted_damages) * p.house.value_usd
        # each element of eads will now be a tuple of the form (damage, construction cost)
        # println(Δh)
        # println(house_height)
        tup = tuple(ead, construction_cost)
    end
    # we need to unzip the array of tuples into an array of EAD and construction costs
    eads, costs = (first.(eads), last.(eads))
    years_idx = p.years .- minimum(p.years)
    discount_fracs = (1 - sow.discount_rate) .^ years_idx
    ead_npv = sum(eads .* discount_fracs)
    costs_npv = sum(costs .* discount_fracs)
    return -(ead_npv + costs_npv)
end