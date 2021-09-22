##################################################################
# function to simulate one iteration of the experiment. For a    #
# given number of beetles, function determines if and when each  #
# beetle dies.                                                   #
##################################################################

# total_beetle: positive integer, refers to total number of beetles in experiment
# α_i, β_i, γ_i, p_i: positive reals, parameter values
# z_0: positive real, initial dose of insecticide

# returns length 14 vector where first 13 values correspond
# to the number of beetles that died on that day and final
# value corresponds to the number of beetles that survived.
function beetle_sim(total_beetle, α_i, β_i, γ_i, p_i, z_0)

    results = zeros(14)

    S_prev = 1

    for i = 1:total_beetle
        # number of releases occuring up to or on 13th day
        V_count = rand(Poisson(α_i * 13))

        # if no releases, count beetle as survived
        if V_count == 0
            results[14] = results[14] + 1
            continue
        end

        # draw inter-arrival times for releases
        T_inter = rand(Exponential(1/α_i), V_count)

        # set amount of insecticide
        z_i = z_0

        # death status of beetle
        death = false

        # time in loop begins at 0
        current_t = 0

        # calculate probability of death for each inter-release
        # period. Then remove amount of toxin released at end of
        # period and move on to next period
        for j = 1:(V_count + 1)

            # last iteration, for time between last release and
            # end of 13th day
            if j == V_count + 1

                # prob of death in inter-release interval
                death_p = z_i * γ_i + p_i * z_i *
                        (2*current_t*(13 - T_inter[j-1]) +
                        (13 - T_inter[j-1])^2) / 2

                if death_p > rand(Uniform(0,1))
                    # if the beetle dies, set death true and
                    # break out of for loop
                    death = true
                    # add to death counter of final interval
                    results[13] = results[13] + 1
                end

                break
            end

            # update time for next loop
            current_t = sum(T_inter[1:j])

            # make sure we haven't gone above 13 days
            if current_t > 13
                break
            end

            S_cur = exp(-(γ_i * current_t * z_i + z_i * p_i * current_t^2 / 2))

            death_p = S_prev - S_cur

            # set S_prev for next iteration
            S_prev = S_cur

            # use Uniform to generate random deaths
            if death_p > rand(Uniform(0,1))
                # if the beetle dies, set death true and
                # break out of for loop
                death = true
                # add to death counter
                index = convert(Int,round(current_t, RoundUp))

                results[index] = results[index] + 1
                break
            end

            # remove insecticide released at end of this period
            z_i = z_i - rand(Exponential(1 / β_i))

            # stop loop if negative insecticide
            if z_i <= 0
                break
            end

        end

        if death == false
            results[14] = results[14] + 1
        end


    end

    return results

end


function beetle_sim_ext(total_beetle, α_i, β_i, γ_i, p_i, z_0)

    results = zeros(14)

    for i = 1:total_beetle
        # number of releases occuring up to or on 13th day
        V_count = rand(Poisson(α_i * 13))
        if V_count == 0
            results[14] = results[14] + 1
            continue
        end

        # draw inter-arrival times for first release
        T_inter = rand(Exponential(1/α_i))

        # prob of death in inter-release interval
        death_p = p_i * z_0 *
                        (T_inter^2) / 2

        if death_p > rand(Uniform(0,1))
            # if the beetle dies, set death true and
            # break out of for loop
            death = true
            # add to death counter
            index = convert(Int,round(T_inter, RoundUp))
            if index > 14
                results[14] = results[14] + 1
                continue
            end

            results[index] = results[index] + 1
        else
            results[14] = results[14] + 1
        end

    end

    return results

end


# ABC simulation function for the beetle experiment. Takes as inputs
# original data from experiment, number of iterations, and initial
# dose of insecticide. For each iteration it records the parameter
# values drawn from the priors and the distance of the simulated
# data from the observed data. Returns a dataframe with those values
# for each iteration.


function abc_sim(data, z_0, n, beetle_n, α_prior, β_prior, γ_prior, p_prior)
    results = DataFrame( α = Float64[], β = Float64[], γ = Float64[],
                p = Float64[], dist = Float64[])

    # these lines put data in same format given in Diggle 1984
    sum_1 = sum(data[6:7])
    sum_2 = sum(data[8:13])
    data = vcat(data[1:5], sum_1, sum_2, data[14])
    for i = 1:n

        # draw instances of parameters from prior distributions
        α_i = rand(α_prior)
        β_i = rand(β_prior) + 5
        γ_i = rand(γ_prior)
        p_i = rand(p_prior)

        # run the experiment simulation
        sim_i = beetle_sim(beetle_n, α_i, β_i, γ_i, p_i,z_0)

        # convert data into form from Diggle 1984
        sum_i1 = sum(sim_i[6:7])
        sum_i2 = sum(sim_i[8:13])
        sim_i  = vcat(sim_i[1:5], sum_i1, sum_i2, sim_i[14])

        # calculate distance b/t simulated and observed data
        dist = norm(sim_i - data, 1)

        # record results
        push!(results, [α_i, β_i, γ_i, p_i, dist])

    end
    return results
end
