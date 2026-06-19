using Plots
using GLMakie
using Statistics
using Distributions
using Random
"""
    Molecule

Class for storing information about a molecule.

### Fields
- `position::Float64`: The molecule's position [m].
- `velocity::Float64`: The molecule's velocity [m/s].

- `mass::Float64`: The mass of the molecule [kg].
- `radius::Float64`: The radius of the molecule [m].
- `chemical_formula::String`: The chemical formula of the molecule.

- `positions_history::Vector{Vector{Float64}}`: A vector storing the history of the molecule's position.
- `velocities_history::Vector{Vector{Float64}}`: A vector storing the history of the molecule's velocity.
"""
mutable struct Molecule
    position::Vector{Float64}
    velocity::Vector{Float64}

    mass::Float64
    radius::Float64
    chemical_formula::String

    g::Vector{Float64}

    positions_history::Vector{Vector{Float64}}
    velocities_history::Vector{Vector{Float64}}
end

"""
    Domain

Class for storing information about the simulation domain.

### Fields
- `lx::Tuple{Float64, Float64}`: The bounds of the domain in the x-direction [m].
- `ly::Tuple{Float64, Float64}`: The bounds of the domain in the y-direction [m].
- `lz::Tuple{Float64, Float64}`: The bounds of the domain in the z-direction [m].
"""
mutable struct Domain
    lx::Tuple{Float64, Float64}
    ly::Tuple{Float64, Float64}
    lz::Tuple{Float64, Float64}
end

"""
    domainVolume(domain::Domain)

Calculate the volume of the simulation domain.

### Fields
- `domain::Domain`: The simulation domain for which to calculate the volume.
"""
function domainVolume(domain::Domain)
    return abs(domain.lx[2] - domain.lx[1]) * abs(domain.ly[2] - domain.ly[1]) * abs(domain.lz[2] - domain.lz[1])
end

"""
    checkDomain(molecule::Molecule, domain::Domain, temperature_top::Float64, temperature_bottom::Float64, add_temperature_gradient::Bool)

Check if a molecule is within the bounds of the simulation domain and update its position and velocity if it collides with the walls.

### Fields
- `molecule::Molecule`: The molecule to check.
- `domain::Domain`: The simulation domain to check against.
- `temperature_top::Float64`: The temperature of the top wall (used if `add_temperature_gradient` is true) [K].
- `temperature_bottom::Float64`: The temperature of the bottom wall (used if `add_temperature_gradient` is true) [K].
- `add_temperature_gradient::Bool`: A boolean indicating whether to add a temperature gradient to the simulation. If true, the molecule's velocity will be updated based on the temperature of the wall it collides with.
"""
function checkDomain(molecule::Molecule,domain::Domain, temperature_top::Float64, temperature_bottom::Float64, add_temperature_gradient::Bool)
    domain_pos = [domain.lx, domain.ly, domain.lz]

    if add_temperature_gradient
        for dim in eachindex(molecule.position)
            if domain_pos[dim][1] == domain_pos[dim][2]
                continue
            end

            if dim != 3
                # left,right,front and back walls
                # left wall
                if molecule.position[dim] - molecule.radius < domain_pos[dim][1]
                    dist = domain_pos[dim][1] - (molecule.position[dim] - molecule.radius)
                    
                    molecule.position[dim] = molecule.position[dim] + 2*dist
                    molecule.velocity[dim] = molecule.velocity[dim] * -1.0

                # right wall
                elseif molecule.position[dim] + molecule.radius > domain_pos[dim][2]
                    dist = (molecule.position[dim] + molecule.radius) - domain_pos[dim][2]

                    molecule.position[dim] = molecule.position[dim] - 2*dist
                    molecule.velocity[dim] = molecule.velocity[dim] * -1.0
                end
            else
                # top and bottom walls
                kb = 1.380649e-23
                # bottom wall
                if molecule.position[dim] - molecule.radius < domain_pos[dim][1]
                    dist = domain_pos[dim][1] - (molecule.position[dim] - molecule.radius)
                    
                    sigma = sqrt((temperature_bottom * kb) / molecule.mass)

                    vx = rand(Normal(0, sigma))
                    vy = rand(Normal(0, sigma))
                    vz = sigma * sqrt(-2 * log(rand()))

                    molecule.position[dim] = molecule.position[dim] + 2*dist
                    molecule.velocity = [vx, vy, vz]

                # top wall
                elseif molecule.position[dim] + molecule.radius > domain_pos[dim][2]
                    dist = (molecule.position[dim] + molecule.radius) - domain_pos[dim][2]

                    sigma = sqrt((temperature_top * kb) / molecule.mass)

                    vx = rand(Normal(0, sigma))
                    vy = rand(Normal(0, sigma))
                    vz = -sigma * sqrt(-2 * log(rand()))

                    molecule.position[dim] = molecule.position[dim] - 2*dist
                    molecule.velocity = [vx, vy, vz]
                end
            end
        end
    else
        for dim in eachindex(molecule.position)
            if domain_pos[dim][1] == domain_pos[dim][2]
                continue
            end
            # left wall
            if molecule.position[dim] - molecule.radius < domain_pos[dim][1]
                dist = domain_pos[dim][1] - (molecule.position[dim] - molecule.radius)
                
                molecule.position[dim] = molecule.position[dim] + 2*dist
                molecule.velocity[dim] = molecule.velocity[dim] * -1.0

            # right wall
            elseif molecule.position[dim] + molecule.radius > domain_pos[dim][2]
                dist = (molecule.position[dim] + molecule.radius) - domain_pos[dim][2]

                molecule.position[dim] = molecule.position[dim] - 2*dist
                molecule.velocity[dim] = molecule.velocity[dim] * -1.0
            end
        end
    end
end

"""
    ComputeNextPosition(molecule::Molecule, dt::Float64)

Update the position of a molecule based on its current velocity and a time step.

# Fields
- `molecule::Molecule`: The molecule whose position will be updated.
- `dt::Float64`: The time step for the position update.
"""
function ComputeNextPosition(molecule::Molecule, dt::Float64)
    molecule.velocity = molecule.velocity + molecule.g .* dt
    molecule.position = molecule.position + dt .* molecule.velocity
end

function simulation(position::Vector{Vector{Float64}}, velocity::Vector{Vector{Float64}}, mass::Vector{Float64}, radius::Vector{Float64}, chemical_formula::Vector{String}, number_of_steps::Int64, dt::Float64, domain::Domain, g::Vector{Float64}, remove_wall::Bool, new_domain::Domain, temperature_top::Float64, temperature_bottom::Float64, add_temperature_gradient::Bool)
    molecules::Vector{Molecule} = Molecule[]    
    current_domain = domain
    
    for i in 1:length(position)
        push!(molecules,Molecule(position[i], velocity[i], mass[i], radius[i], chemical_formula[i], g, [zeros(Float64,length(position[i])) for _ in 1:number_of_steps], [zeros(Float64,length(position[i])) for _ in 1:number_of_steps]))
    end

    for m in molecules
        m.positions_history[1] .= m.position
        m.velocities_history[1] .= m.velocity
    end

    for t in 2:number_of_steps

        if t % 1000 == 0
            println("Step: ", t, "/", number_of_steps)
        end
        
        for m in molecules
            ComputeNextPosition(m, dt)
        end

        for m in molecules
            if remove_wall
                if t > div(number_of_steps, 2)
                    checkDomain(m, new_domain, temperature_top, temperature_bottom, add_temperature_gradient)
                else 
                    checkDomain(m, current_domain, temperature_top, temperature_bottom, add_temperature_gradient)
                end
            else
                checkDomain(m, current_domain, temperature_top, temperature_bottom, add_temperature_gradient)
            end
        end

        checkCollision(molecules)

        for m in molecules
            m.positions_history[t] .= m.position
            m.velocities_history[t] .= m.velocity 
        end
    end

    return molecules
end

function computeNextPositionLennardJones(molecule::Molecule, force::Vector{Float64}, dt::Float64)
    molecule.velocity = molecule.velocity + (force ./ molecule.mass) .* dt # F = m*a => a = F/m
    molecule.position = molecule.position + molecule.velocity .* dt
end

function simulationLennardJones(position::Vector{Vector{Float64}}, velocity::Vector{Vector{Float64}}, mass::Vector{Float64}, radius::Vector{Float64}, chemical_formula::Vector{String}, number_of_steps::Int64, dt::Float64, domain::Domain, g::Vector{Float64}, temperatures_at_time::Vector{Float64}, sigma::Float64, epsilon::Float64)
    molecules::Vector{Molecule} = Molecule[]

    for i in 1:length(position)
        push!(molecules, Molecule(position[i], velocity[i], mass[i], radius[i], chemical_formula[i], g, [zeros(Float64, length(position[i])) for _ in 1:number_of_steps], [zeros(Float64, length(position[i])) for _ in 1:number_of_steps]))
    end

    for m in molecules
        m.positions_history[1] .= m.position
        m.velocities_history[1] .= m.velocity
    end

    for t in 2:number_of_steps
        if t % 1000 == 0
            println("Step: ", t, "/", number_of_steps)
        end

        forces = [lennardJonesForce(molecules, i, sigma, epsilon) for i in eachindex(molecules)]

        for i in eachindex(molecules)
            computeNextPositionLennardJones(molecules[i], forces[i], dt)
        end

        for m in molecules
            checkDomain(m, domain, 0.0, 0.0, false)
        end

        for m in molecules
            m.positions_history[t] .= m.position
            m.velocities_history[t] .= m.velocity
        end

        velocityRescaling(molecules, temperatures_at_time[t], t)

        for m in molecules
            m.velocities_history[t] .= m.velocity
        end
    end

    return molecules
end


function detectCollision(molecule_a::Molecule, molecule_b::Molecule)
    dist_min = molecule_a.radius + molecule_b.radius

    dist = molecule_a.position .- molecule_b.position
    dist_norme = sqrt(sum(dist .^ 2))

    return dist_norme <= dist_min
end

function checkCollision(molecules::Vector{Molecule})
    for i in eachindex(molecules)
        for j in eachindex(molecules)
            if i > j && detectCollision(molecules[i],molecules[j])
                normal_vect = (molecules[i].position .- molecules[j].position) / sqrt(sum((molecules[i].position .- molecules[j].position).^2))
                factor_1 = (2*molecules[j].mass / (molecules[i].mass + molecules[j].mass)).* (sum((molecules[i].velocity .- molecules[j].velocity) .* normal_vect)) .* normal_vect
                factor_2 = (2*molecules[i].mass / (molecules[i].mass + molecules[j].mass)).* (sum((molecules[i].velocity .- molecules[j].velocity) .* normal_vect)) .* normal_vect
                
                molecules[i].velocity = molecules[i].velocity .- factor_1
                molecules[j].velocity = molecules[j].velocity .+ factor_2
            end
        end
    end
end

function lennardJonesForce(molecules::Vector{Molecule}, i::Int64, sigma::Float64, epsilon::Float64)
    force = zeros(Float64, length(molecules[i].position))

    for j in eachindex(molecules)
        if i == j
            continue
        end

        r_ij = molecules[i].position .- molecules[j].position
        r = sqrt(sum(r_ij .^ 2))

        if r < 3 * sigma
            force .+= 24 * epsilon / r^2 * (2 * (sigma/r)^12 - (sigma/r)^6) .* r_ij
        end
    end

    return force
end

function velocityRescaling(molecules::Vector{Molecule}, T_ref::Float64, t::Int64)
    T_cal = calcTemperature(molecules, t, 2)
    factor = sqrt(T_ref / T_cal)
    for m in molecules
        m.velocity .*= factor
    end
end

function energyStableValue(molecules::Vector{Molecule}, t::Int64, window_size::Int64)
    if t <= window_size
        return false
    end

    energies = [calcEmec(molecules, t_i) for t_i in (t-window_size):t]
    energies_mean = mean(energies)
    energies_std = std(energies)

    energies_std / energies_mean
end

function plotSystem(molecules::Vector{Molecule}, t::Int64, domain::Domain)
    # camera=(0, 0)
    Plots.plot(legend=false,xlims=(domain.lx[1],domain.lx[2]),ylims=(domain.ly[1],domain.ly[2]),zlims=(domain.lz[1],domain.lz[2]))

    for m in molecules
        x = [m.positions_history[t][1]]
        y = [m.positions_history[t][2]]
        z = [m.positions_history[t][3]]
        Plots.scatter!(x,y,z,markersize=3)
    end
end

function plotSystem2D(molecules::Vector{Molecule}, t::Int64, domain::Domain, number_of_steps::Int64, temperature::Float64)
    Plots.plot(
        legend=false,
        xlims=(domain.lx[1],domain.lx[2]), 
        ylims=(domain.ly[1],domain.ly[2]),
        aspect_ratio=:equal,
        size=(1000, 1000),
        title = "Step $t / $number_of_steps | T = $(round(temperature, digits=2)) K",
        background_color=:black,
        foreground_color=:white)

    for m in molecules
        x = [m.positions_history[t][1]]
        y = [m.positions_history[t][2]]
        Plots.scatter!(x,y,markersize=13,color=:cyan)
    end
end

function calcEmec(molecules::Vector{Molecule}, t::Int64)
    emec = 0

    for m in molecules
        emec += 1/2 * m.mass * (sum(m.velocities_history[t] .^ 2))
    end

    return emec
end

function plotEmec(molecules::Vector{Molecule})
    val = calcEmec(molecules, 1)
    y = [calcEmec(molecules, t) for t in 1:length(molecules[1].velocities_history)]
    p = Plots.plot([1:length(molecules[1].velocities_history)], y, title="mechanical energy of the system", grid=false, legend=false, xlabel="time [s]", ylabel="Mechanical energy [J]", ylims=:auto, xlims=:auto, xticks=:auto, yticks=range(minimum(y), maximum(y), length=5))
    display(p)
end

function calcQuantityOfMovement(molecules::Vector{Molecule}, t::Int64)
    p::Vector{Float64} = zeros(Float64,length(molecules[1].velocity))

    for m in molecules
        p .+= m.mass .* m.velocities_history[t]
    end

    return p
end

function plotQuantityOfMovement(molecules::Vector{Molecule})
    val = calcQuantityOfMovement(molecules, 1)

    for dim in 1:length(molecules[1].velocity)
        y = [calcQuantityOfMovement(molecules, t)[dim] for t in 1:length(molecules[1].velocities_history)]
        p = Plots.plot([1:length(molecules[1].velocities_history)], y, title="quantity of movement for axis $dim", grid=false, legend=false, xlabel="time [s]", ylabel="Quantity of movement (axis $dim) [kg*m/s]", ylims=:auto, xlims=:auto, xticks=:auto, yticks=range(minimum(y), maximum(y), length=5))
        display(p)
    end
end

function calcMeanVelocity(molecules::Vector{Molecule}, t::Int64)
    velocities::Vector{Float64} = []

    for m in molecules
        push!(velocities, sqrt(sum(m.velocities_history[t] .^ 2)))
    end

    return mean(velocities)
end

function plotMeanVelocity(molecules::Vector{Molecule})
    y = [calcMeanVelocity(molecules, t) for t in 1:length(molecules[1].velocities_history)]
    p = Plots.plot([1:length(molecules[1].velocities_history)], y, title="mean velocity over time", grid=false, legend=false, xlabel="time [s]", ylabel="velocity [m/s]", ylims=:auto, xlims=:auto, xticks=:auto, yticks=range(minimum(y), maximum(y), length=5))
    display(p)
end

function plotVelocityDistributionFinal(molecules::Vector{Molecule})
    t_final = length(molecules[1].velocities_history)

    val = [sqrt(sum(m.velocities_history[t_final] .^ 2)) for m in molecules]
    p = Plots.histogram(val, bins = 30, title="final velocity magnitude distribution", grid=false, legend=false, xlabel="velocity value [m/s]", ylabel="number of molecules", ylims=:auto, xlims=:auto)
    display(p)
end

function calcMeanVelocitySquare(molecules::Vector{Molecule}, t::Int64)
    velocities::Vector{Float64} = []

    for m in molecules
        push!(velocities, m.mass * sum(m.velocities_history[t] .^ 2))
    end

    return mean(velocities)
end

function calcTemperature(molecules::Vector{Molecule}, t::Int64, dims::Int64 = 3)
    kb = 1.380649e-23

    # mass = molecules[1].mass
    mean_velocity = calcMeanVelocitySquare(molecules, t)

    # temperature = (mass * mean_velocity) /  (dims*kb)
    temperature = (mean_velocity) /  (dims*kb)

    return temperature
end

function plotTemperature(molecules::Vector{Molecule})
    y = [calcTemperature(molecules, t) for t in 1:length(molecules[1].velocities_history)]
    p = Plots.plot([1:length(molecules[1].velocities_history)], y, title="temperature over time", grid=false, legend=false, xlabel="time [s]", ylabel="Temperature [K]", ylims=:auto, xlims=:auto, xticks=:auto, yticks=range(minimum(y), maximum(y), length=5))
    display(p)
end

function calcPressure(molecules::Vector{Molecule}, t::Int64, domain::Domain)
    number_molecules = length(molecules) 
    # mass = molecules[1].mass
    mean_velocity = calcMeanVelocitySquare(molecules, t)
    domain_volume = domainVolume(domain)

    # pressure = (number_molecules * mass * mean_velocity) /  (3 * domain_volume)
    pressure = (number_molecules * mean_velocity) /  (3 * domain_volume)

    return pressure
end

function plotPressure(molecules::Vector{Molecule}, domain::Domain)
    y = [calcPressure(molecules, t, domain) for t in 1:length(molecules[1].velocities_history)]
    p = Plots.plot([1:length(molecules[1].velocities_history)], y, title="pressure over time", grid=false, legend=false, xlabel="time [s]", ylabel="Pressure [Pa]", ylims=:auto, xlims=:auto, xticks=:auto, yticks=range(minimum(y), maximum(y), length=5))
    display(p)
end


function plotPositionZDistribution(molecules::Vector{Molecule}, t::Int64)
    val = [m.positions_history[t][3] for m in molecules]
    p = Plots.histogram(val, bins = 20, title="final position z distribution", grid=false, legend=false, normalize=:probability, xlabel="position z [m]", ylabel="percentage", ylims=:auto, xlims=:auto)
    display(p)
end

function plotPositionZDistributionMultiSpecies(molecules::Vector{Molecule}, t::Int64)
    species = unique([m.chemical_formula for m in molecules])
    colors = [:red, :blue, :green, :orange, :purple, :cyan, :magenta]

    p = Plots.plot(title="final position z distribution by species", grid=false, legend=:topright, normalize=:probability, xlabel="position z [m]", ylabel="percentage", ylims=:auto, xlims=:auto)

    for (i, s) in enumerate(species)
        val = [m.positions_history[t][3] for m in molecules if m.chemical_formula == s]
        Plots.histogram!(p, val, bins=30, label=s, color=colors[i], alpha=0.5)
    end

    display(p)
end

function plotSpatialDistributionXY(molecules::Vector{Molecule}, t_range::AbstractRange, domain::Domain, phase_name::String)
    number_bins = 25

    dx = (domain.lx[2] - domain.lx[1]) / number_bins
    dy = (domain.ly[2] - domain.ly[1]) / number_bins

    counts = zeros(Float64, number_bins, number_bins)
    for t in t_range
        for m in molecules
            x = m.positions_history[t][1]
            y = m.positions_history[t][2]
            i = clamp(Int(floor((x - domain.lx[1]) / dx)) + 1, 1, number_bins)
            j = clamp(Int(floor((y - domain.ly[1]) / dy)) + 1, 1, number_bins)
            counts[j, i] += 1
        end
    end
    counts ./= sum(counts)

    x_centers = [domain.lx[1] + (i - 0.5) * dx for i in 1:number_bins]
    y_centers = [domain.ly[1] + (j - 0.5) * dy for j in 1:number_bins]

    ph = Plots.heatmap(x_centers, y_centers, counts,
        color = :inferno, background_color_inside = :black,
        title = "spatial distribution in xy plane - $phase_name",
        grid = false, xlabel = "position x [m]", ylabel = "position y [m]",
        xlims = (domain.lx[1], domain.lx[2]), ylims = (domain.ly[1], domain.ly[2]))
    display(ph)
end

function calcPressureZDistribution(molecules::Vector{Molecule}, t::Int64, domain::Domain, number_bins::Int64)
    bin_size = abs(domain.lz[2] - domain.lz[1]) / number_bins

    molecules_in_bin::Vector{Vector{Molecule}} = [Molecule[] for _ in 1:number_bins]
    z_values = [domain.lz[1] + (i-0.5)*bin_size for i in 1:number_bins]
    pressure_values = zeros(Float64, number_bins)

    for i in 1:number_bins
        for m in molecules
            if m.positions_history[t][3] >= domain.lz[1] + bin_size*(i-1) && m.positions_history[t][3] < domain.lz[1] + bin_size*i
                push!(molecules_in_bin[i], m)
            end
        end

        if length(molecules_in_bin[i]) > 0
            pressure_values[i] = calcPressure(molecules_in_bin[i], t, Domain(domain.lx, domain.ly, (domain.lz[1] + bin_size*(i-1), domain.lz[1] + bin_size*i)))
        end
    end
    
    return z_values, pressure_values
end

function calcTemperatureZDistribution(molecules::Vector{Molecule}, t::Int64, domain::Domain, number_bins::Int64)
    bin_size = abs(domain.lz[2] - domain.lz[1]) / number_bins

    molecules_in_bin::Vector{Vector{Molecule}} = [Molecule[] for _ in 1:number_bins]
    z_values = [domain.lz[1] + (i-0.5)*bin_size for i in 1:number_bins]
    temperature_values = zeros(Float64, number_bins)

    for i in 1:number_bins
        for m in molecules
            if m.positions_history[t][3] >= domain.lz[1] + bin_size*(i-1) && m.positions_history[t][3] < domain.lz[1] + bin_size*i
                push!(molecules_in_bin[i], m)
            end
        end

        if length(molecules_in_bin[i]) > 0
            temperature_values[i] = calcTemperature(molecules_in_bin[i], t)
        end
    end
    
    return z_values, temperature_values
end

function calcMeanVelocitySquareZDistribution(molecules::Vector{Molecule}, t::Int64, domain::Domain, number_bins::Int64)
    bin_size = abs(domain.lz[2] - domain.lz[1]) / number_bins

    molecules_in_bin::Vector{Vector{Molecule}} = [Molecule[] for _ in 1:number_bins]
    z_values = [domain.lz[1] + (i-0.5)*bin_size for i in 1:number_bins]
    mean_velocity_square_values = zeros(Float64, number_bins)

    for i in 1:number_bins
        for m in molecules
            if m.positions_history[t][3] >= domain.lz[1] + bin_size*(i-1) && m.positions_history[t][3] < domain.lz[1] + bin_size*i
                push!(molecules_in_bin[i], m)
            end
        end

        if length(molecules_in_bin[i]) > 0
            mean_velocity_square_values[i] = calcMeanVelocitySquare(molecules_in_bin[i], t)
        end
    end
    
    return z_values, mean_velocity_square_values
end

function plotPressureZDistribution(molecules::Vector{Molecule}, domain::Domain, t::Int64)
    z_values, pressure_values = calcPressureZDistribution(molecules, t, domain, 20)

    p = Plots.bar(z_values, pressure_values, title="pressure distribution along z axis", grid=false, legend=false, xlabel="position z [m]", ylabel="Pressure [Pa]", ylims=:auto, xlims=:auto, xticks=:auto, yticks=range(minimum(pressure_values), maximum(pressure_values), length=5))
    display(p)
end

function plotTemperatureZDistribution(molecules::Vector{Molecule}, domain::Domain, t::Int64)
    z_values, temperature_values = calcTemperatureZDistribution(molecules, t, domain, 20)

    p = Plots.bar(z_values, temperature_values, title="temperature distribution along z axis", grid=false, legend=false, xlabel="position z [m]", ylabel="Temperature [K]", ylims=:auto, xlims=:auto, xticks=:auto, yticks=range(minimum(temperature_values), maximum(temperature_values), length=5))
    display(p)
end

function plotMeanVelocityZDistribution(molecules::Vector{Molecule}, domain::Domain, t::Int64)
    z_values, mean_velocity_square_values = calcMeanVelocitySquareZDistribution(molecules, t, domain, 20)

    p = Plots.bar(z_values, mean_velocity_square_values, title="mean velocity square distribution along z axis", grid=false, legend=false, xlabel="position z [m]", ylabel="mean velocity square [m^2/s^2]", ylims=:auto, xlims=:auto, xticks=:auto, yticks=range(minimum(mean_velocity_square_values), maximum(mean_velocity_square_values), length=5))
    display(p)
end

function calcEntropie(molecules::Vector{Molecule}, domain::Domain, t::Int64, domain_v_square::Vector{Float64} = [0, 200000], number_bins_x::Int64 = 10, number_bins_y::Int64 = 10, number_bins_z::Int64 = 10, number_bins_v::Int64 = 200)
    bin_x_size = abs(domain.lx[2] - domain.lx[1]) / number_bins_x
    bin_y_size = abs(domain.ly[2] - domain.ly[1]) / number_bins_y
    bin_z_size = abs(domain.lz[2] - domain.lz[1]) / number_bins_z
    bin_v_size = (domain_v_square[2] - domain_v_square[1]) / number_bins_v

    number_molecules = length(molecules)
    molecules_x = zeros(Float64, number_bins_x)
    molecules_y = zeros(Float64, number_bins_y)
    molecules_z = zeros(Float64, number_bins_z)
    molecules_v = zeros(Float64, number_bins_v)

    for i in 1:number_bins_v
        for m in molecules
            v_square = sqrt(sum(m.velocities_history[t] .^ 2))
            if v_square >= domain_v_square[1] + bin_v_size*(i-1) && v_square < domain_v_square[1] + bin_v_size*i
                molecules_v[i] += 1
            end
        end
    end

    proba_v = molecules_v ./ number_molecules
    
    for i in 1:number_bins_x
        for m in molecules
            if m.positions_history[t][1] >= domain.lx[1] + bin_x_size*(i-1) && m.positions_history[t][1] < domain.lx[1] + bin_x_size*i
                molecules_x[i] += 1
            end
        end
    end

    for i in 1:number_bins_y
        for m in molecules
            if m.positions_history[t][2] >= domain.ly[1] + bin_y_size*(i-1) && m.positions_history[t][2] < domain.ly[1] + bin_y_size*i
                molecules_y[i] += 1
            end
        end
    end

    for i in 1:number_bins_z
        for m in molecules
            if m.positions_history[t][3] >= domain.lz[1] + bin_z_size*(i-1) && m.positions_history[t][3] < domain.lz[1] + bin_z_size*i
                molecules_z[i] += 1
            end
        end
    end

    proba_x = molecules_x ./ number_molecules
    proba_y = molecules_y ./ number_molecules   
    proba_z = molecules_z ./ number_molecules

    entropie_x = -sum([p > 0 ? p * log(p) : 0.0 for p in proba_x])
    entropie_y = -sum([p > 0 ? p * log(p) : 0.0 for p in proba_y])
    entropie_z = -sum([p > 0 ? p * log(p) : 0.0 for p in proba_z])

    entropie_v = -sum([p > 0 ? p * log(p) : 0.0 for p in proba_v])

    return entropie_x + entropie_y + entropie_z + entropie_v
end


function plotEntropie(molecules::Vector{Molecule}, domain::Domain)
    y = [calcEntropie(molecules, domain, t, [0.0, 200000.0], 10, 10, 10, 200) for t in 1:length(molecules[1].velocities_history)]
    p = Plots.plot([1:length(molecules[1].velocities_history)], y, title="entropie of the system", grid=false, legend=false, xlabel="time [s]", ylabel="Entropie ", ylims=:auto, xlims=:auto, xticks=:auto, yticks=range(minimum(y), maximum(y), length=5))
    display(p)
end

function plotEntropieRemoveWall(molecules::Vector{Molecule}, domain::Domain, new_domain::Domain, number_of_steps::Int64)
    y1 = [calcEntropie(molecules, domain, t, [0.0, 200000.0], 10, 10, 10, 200) for t in 1:div(number_of_steps, 2)]
    y2 = [calcEntropie(molecules, new_domain, t, [0.0, 200000.0], 20, 10, 10, 200) for t in div(number_of_steps, 2)+1:number_of_steps]
    p = Plots.plot([1:number_of_steps], vcat(y1, y2), title="entropie of the system with wall removal", grid=false, legend=false, xlabel="time [s]", ylabel="Entropie ", ylims=:auto, xlims=:auto, xticks=:auto, yticks=range(minimum(vcat(y1, y2)), maximum(vcat(y1, y2)), length=5))
    display(p)
end

function makieSystem(molecules::Vector{Molecule}, domain::Domain, number_of_steps::Int64)

    fig = Figure()
    
    ax = Axis3(fig[1,1],limits = (domain.lx[1], domain.lx[2], domain.ly[1], domain.ly[2], domain.lz[1], domain.lz[2]), xgridvisible = false, ygridvisible = false, zgridvisible = false)
    
    ax.azimuth[] = pi/2
    ax.elevation[] = 0.0

    positions = Observable([Point3f(m.position[1], m.position[2], m.position[3]) for m in molecules])

    sizes = [m.radius for m in molecules]

    molecules_colors = []

    for m in molecules
        if m.chemical_formula == "He"
            push!(molecules_colors, :red)
        elseif m.chemical_formula == "Ar"
            push!(molecules_colors, :blue)
        else
            push!(molecules_colors, :green)
        end
    end
    
    meshscatter!(ax, positions, markersize = sizes, color = molecules_colors)

    return fig, positions
end

function makieSystemInteractive(molecules::Vector{Molecule}, domain::Domain, number_of_steps::Int64)

    fig = Figure()
    
    ax = Axis3(fig[1,1],limits = (domain.lx[1], domain.lx[2], domain.ly[1], domain.ly[2], domain.lz[1], domain.lz[2]), xgridvisible = false, ygridvisible = false, zgridvisible = false)
    ax_distribution_z = Axis(fig[3,1], title="Position z distribution", xlabel="position z [m]", ylabel="percentage", xgridvisible = false, ygridvisible = false)
    ax_pressure_z = Axis(fig[4,1], title="Pressure distribution along z axis", xlabel="position z [m]", ylabel="Pressure [Pa]", xgridvisible = false, ygridvisible = false)

    slider = Slider(fig[2,1], range = 1:number_of_steps, startvalue = 1)
    button = Button(fig[2,2], label = "Play")

    molecules_colors = []

    for m in molecules
        if m.chemical_formula == "He"
            push!(molecules_colors, :red)
        elseif m.chemical_formula == "Ar"
            push!(molecules_colors, :blue)
        else
            push!(molecules_colors, :green)
        end
    end
    
    positions = Observable([Point3f(m.position[1], m.position[2], m.position[3]) for m in molecules])
    sizes = [m.radius for m in molecules]
    meshscatter!(ax, positions, markersize = sizes, color = molecules_colors)

    # z_distribution = Observable([m.positions_history[1][3] for m in molecules])
    # hist!(ax_distribution_z, z_distribution, bins = 20, normalization = :probability)

    z_distribution_Ar = Observable([m.positions_history[1][3] for m in molecules if m.chemical_formula == "Ar"])
    z_distribution_He = Observable([m.positions_history[1][3] for m in molecules if m.chemical_formula == "He"])

    z_distribution = [z_distribution_Ar, z_distribution_He]

    chimical_formulas_distribution = ["Ar", "He"]
    hist!(ax_distribution_z, z_distribution[1], bins = 20, normalization = :probability, color = :blue, label = chimical_formulas_distribution[1], alpha = 0.3)
    hist!(ax_distribution_z, z_distribution[2], bins = 20, normalization = :probability, color = :red, label = chimical_formulas_distribution[2], alpha = 0.3)

    z_values = Observable(Float64[])
    z_pressure_values = Observable(Float64[])

    z_init, p_init = calcPressureZDistribution(molecules, 1, domain, 20)

    z_values[] = z_init
    z_pressure_values[] = p_init

    barplot!(ax_pressure_z, z_values, z_pressure_values)

    return fig, positions, slider, button, z_distribution, z_values, z_pressure_values, chimical_formulas_distribution
end

function makieGetPositions(molecules, t)
    return [Point3f(m.positions_history[t][1],m.positions_history[t][2],m.positions_history[t][3]) for m in molecules]
end

function newMoleculePosition(domain::Domain, existing_positions::Vector{Vector{Float64}}, sigma::Float64)::Vector{Float64}
    min_dist = 0.9 * sigma

    while true
        pos = [
            rand() * (domain.lx[2] - domain.lx[1]) + domain.lx[1],
            rand() * (domain.ly[2] - domain.ly[1]) + domain.ly[1],
            rand() * (domain.lz[2] - domain.lz[1]) + domain.lz[1]
        ]

        valid = true
        for existing in existing_positions
            if sqrt(sum((pos .- existing) .^ 2)) <= min_dist
                valid = false
                break
            end
        end

        if valid
            return pos
        end
    end
end

function temperatureProfile(number_of_steps::Int64, t1::Float64, t2::Float64, n1::Int64, n2::Int64)                                                                            
    temperatures = zeros(Float64, number_of_steps)                                                                                                                             
                                                                                                                                                                               
    for t in 1:number_of_steps                                                                                                                                                 
        if t <= n1                                                                                                                                                             
            temperatures[t] = t1                                                                                                                                               
        elseif t <= n1 + n2                                                                                                                                                    
            temperatures[t] = t1 - (t1 - t2) * (t - n1) / n2                                                                                                                   
        else                                                                                                                                                                   
            temperatures[t] = t2                                                                                                                                               
        end                                                                                                                                                                    
    end                                                                                                                                                                        
                                                                                                                                                                               
    return temperatures
end

function main_lennardJones_condensation()
    dt::Float64 = 1.0e-15
    t_final::Float64 = 5.0e-11
    number_of_steps::Int64 = div(t_final, dt) + 1

    FPS = 30

    g = [0.0,0.0,0.0]


    positions::Vector{Vector{Float64}} = []
    velocities::Vector{Vector{Float64}} = []
    masses::Vector{Float64} = []
    radius::Vector{Float64} = []
    chemical_formulas::Vector{String} = []

    domain::Domain = Domain((-5e-9, 5e-9), (-5e-9, 5e-9), (0, 0))

    temperature_ref::Float64 = 40 # [K]
    temperature_final::Float64 = 10 # [K]
    N1_steps::Int64 = 15000
    N2_steps::Int64 = 15000
    temperatures_at_time::Vector{Float64} = temperatureProfile(number_of_steps, temperature_ref, temperature_final, N1_steps, N2_steps)
    
    sigma::Float64 = 2.74e-10 # [m]
    epsilon::Float64 = 4.91511044e-22 # [J]

    # Random generation

    kb = 1.380649e-23
    sigma_v = sqrt(kb * temperature_ref / 3.35105e-26)
    number_atomes::Int64 = 100

    spawn_domain::Domain = Domain((-5e-9, 5e-9), (-5e-9, 5e-9), (0, 0))

    for i in 1:number_atomes
        push!(positions, newMoleculePosition(spawn_domain, positions, sigma))

        velocity::Vector{Float64} = [rand(Normal(0, sigma_v)), rand(Normal(0, sigma_v)), 0.0]

        push!(velocities, velocity)
        push!(masses, 3.35105e-26)
        push!(radius, 1.37e-10)
        push!(chemical_formulas, "Ne")
    end

    @assert length(positions) == length(velocities) == length(masses) == length(radius) == length(chemical_formulas)

    # Simulation 2D - Lennard-Jones
    molecules::Vector{Molecule} = simulationLennardJones(positions,velocities,masses,radius,chemical_formulas, 
                                            number_of_steps, dt, domain, g, temperatures_at_time, sigma, epsilon)

    filename = "results/molecule_lennardJones_$temperature_ref - $temperature_final.mp4"

    step_anim = FPS * 10
    frames = unique(vcat(collect(1:step_anim:number_of_steps), number_of_steps))
    animation = @animate for t in frames
        plotSystem2D(molecules, t, domain, number_of_steps, temperatures_at_time[t])
    end

    mp4(animation, filename, fps = FPS)

    phase_gas_range = 1:N1_steps
    phase_transition_range = (N1_steps + 1):(N1_steps + N2_steps)
    phase_solid_range = (N1_steps + N2_steps + 1):number_of_steps

    plotSpatialDistributionXY(molecules, phase_gas_range, domain, "phase gaz")
    plotSpatialDistributionXY(molecules, phase_transition_range, domain, "phase transition")
    plotSpatialDistributionXY(molecules, phase_solid_range, domain, "phase solide")
end

function main_lennardJones()
    dt::Float64 = 1.0e-15
    t_final::Float64 = 5.0e-11
    number_of_steps::Int64 = div(t_final, dt) + 1

    FPS = 30

    g = [0.0,0.0,0.0]


    positions::Vector{Vector{Float64}} = []
    velocities::Vector{Vector{Float64}} = []
    masses::Vector{Float64} = []
    radius::Vector{Float64} = []
    chemical_formulas::Vector{String} = []

    domain::Domain = Domain((-5e-9, 5e-9), (-5e-9, 5e-9), (0, 0))

    temperature_ref::Float64 = 40 # [K]
    temperature_final::Float64 = 40 # [K]
    N1_steps::Int64 = 15000
    N2_steps::Int64 = 15000
    temperatures_at_time::Vector{Float64} = temperatureProfile(number_of_steps, temperature_ref, temperature_final, N1_steps, N2_steps)
    
    sigma::Float64 = 2.74e-10 # [m]
    epsilon::Float64 = 4.91511044e-22 # [J]

    # Random generation

    kb = 1.380649e-23
    sigma_v = sqrt(kb * temperature_ref / 3.35105e-26)
    number_atomes::Int64 = 100

    spawn_domain::Domain = Domain((-5e-9, 5e-9), (-5e-9, 5e-9), (0, 0))

    for i in 1:number_atomes
        push!(positions, newMoleculePosition(spawn_domain, positions, sigma))

        velocity::Vector{Float64} = [rand(Normal(0, sigma_v)), rand(Normal(0, sigma_v)), 0.0]

        push!(velocities, velocity)
        push!(masses, 3.35105e-26)
        push!(radius, 1.37e-10)
        push!(chemical_formulas, "Ne")
    end

    @assert length(positions) == length(velocities) == length(masses) == length(radius) == length(chemical_formulas)

    # Simulation 2D - Lennard-Jones
    molecules::Vector{Molecule} = simulationLennardJones(positions,velocities,masses,radius,chemical_formulas, 
                                            number_of_steps, dt, domain, g, temperatures_at_time, sigma, epsilon)

    filename = "results/molecule_lennardJones_$temperature_ref - $temperature_final.mp4"

    step_anim = FPS * 10
    frames = unique(vcat(collect(1:step_anim:number_of_steps), number_of_steps))
    animation = @animate for t in frames
        plotSystem2D(molecules, t, domain, number_of_steps, temperatures_at_time[t])
    end

    mp4(animation, filename, fps = FPS)

    phase_gas_range = 1:N1_steps
    phase_transition_range = (N1_steps + 1):(N1_steps + N2_steps)
    phase_solid_range = (N1_steps + N2_steps + 1):number_of_steps

    plotSpatialDistributionXY(molecules, phase_gas_range, domain, "phase gaz")
    plotSpatialDistributionXY(molecules, phase_transition_range, domain, "phase transition")
    plotSpatialDistributionXY(molecules, phase_solid_range, domain, "phase solide")
end

function main_entropy(remove_wall::Bool = true)
    dt::Float64 = 1.0e-14
    t_final::Float64 = 10.0e-11
    number_of_steps::Int64 = div(t_final, dt)

    FPS = 30

    g = [0.0,0.0,0.0]


    positions::Vector{Vector{Float64}} = []
    velocities::Vector{Vector{Float64}} = []
    masses::Vector{Float64} = []
    radius::Vector{Float64} = []
    chemical_formulas::Vector{String} = []

    domain::Domain = Domain((-5e-9, 5e-9), (-5e-9, 5e-9), (-5e-9, 5e-9))
    new_domain::Domain = Domain((-5e-9, 15e-9), (-5e-9, 5e-9), (-5e-9, 5e-9))


    # Random generation

    velocityValue::Float64 = 1400 # [m/s]
    number_atomes::Int64 = 400

    spawn_domain::Domain = Domain((-5e-9, 0), (-5e-9, 0), (-5e-9, 0))

    for i in 1:number_atomes
        push!(positions, [
                rand() * (spawn_domain.lx[2] - spawn_domain.lx[1]) + spawn_domain.lx[1],
                rand() * (spawn_domain.ly[2] - spawn_domain.ly[1]) + spawn_domain.ly[1],
                rand() * (spawn_domain.lz[2] - spawn_domain.lz[1]) + spawn_domain.lz[1]
            ])

        velocity::Vector{Float64} = [rand()*10-5,rand()*10-5,rand()*10-5]
        velocity = velocity ./ sqrt(sum(velocity .^2))
        velocity = velocity .* velocityValue

        @assert isapprox(sqrt(sum(velocity .^ 2)), velocityValue; atol=1e-6)

        push!(velocities, velocity)
        push!(masses, 6.646e-27)
        push!(radius, 1.1e-10)
        push!(chemical_formulas, "He")
    end

    @assert length(positions) == length(velocities) == length(masses) == length(radius) == length(chemical_formulas)

    # Simulation 3D
    molecules::Vector{Molecule} = simulation(positions,velocities,masses,radius,chemical_formulas, 
                                            number_of_steps, dt, domain, g, remove_wall, new_domain, 
                                            0.0, 0.0, false)

    if remove_wall
        plotEntropieRemoveWall(molecules, domain, new_domain, number_of_steps)
    else
        plotEntropie(molecules, domain)
    end

    # Save the video
    fig, pos  = makieSystem(molecules, new_domain, number_of_steps)

    record(fig, "results/molecule.mp4", 1:number_of_steps) do t
        pos[] = makieGetPositions(molecules, t)
    end
end

function main_temp_gradient()
    dt::Float64 = 1.0e-14
    t_final::Float64 = 1.0e-10
    number_of_steps::Int64 = div(t_final, dt)

    FPS = 30

    g = [0.0,0.0,0.0]


    positions::Vector{Vector{Float64}} = []
    velocities::Vector{Vector{Float64}} = []
    masses::Vector{Float64} = []
    radius::Vector{Float64} = []
    chemical_formulas::Vector{String} = []

    domain::Domain = Domain((-2e-9, 2e-9), (-2e-9, 2e-9), (-2e-9, 2e-9))
    new_domain::Domain = Domain((-2e-9, 2e-9), (-2e-9, 2e-9), (-2e-9, 2e-9))
    
    temperature_top::Float64 = 700 # [K]
    temperature_bottom::Float64 = 300 # [K]

    # Random generation
    velocityValue::Float64 = 1400 # [m/s]
    number_atomes::Int64 = 500

    spawn_domain::Domain = Domain((-2e-9, 2e-9), (-2e-9, 2e-9), (-2e-9, 2e-9))

    for i in 1:number_atomes
        push!(positions, [
                rand() * (spawn_domain.lx[2] - spawn_domain.lx[1]) + spawn_domain.lx[1],
                rand() * (spawn_domain.ly[2] - spawn_domain.ly[1]) + spawn_domain.ly[1],
                rand() * (spawn_domain.lz[2] - spawn_domain.lz[1]) + spawn_domain.lz[1]
            ])

        velocity::Vector{Float64} = [rand()*10-5,rand()*10-5,rand()*10-5]
        velocity = velocity ./ sqrt(sum(velocity .^2))
        velocity = velocity .* velocityValue

        @assert isapprox(sqrt(sum(velocity .^ 2)), velocityValue; atol=1e-6)

        push!(velocities, velocity)
        push!(masses, 6.646e-27)
        push!(radius, 1.1e-10)
        push!(chemical_formulas, "He")
    end

    @assert length(positions) == length(velocities) == length(masses) == length(radius) == length(chemical_formulas)

    # Simulation 3D
    molecules::Vector{Molecule} = simulation(positions,velocities,masses,radius,chemical_formulas, 
                                            number_of_steps, dt, domain, g, false, new_domain, 
                                            temperature_top, temperature_bottom, true)

    plotEntropie(molecules, domain)
    plotTemperatureZDistribution(molecules, domain, length(molecules[1].positions_history))
    
    # Save the video
    fig, pos  = makieSystem(molecules, domain, number_of_steps)

    record(fig, "results/molecule.mp4", 1:number_of_steps) do t
        pos[] = makieGetPositions(molecules, t)
    end
end

function main_multi_species(check_stability::Bool = true)
    dt::Float64 = 1.0e-14
    t_final::Float64 = 10.0e-11
    number_of_steps::Int64 = div(t_final, dt)

    FPS = 30

    g = [0.0,0.0,-9.81e13]


    positions::Vector{Vector{Float64}} = []
    velocities::Vector{Vector{Float64}} = []
    masses::Vector{Float64} = []
    radius::Vector{Float64} = []
    chemical_formulas::Vector{String} = []

    domain::Domain = Domain((-1e-8, 1e-8), (-1e-8, 1e-8), (-1e-8, 1e-8))
    new_domain::Domain = Domain((-1e-8, 1e-8), (-1e-8, 1e-8), (-1e-8, 1e-8))
    
    # Random generation

    velocityValue::Float64 = 789.45 # [m/s]
    number_atomes::Int64 = 400
    spawn_domain::Domain = Domain((-1e-8, 1e-8), (-1e-8, 1e-8), (-1e-8, 1e-8))

    for i in 1:number_atomes
        push!(positions, [
                rand() * (spawn_domain.lx[2] - spawn_domain.lx[1]) + spawn_domain.lx[1],
                rand() * (spawn_domain.ly[2] - spawn_domain.ly[1]) + spawn_domain.ly[1],
                rand() * (spawn_domain.lz[2] - spawn_domain.lz[1]) + spawn_domain.lz[1]
            ])

        velocity::Vector{Float64} = [rand()*10-5,rand()*10-5,rand()*10-5]
        velocity = velocity ./ sqrt(sum(velocity .^2))
        velocity = velocity .* velocityValue

        @assert isapprox(sqrt(sum(velocity .^ 2)), velocityValue; atol=1e-6)

        push!(velocities, velocity)
        push!(masses, 6.646e-27)
        push!(radius, 1.1e-10)
        push!(chemical_formulas, "He")
    end

    velocityValue = 249.88 # [m/s]
    number_atomes = 200
    spawn_domain = Domain((-1e-8, 1e-8), (-1e-8, 1e-8), (-1e-8, 1e-8))

    for i in 1:number_atomes
        push!(positions, [
                rand() * (spawn_domain.lx[2] - spawn_domain.lx[1]) + spawn_domain.lx[1],
                rand() * (spawn_domain.ly[2] - spawn_domain.ly[1]) + spawn_domain.ly[1],
                rand() * (spawn_domain.lz[2] - spawn_domain.lz[1]) + spawn_domain.lz[1]
            ])

        velocity::Vector{Float64} = [rand()*10-5,rand()*10-5,rand()*10-5]
        velocity = velocity ./ sqrt(sum(velocity .^2))
        velocity = velocity .* velocityValue

        @assert isapprox(sqrt(sum(velocity .^ 2)), velocityValue; atol=1e-6)

        push!(velocities, velocity)
        push!(masses, 6.634e-26)
        push!(radius, 1.88e-10)
        push!(chemical_formulas, "Ar")
    end

    @assert length(positions) == length(velocities) == length(masses) == length(radius) == length(chemical_formulas)

    # Simulation 3D
    molecules::Vector{Molecule} = simulation(positions,velocities,masses,radius,chemical_formulas, 
                                            number_of_steps, dt, domain, g, false, new_domain, 
                                            0.0, 0.0, false)

    # Verification of the stability of the simulation
    if check_stability
        window_size = div(number_of_steps, 10)
        stability_value = 6e-3

        stability_values = []
        times_stability = []

        t_stable = 0

        for t in window_size+1:100:number_of_steps        
            push!(stability_values, energyStableValue(molecules, t, window_size))
            push!(times_stability, t)

            if stability_values[end] < stability_value && t_stable == 0
                t_stable = t
            end
        end

        p = Plots.plot(times_stability, stability_values, title="stability of the simulation", grid=false, legend=false, xlabel="time [s]", ylabel="energy stability value", ylims=:auto, xlims=:auto, xticks=:auto, yticks=range(minimum(stability_values), maximum(stability_values), length=5))
        Plots.plot!(p, times_stability, [stability_value for _ in times_stability], color=:red, label="stability threshold")
        display(p)
        println("The simulation is stable at time: ", t_stable > 0 ? t_stable : "not stable during the simulation")
    end


    plotEmec(molecules)
    plotQuantityOfMovement(molecules)
    plotMeanVelocity(molecules)
    plotVelocityDistributionFinal(molecules)
    plotTemperature(molecules)
    plotPressure(molecules, domain)
    plotPositionZDistribution(molecules, length(molecules[1].positions_history))
    plotPositionZDistributionMultiSpecies(molecules, length(molecules[1].positions_history))
    plotPressureZDistribution(molecules, domain, length(molecules[1].positions_history))
    plotTemperatureZDistribution(molecules, domain, length(molecules[1].positions_history))
    plotMeanVelocityZDistribution(molecules, domain, length(molecules[1].positions_history))
    plotEntropie(molecules, domain)

    # filename = "results/molecule.mp4"

    # animation = @animate for t in 1:number_of_steps
    #     plotSystem(molecules, t, domain)
    # end

    # mp4(animation, filename, fps = FPS)


    # With makie

    # Save the video
    fig, pos  = makieSystem(molecules, domain, number_of_steps)
    
    record(fig, "results/molecule.mp4", 1:300:number_of_steps) do t
        pos[] = makieGetPositions(molecules, t)
    end


    # Interactive animation with Makie

    # fig, pos, slider, button, z_distribution, z_values, z_pressure_values, chemical_formulas_distribution = makieSystemInteractive(molecules, domain, number_of_steps)

    # playing = Observable(false)
    # stop_animation = true

    # on(slider.value) do t
    #     pos[] = makieGetPositions(molecules, Int(t))
    #     z_distribution[1][] = [m.positions_history[Int(t)][3] for m in molecules if m.chemical_formula == chemical_formulas_distribution[1]]
    #     z_distribution[2][] = [m.positions_history[Int(t)][3] for m in molecules if m.chemical_formula == chemical_formulas_distribution[2]]
    #     z_values[], z_pressure_values[] = calcPressureZDistribution(molecules, Int(t), domain, 20)
    # end

    # on(button.clicks) do _
    #     playing[] = !playing[]
    # end

    # @async while stop_animation
    #     if playing[]
    #         slider.value[] = mod(slider.value[] , number_of_steps) + 1
    #     end
    #     sleep(1/FPS)
    # end
    

    # wait(display(fig))

    # stop_animation = false
end

function main_standard(check_stability::Bool = true)
    dt::Float64 = 1.0e-14
    t_final::Float64 = 10.0e-11
    number_of_steps::Int64 = div(t_final, dt)

    FPS = 30

    g = [0.0,0.0,-9.81e13]


    positions::Vector{Vector{Float64}} = []
    velocities::Vector{Vector{Float64}} = []
    masses::Vector{Float64} = []
    radius::Vector{Float64} = []
    chemical_formulas::Vector{String} = []

    domain::Domain = Domain((-1e-8, 1e-8), (-1e-8, 1e-8), (-1e-8, 1e-8))
    new_domain::Domain = Domain((-1e-8, 1e-8), (-1e-8, 1e-8), (-1e-8, 1e-8))
    
    # Random generation

    velocityValue::Float64 = 1367 # [m/s]
    number_atomes::Int64 = 500

    spawn_domain::Domain = Domain((-1e-8, 1e-8), (-1e-8, 1e-8), (-1e-8, 1e-8))

    for i in 1:number_atomes
        push!(positions, [
                rand() * (spawn_domain.lx[2] - spawn_domain.lx[1]) + spawn_domain.lx[1],
                rand() * (spawn_domain.ly[2] - spawn_domain.ly[1]) + spawn_domain.ly[1],
                rand() * (spawn_domain.lz[2] - spawn_domain.lz[1]) + spawn_domain.lz[1]
            ])

        velocity::Vector{Float64} = [rand()*10-5,rand()*10-5,rand()*10-5]
        velocity = velocity ./ sqrt(sum(velocity .^2))
        velocity = velocity .* velocityValue

        @assert isapprox(sqrt(sum(velocity .^ 2)), velocityValue; atol=1e-6)

        push!(velocities, velocity)
        push!(masses, 6.646e-27)
        push!(radius, 1.1e-10)
        push!(chemical_formulas, "He")
    end

    @assert length(positions) == length(velocities) == length(masses) == length(radius) == length(chemical_formulas)

    # Simulation 3D
    molecules::Vector{Molecule} = simulation(positions,velocities,masses,radius,chemical_formulas, 
                                            number_of_steps, dt, domain, g, false, new_domain, 
                                            0.0, 0.0, false)

    # Verification of the stability of the simulation
    if check_stability
        window_size = div(number_of_steps, 10)
        stability_value = 6e-3

        stability_values = []
        times_stability = []

        t_stable = 0

        for t in window_size+1:100:number_of_steps        
            push!(stability_values, energyStableValue(molecules, t, window_size))
            push!(times_stability, t)

            if stability_values[end] < stability_value && t_stable == 0
                t_stable = t
            end
        end

        p = Plots.plot(times_stability, stability_values, title="stability of the simulation", grid=false, legend=false, xlabel="time [s]", ylabel="energy stability value", ylims=:auto, xlims=:auto, xticks=:auto, yticks=range(minimum(stability_values), maximum(stability_values), length=5))
        Plots.plot!(p, times_stability, [stability_value for _ in times_stability], color=:red, label="stability threshold")
        display(p)
        println("The simulation is stable at time: ", t_stable > 0 ? t_stable : "not stable during the simulation")
    end


    plotEmec(molecules)
    plotQuantityOfMovement(molecules)
    plotMeanVelocity(molecules)
    plotVelocityDistributionFinal(molecules)
    plotTemperature(molecules)
    plotPressure(molecules, domain)
    plotPositionZDistribution(molecules, length(molecules[1].positions_history))
    plotPressureZDistribution(molecules, domain, length(molecules[1].positions_history))
    plotTemperatureZDistribution(molecules, domain, length(molecules[1].positions_history))
    plotMeanVelocityZDistribution(molecules, domain, length(molecules[1].positions_history))
    plotEntropie(molecules, domain)

    # filename = "results/molecule.mp4"

    # animation = @animate for t in 1:number_of_steps
    #     plotSystem(molecules, t, domain)
    # end

    # mp4(animation, filename, fps = FPS)


    # With makie

    # Save the video
    fig, pos  = makieSystem(molecules, domain, number_of_steps)
    
    record(fig, "results/molecule.mp4", 1:300:number_of_steps) do t
        pos[] = makieGetPositions(molecules, t)
    end


    # Interactive animation with Makie

    # fig, pos, slider, button, z_distribution, z_values, z_pressure_values, chemical_formulas_distribution = makieSystemInteractive(molecules, domain, number_of_steps)

    # playing = Observable(false)
    # stop_animation = true

    # on(slider.value) do t
    #     pos[] = makieGetPositions(molecules, Int(t))
    #     z_distribution[1][] = [m.positions_history[Int(t)][3] for m in molecules if m.chemical_formula == chemical_formulas_distribution[1]]
    #     z_distribution[2][] = [m.positions_history[Int(t)][3] for m in molecules if m.chemical_formula == chemical_formulas_distribution[2]]
    #     z_values[], z_pressure_values[] = calcPressureZDistribution(molecules, Int(t), domain, 20)
    # end

    # on(button.clicks) do _
    #     playing[] = !playing[]
    # end

    # @async while stop_animation
    #     if playing[]
    #         slider.value[] = mod(slider.value[] , number_of_steps) + 1
    #     end
    #     sleep(1/FPS)
    # end
    

    # wait(display(fig))

    # stop_animation = false
end

main_lennardJones()
# main_lennardJones_condensation()
# main_entropy(true)
# main_temp_gradient()
# main_multi_species(false)
# main_standard(false)