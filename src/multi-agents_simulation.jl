using Plots
using GLMakie
using Statistics

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

mutable struct Domain
    lx::Tuple{Float64, Float64}
    ly::Tuple{Float64, Float64}
    lz::Tuple{Float64, Float64}
end

function domainVolume(domain::Domain)
    return abs(domain.lx[2] - domain.lx[1]) * abs(domain.ly[2] - domain.ly[1]) * abs(domain.lz[2] - domain.lz[1])
end

function checkDomain(molecule::Molecule,domain::Domain)
    domain_pos = [domain.lx, domain.ly, domain.lz]

    for dim in eachindex(molecule.position)
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

function simulation(position::Vector{Vector{Float64}}, velocity::Vector{Vector{Float64}}, mass::Vector{Float64}, radius::Vector{Float64}, chemical_formula::Vector{String}, number_of_steps::Int64, dt::Float64, domain::Domain, g::Vector{Float64}, remove_wall::Bool, new_domain::Domain)
    molecules::Vector{Molecule} = Molecule[]    
    current_domain = domain
    
    for i in 1:length(position)
        push!(molecules,Molecule(position[i], velocity[i], mass[i], radius[i], chemical_formula[i], g, [zeros(Float64,3) for _ in 1:number_of_steps], [zeros(Float64,3) for _ in 1:number_of_steps]))
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
                    checkDomain(m, new_domain)
                else 
                    checkDomain(m, current_domain)
                end
            else
                checkDomain(m, current_domain)
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
    p = Plots.histogram(val, bins = 70, title="final velocity magnitude distribution", grid=false, legend=false, xlabel="velocity value [m/s]", ylabel="number of molecules", ylims=:auto, xlims=:auto)
    display(p)
end

function calcMeanVelocitySquare(molecules::Vector{Molecule}, t::Int64)
    velocities::Vector{Float64} = []

    for m in molecules
        push!(velocities, m.mass * sum(m.velocities_history[t] .^ 2))
    end

    return mean(velocities)
end

function calcTemperature(molecules::Vector{Molecule}, t::Int64)
    kb = 1.380649e-23
    # mass = molecules[1].mass
    mean_velocity = calcMeanVelocitySquare(molecules, t)

    # temperature = (mass * mean_velocity) /  (3*kb)
    temperature = (mean_velocity) /  (3*kb)

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

function main(remove_wall::Bool = true)
    number_of_steps::Int64 = 10000
    FPS = 240

    g = [0.0,0.0,0.0]

    dt::Float64 = 1.0e-14

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

    spawn_domain::Domain = Domain((-5e-9,0.0), (-5e-9, 0.0), (-5e-9, 0.0))

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

    molecules::Vector{Molecule} = simulation(positions,velocities,masses,radius,chemical_formulas, number_of_steps, dt, domain, g, remove_wall, new_domain)


    # Verification of the stability of the simulation

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

    # With Plots

    if remove_wall
        plotEntropieRemoveWall(molecules, domain, new_domain, number_of_steps)
    else
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

        # # Save the video
        # fig, pos  = makieSystem(molecules, domain, number_of_steps)

        # record(fig, "results/molecule.mp4", 1:number_of_steps) do t
        #     pos[] = makieGetPositions(molecules, t)
        # end

        fig, pos, slider, button, z_distribution, z_values, z_pressure_values, chemical_formulas_distribution = makieSystemInteractive(molecules, domain, number_of_steps)

        playing = Observable(false)
        stop_animation = true

        on(slider.value) do t
            pos[] = makieGetPositions(molecules, Int(t))
            z_distribution[1][] = [m.positions_history[Int(t)][3] for m in molecules if m.chemical_formula == chemical_formulas_distribution[1]]
            z_distribution[2][] = [m.positions_history[Int(t)][3] for m in molecules if m.chemical_formula == chemical_formulas_distribution[2]]
            z_values[], z_pressure_values[] = calcPressureZDistribution(molecules, Int(t), domain, 20)
        end

        on(button.clicks) do _
            playing[] = !playing[]
        end

        @async while stop_animation
            if playing[]
                slider.value[] = mod(slider.value[] , number_of_steps) + 1
            end
            sleep(1/FPS)
        end
        

        wait(display(fig))

        stop_animation = false
    end
end

main(false)