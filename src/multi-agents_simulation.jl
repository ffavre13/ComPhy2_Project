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
    lx::Float64
    ly::Float64
    lz::Float64
end

function domainVolume(cuboid::Domain)
    return cuboid.lx * cuboid.ly * cuboid.lz
end

function checkDomain(molecule::Molecule,domain::Domain)
    domain_pos = [domain.lx, domain.ly, domain.lz]

    for dim in eachindex(molecule.position)
        # left wall
        if molecule.position[dim] - molecule.radius < -domain_pos[dim]/2
            dist = -domain_pos[dim]/2 - (molecule.position[dim] - molecule.radius)
            
            molecule.position[dim] = molecule.position[dim] + 2*dist
            molecule.velocity[dim] = molecule.velocity[dim] * -1.0

        # right wall
        elseif molecule.position[dim] + molecule.radius > domain_pos[dim]/2
            dist = (molecule.position[dim] + molecule.radius) - domain_pos[dim]/2

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

function simulation(position::Vector{Vector{Float64}}, velocity::Vector{Vector{Float64}}, mass::Vector{Float64}, radius::Vector{Float64}, chemical_formula::Vector{String}, number_of_steps::Int64, dt::Float64, domain::Domain, g::Vector{Float64})
    molecules::Vector{Molecule} = Molecule[]

    for i in 1:length(position)
        push!(molecules,Molecule(position[i], velocity[i], mass[i], radius[i], chemical_formula[i], g, [zeros(Float64,3) for _ in 1:number_of_steps], [zeros(Float64,3) for _ in 1:number_of_steps]))
    end

    for m in molecules
        m.positions_history[1] .= m.position
        m.velocities_history[1] .= m.velocity
    end


    for t in 2:number_of_steps
        for m in molecules
            ComputeNextPosition(m, dt)
        end

        for m in molecules
            checkDomain(m, domain)
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

function plotSystem(molecules::Vector{Molecule}, t::Int64, domain::Domain)
    # camera=(0, 0)
    Plots.plot(legend=false,xlims=(-domain.lx/2,domain.lx/2),ylims=(-domain.ly/2,domain.ly/2),zlims=(-domain.lz/2,domain.lz/2))

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

    p = Plots.plot([1:length(molecules[1].velocities_history)],[calcEmec(molecules, t) for t in 1:length(molecules[1].velocities_history)], title="mechanical energy of the system", grid=false, legend=false, xlabel="time [s]", ylabel="Mechanical energy [J]")
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
        p = Plots.plot([1:length(molecules[1].velocities_history)],[calcQuantityOfMovement(molecules, t)[dim] for t in 1:length(molecules[1].velocities_history)], title="quantity of movement for axis $dim", grid=false, legend=false, xlabel="time [s]", ylabel="Quantity of movement (axis $dim) [kg*m/s]")
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
    p = Plots.plot([1:length(molecules[1].velocities_history)],[calcMeanVelocity(molecules, t) for t in 1:length(molecules[1].velocities_history)], title="mean velocity over time", grid=false, legend=false, xlabel="time [s]", ylabel="velocity [m/s]")
    display(p)
end

function plotVelocityDistributionFinal(molecules::Vector{Molecule})
    t_final = length(molecules[1].velocities_history)

    p = Plots.histogram([sqrt(sum(m.velocities_history[t_final] .^ 2)) for m in molecules], bins = 70, title="final velocity magnitude distribution", grid=false, legend=false, xlabel="velocity value [m/s]", ylabel="number of molecules")
    display(p)
end

function calcMeanVelocitySquare(molecules::Vector{Molecule}, t::Int64)
    velocities::Vector{Float64} = []

    for m in molecules
        push!(velocities, sqrt(sum(m.velocities_history[t] .^ 2)))
    end

    return mean(velocities .^ 2)
end

function calcTemperature(molecules::Vector{Molecule}, t::Int64)
    kb = 1.380649e-23
    mass = molecules[1].mass
    mean_velocity = calcMeanVelocitySquare(molecules, t)

    temperature = (mass * mean_velocity) /  (3*kb)

    return temperature
end

function plotTemperature(molecules::Vector{Molecule})
    p = Plots.plot([1:length(molecules[1].velocities_history)],[calcTemperature(molecules, t) for t in 1:length(molecules[1].velocities_history)], title="temperature over time", grid=false, legend=false, xlabel="time [s]", ylabel="Temperature [K]")
    display(p)
end

function calcPressure(molecules::Vector{Molecule}, t::Int64, cuboid::Domain)
    number_molecules = length(molecules) 
    mass = molecules[1].mass
    mean_velocity = calcMeanVelocitySquare(molecules, t)
    domain_volume = domainVolume(cuboid)

    pressure = (number_molecules * mass * mean_velocity) /  (3 * domain_volume)

    return pressure
end

function plotPressure(molecules::Vector{Molecule}, cuboid::Domain)
    p = Plots.plot([1:length(molecules[1].velocities_history)],[calcPressure(molecules, t, cuboid) for t in 1:length(molecules[1].velocities_history)], title="pressure over time", grid=false, legend=false, xlabel="time [s]", ylabel="Pressure [Pa]")
    display(p)
end


function plotPositionZDistribution(molecules::Vector{Molecule}, t::Int64)
    p = Plots.histogram([m.positions_history[t][3] for m in molecules], bins = 20, title="final position z distribution", grid=false, legend=false, normalize=:probability, xlabel="position z [m]", ylabel="percentage")
    display(p)
end

function calcPressureZDistribution(molecules::Vector{Molecule}, t::Int64, cuboid::Domain, number_bins::Int64)
    bin_size = cuboid.lz / number_bins

    molecules_in_bin::Vector{Vector{Molecule}} = [Molecule[] for _ in 1:number_bins]
    z_values = [-cuboid.lz/2 + (i-0.5)*bin_size for i in 1:number_bins]
    pressure_values = zeros(Float64, number_bins)

    for i in 1:number_bins
        for m in molecules
            if m.positions_history[t][3] >= -cuboid.lz/2 + bin_size*(i-1) && m.positions_history[t][3] < -cuboid.lz/2 + bin_size*i
                push!(molecules_in_bin[i], m)
            end
        end

        if length(molecules_in_bin[i]) > 0
            pressure_values[i] = calcPressure(molecules_in_bin[i], t, Domain(cuboid.lx, cuboid.ly, bin_size))
        end
    end
    
    return z_values, pressure_values
end

function plotPressureZDistribution(molecules::Vector{Molecule}, cuboid::Domain, t::Int64)
    z_values, pressure_values = calcPressureZDistribution(molecules, t, cuboid, 20)

    p = Plots.bar(z_values, pressure_values, title="pressure distribution along z axis", grid=false, legend=false, xlabel="position z [m]", ylabel="Pressure [Pa]")
    display(p)
end

function makieSystem(molecules, domain, number_of_steps)

    fig = Figure()
    
    ax = Axis3(fig[1,1],limits = (-domain.lx/2, domain.lx/2, -domain.ly/2, domain.ly/2, -domain.lz/2, domain.lz/2), xgridvisible = false, ygridvisible = false, zgridvisible = false)

    positions = Observable([Point3f(m.position[1], m.position[2], m.position[3]) for m in molecules])

    sizes = [m.radius for m in molecules]
    meshscatter!(ax, positions, markersize = sizes)

    return fig, positions
end

function makieSystemInteractive(molecules, domain, number_of_steps)

    fig = Figure()
    
    ax = Axis3(fig[1,1],limits = (-domain.lx/2, domain.lx/2, -domain.ly/2, domain.ly/2, -domain.lz/2, domain.lz/2), xgridvisible = false, ygridvisible = false, zgridvisible = false)
    ax_distribution_z = Axis(fig[3,1], title="Position z distribution", xlabel="position z [m]", ylabel="percentage", xgridvisible = false, ygridvisible = false)
    ax_pressure_z = Axis(fig[4,1], title="Pressure distribution along z axis", xlabel="position z [m]", ylabel="Pressure [Pa]", xgridvisible = false, ygridvisible = false)

    slider = Slider(fig[2,1], range = 1:number_of_steps, startvalue = 1)
    button = Button(fig[2,2], label = "Play")
    
    positions = Observable([Point3f(m.position[1], m.position[2], m.position[3]) for m in molecules])
    sizes = [m.radius for m in molecules]
    meshscatter!(ax, positions, markersize = sizes)

    z_distribution = Observable([m.positions_history[1][3] for m in molecules])

    hist!(ax_distribution_z, z_distribution, bins = 20, normalization = :probability)

    z_values = Observable(Float64[])
    z_pressure_values = Observable(Float64[])

    z_init, p_init = calcPressureZDistribution(molecules, 1, domain, 20)

    z_values[] = z_init
    z_pressure_values[] = p_init

    barplot!(ax_pressure_z, z_values, z_pressure_values)

    return fig, positions, slider, button, z_distribution, z_values, z_pressure_values
end

function makieGetPositions(molecules, t)
    return [Point3f(m.positions_history[t][1],m.positions_history[t][2],m.positions_history[t][3]) for m in molecules]
end

function main()
    number_of_steps::Int64 = 10000
    FPS = 60

    g = [0,0,-9.81 * 10^13]

    dt::Float64 = 1.0e-14

    positions::Vector{Vector{Float64}} = []
    velocities::Vector{Vector{Float64}} = []
    masses::Vector{Float64} = []
    radius::Vector{Float64} = []
    chemical_formulas::Vector{String} = []

    domain::Domain = Domain(4e-8,4e-8,4e-8)
    velocityValue::Float64 = 1367 # [m/s]
    
    number_atomes = 500

    # Random generation

    for i in 1:number_atomes
        push!(positions, [rand()*domain.lx-domain.lx/2,rand()*domain.ly-domain.ly/2,rand()*domain.lz-domain.lz/2])

        velocity::Vector{Float64} = [rand()*10-5,rand()*10-5,rand()*10-5]
        velocity = velocity ./ sqrt(sum(velocity .^2))
        velocity = velocity .* velocityValue

        @assert round(Int, sqrt(sum(velocity .^2))) == velocityValue

        push!(velocities, velocity)
        push!(masses, 6.646e-27)
        push!(radius, 1.1e-10)
        push!(chemical_formulas, "He")
    end



    @assert length(positions) == length(velocities) == length(masses) == length(radius) == length(chemical_formulas)

    molecules::Vector{Molecule} = simulation(positions,velocities,masses,radius,chemical_formulas, number_of_steps, dt, domain, g)


    # With Plots

    plotEmec(molecules)
    plotQuantityOfMovement(molecules)
    plotMeanVelocity(molecules)
    plotVelocityDistributionFinal(molecules)
    plotTemperature(molecules)
    plotPressure(molecules, domain)
    plotPositionZDistribution(molecules, length(molecules[1].positions_history))
    plotPressureZDistribution(molecules, domain, length(molecules[1].positions_history))

    # filename = "results/molecule.mp4"

    # animation = @animate for t in 1:number_of_steps
    #     plotSystem(molecules, t, domain)
    # end

    # mp4(animation, filename, fps = FPS)


    # With makie

    fig, pos  = makieSystem(molecules, domain, number_of_steps)

    record(fig, "results/molecule.mp4", 1:number_of_steps) do t
        pos[] = makieGetPositions(molecules, t)
    end

    fig, pos, slider, button, z_distribution, z_values, z_pressure_values = makieSystemInteractive(molecules, domain, number_of_steps)

    playing = Observable(false)
    stop_animation = true

    on(slider.value) do t
        pos[] = makieGetPositions(molecules, Int(t))
        z_distribution[] = [m.positions_history[Int(t)][3] for m in molecules]
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

main()