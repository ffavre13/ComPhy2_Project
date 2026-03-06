using Plots
using GLMakie

"""
    Molecule

Class for storing information about a molecule.

### Fields
- `position::Float64`: The molecule's position [m].
- `velocity::Float64`: The molecule's velocity [m/s].

- `mass::Float64`: The mass of the molecule [kg].
- `radius::Float64`: The radius of the molecule [m].
- `chimical_formula::String`: The chemical formula of the molecule.

- `positions_history::Vector{Vector{Float64}}`: A vector storing the history of the molecule's position.
- `velocities_history::Vector{Vector{Float64}}`: A vector storing the history of the molecule's velocity.
"""
mutable struct Molecule
    position::Vector{Float64}
    velocity::Vector{Float64}

    mass::Float64
    radius::Float64
    chimical_formula::String

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
        if molecule.position[dim] - molecule.radius <= -domain_pos[dim]/2
            dist = abs(-domain_pos[dim]/2 - (molecule.position[dim] - molecule.radius))
            
            molecule.position[dim] = -domain_pos[dim]/2 + dist
            molecule.velocity[dim] *= -1

        elseif molecule.position[dim] + molecule.radius >= domain_pos[dim]/2
            dist = abs((molecule.position[dim] + molecule.radius) - domain_pos[dim]/2)

            molecule.position[dim] = domain_pos[dim]/2 - dist
            molecule.velocity[dim] *= -1
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
    molecule.position = molecule.position + dt .* molecule.velocity
end

function simulation(position::Vector{Vector{Float64}}, velocity::Vector{Vector{Float64}}, mass::Vector{Float64}, radius::Vector{Float64}, chimical_formula::Vector{String}, number_of_steps::Int64, dt::Float64, domain::Domain)
    molecules::Vector{Molecule} = Molecule[]

    for i in 1:length(position)
        push!(molecules,Molecule(position[i], velocity[i], mass[i], radius[i], chimical_formula[i], [zeros(Float64,3) for _ in 1:number_of_steps], [zeros(Float64,3) for _ in 1:number_of_steps]))
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
        emec += 1/2 * m.mass * (sum(m.velocity .^ 2))
    end

    return emec
end

function plotEmec(molecules::Vector{Molecule})
    val = calcEmec(molecules, 1)

    p = Plots.plot([1:length(molecules[1].velocities_history)],[calcEmec(molecules, t) for t in 1:length(molecules[1].velocities_history)], title="mechanical energy of the system", grid=false, legend=false, xlabel="time [s]", ylabel="Mechanical energy [J]",ylims=(val-0.001,val+0.001))
    display(p)
end

function calcQuantityOfMovement(molecules::Vector{Molecule}, t::Int64)
    p::Vector{Float64} = zeros(Float64,length(molecules[1].velocity))

    for m in molecules
        p .+= m.mass .* m.velocity
    end

    return p
end

function PlotQuantityOfMovement(molecules::Vector{Molecule})
    val = calcQuantityOfMovement(molecules, 1)

    for dim in 1:length(molecules[1].velocity)
        p = Plots.plot([1:length(molecules[1].velocities_history)],[calcQuantityOfMovement(molecules, t)[dim] for t in 1:length(molecules[1].velocities_history)], title="quantity of movement for axis $dim", grid=false, legend=false, xlabel="time [s]", ylabel="Quantity of movement (axis $dim) [kg*m/s]", ylims=(val[dim]-0.0001,val[dim]+0.0001))
        display(p)
    end
end

function makieSystem(molecules, domain)

    fig = Figure()

    ax = Axis3(fig[1,1],limits = (-domain.lx/2, domain.lx/2, -domain.ly/2, domain.ly/2, -domain.lz/2, domain.lz/2))

    positions = Observable([Point3f(m.position[1], m.position[2], m.position[3]) for m in molecules])

    meshscatter!(ax, positions, markersize = 0.02)

    return fig, positions
end

function makieGetPositions(molecules, t)
    return [Point3f(m.positions_history[t][1],m.positions_history[t][2],m.positions_history[t][3]) for m in molecules]
end

function main()
    number_of_steps::Int64 = 200
    FPS = 30

    dt::Float64 = 0.001

    positions::Vector{Vector{Float64}} = []
    velocities::Vector{Vector{Float64}} = []
    masses::Vector{Float64} = []
    radius::Vector{Float64} = []
    chimical_formulas::Vector{String} = []

    domain::Domain = Domain(1.0,1.0,1.0)
    
    for i in 1:10
        push!(positions, [rand()*domain.lx-domain.lx/2,rand()*domain.ly-domain.ly/2,rand()*domain.lz-domain.lz/2])
        push!(velocities, [rand()*60.0-30.0,rand()*60.0-30.0,rand()*60.0-30.0])
        push!(masses, rand()*5.0)
        push!(radius, 0.01*rand())
        push!(chimical_formulas, "TEST")
    end



    @assert length(positions) == length(velocities) == length(masses) == length(radius) == length(chimical_formulas)

    molecules::Vector{Molecule} = simulation(positions,velocities,masses,radius,chimical_formulas, number_of_steps, dt, domain)

    # With makie

    fig, pos = makieSystem(molecules, domain)

    display(fig)

    record(fig, "results/molecule.mp4", 1:number_of_steps) do t
        pos[] = makieGetPositions(molecules, t)
    end

    # With Plots

    # plotEmec(molecules)
    # PlotQuantityOfMovement(molecules)

    # filename = "results/molecule.mp4"

    # animation = @animate for t in 1:number_of_steps
    #     plotSystem(molecules, t, domain)
    # end

    # mp4(animation, filename, fps = FPS)
end

main()