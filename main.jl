## run at julia start; to be ran once.
# Tested on Julia 1.11.3
working_directory = @__DIR__
cd(working_directory)
if working_directory ∉ LOAD_PATH
    push!(LOAD_PATH, working_directory) 
end

using Pkg
Pkg.activate("./Env")  

using Revise

## Main code
using Points
using Particles

using SimpleCanvas

w,h = 1080,1080
c = canvas(zeros(h,w))
name!(c,"Gravity Simulation")

function cmap(value::Float64)
    if value < 0.0; r = g = b = 0.0;
    elseif value <= 0.25; r = 0.0; g = 4 * value; b = 1.0; 
    elseif value <= 0.5; r = 0.0; g = 1.0; b = 1.0 - 4 * (value - 0.25); 
    elseif value <= 0.75; r = 4 * (value - 0.5); g = 1.0; b = 0.0
    elseif value <= 1.0; r = 1.0; g = 1.0 - 4 * (value - 0.75); b = 0.0
    else; r = g = b = 1.0; end
    round.(UInt8, 255 .*(r,g,b))
end
colormap!(c, cmap) # c::Canvas

##

# Number of particles
const N = 2000

particles = [gen_particle02(1.0,20.0,2.0) for _ in 1:N]
particles[1].pos = Point(0.0,0.0)
particles[1].vel = Point(0.0,0.0)
particles[1].mass = 100.00
set_zero_momentum!(particles)

d = zeros(h,w).-1.0
particles_plot = particles



## Run the simulation
zoomlvl = 10
d.=-1.0
for step in 1:100000
    # println(step)
    # update position of particles
    # println("point 3 :",total_mass(particles))
    particles_new = new_position(particles,1.0)
    # center_particle_1!(particles_new)
    # println("point 1 :",total_mass(particles_new))

    # Update visualization
    # plt[1] = [p.pos.x for p in particles]
    # plt[2] = [p.pos.y for p in particles]
    if mod(step,1) == 0
        for p in particles_plot
            size = round(Int,(p.mass*100)^(1/3))
            # plot_point(d,p.pos,0.0,-2.0,+2.0,-2.0,+2.0,size)
            plot_point(d,p.pos,-1.0,-zoomlvl,+zoomlvl,-zoomlvl,+zoomlvl,size)
        end
        
        particles_plot = particles
        for p in particles
            size = round(Int,(p.mass*100)^(1/3))
            # plot_point(d,p.pos,1.0,-2.0,+2.0,-2.0,+2.0,size)
            # color = norm(p.vel)/100
            color = p.E/p.mass * (1/100) 
            plot_point(d,p.pos,color,-zoomlvl,+zoomlvl,-zoomlvl,+zoomlvl,size)
        end
        c[:,:] = d
        # sleep(0.00001)
    end
    particles = particles_new
    # println("point 2 :",total_mass(particles_new))
    # sleep(0.000001)
    if mod(step,10) == 0
        max_E = 0.0
        max_mass = 0.0
        max_E_per_mass = 0.0
        for p in particles
            if p.E > max_E; max_E = p.E; end
            if p.mass > max_mass; max_mass = p.mass; end
            if p.E/p.mass > max_E_per_mass; max_E_per_mass = p.E/p.mass; end
        end
        println("step: $step, number of particles: $(length(particles)), largest E/mass: $max_E_per_mass")
    end
end

##

close(c)