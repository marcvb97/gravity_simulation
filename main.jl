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

##

# Number of particles
const N = 2000

particles = [gen_particle02(1.0,20.0,2.0) for _ in 1:N]
particles[1].pos = Point(0.0,0.0)
particles[1].vel = Point(0.0,0.0)
particles[1].mass = 100.00
set_zero_momentum!(particles)

c[1:h,1:w] .= 0.0
d = zeros(h,w)
particles_plot = particles

## Run the simulation
for step in 1:10000
    # println(step)
    # update position of particles
    # println("point 3 :",total_mass(particles))
    particles_new = new_position(particles,1.0)
    center_particle_1!(particles_new)
    # println("point 1 :",total_mass(particles_new))

    # Update visualization
    # plt[1] = [p.pos.x for p in particles]
    # plt[2] = [p.pos.y for p in particles]
    if mod(step,1) == 0
        for p in particles_plot
            size = round(Int,(p.mass*100)^(1/3))
            # plot_point(d,p.pos,0.0,-2.0,+2.0,-2.0,+2.0,size)
            plot_point(d,p.pos,0.0,-10.0,+10.0,-10.0,+10.0,size)
        end
        
        particles_plot = particles
        for p in particles
            size = round(Int,(p.mass*100)^(1/3))
            # plot_point(d,p.pos,1.0,-2.0,+2.0,-2.0,+2.0,size)
            plot_point(d,p.pos,1.0,-10.0,+10.0,-10.0,+10.0,size)
        end
        c[:,:] = d
        # sleep(0.00001)
    end
    particles = particles_new
    # println("point 2 :",total_mass(particles_new))
    # sleep(0.000001)
    if mod(step,10) == 0
        max_mass = 0.0
        for p in particles
            if p.mass > max_mass
                max_mass = p.mass
            end
        end
        println("step: $step, number of particles: $(length(particles)), heaviest mass: $max_mass")
    end
end

##

close(c)