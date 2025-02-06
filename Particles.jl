module Particles

using Points

import Base.copy 

export Particle
export gen_particle, gen_particle02, set_zero_momentum!
export compute_force, check_mass, total_mass
export collision_particles, new_position, center_particle_1!

# Create particle data in Julia
mutable struct Particle
    pos::Point{Float64}
    vel::Point{Float64}
    mass::Float64
    E::Float64
end

# Initialize particles with random positions
function gen_particle()
    speed = 0.0
    T = Float64
    pos = Point(rand(T) * 2 - 1, rand(T) * 2 - 1) * 0.50
    vel = Point(rand(T) * 2 - 1, rand(T) * 2 - 1) * speed
    E = 0.0;
    mass = 0.01;
    return Particle(pos,vel,mass,E)
end

# Initialize particles in spiralform
function gen_particle02(dist,speed,e)
    # speed = 1.5
    T = Float64
    pos = Point(rand(T) * 2 - 1, rand(T) * 2 - 1) * dist
    vel = Point(-pos.y,pos.x)
    vel = vel * (1/norm(vel)^e) * speed
    E = 0.0;
    mass = 0.01;
    if norm(pos)<1.0
        pos *= 1/norm(pos)
    end
    return Particle(pos,vel,mass,E)
end


# Assumes all particles have the same mass
function set_zero_momentum!(particles::Array{Particle})
    momentum = sum(p.vel*p.mass for p in particles)
    mass = sum(p.mass for p in particles)
    for p in particles 
        p.vel -= momentum * (1/mass)
    end
end



copy(p::Particle) = Particle(p.pos, p.vel, p.mass, p.E)

function compute_force(p, particles,j)
    eps = 0.01
    G = 10.0
    force = Point(0.0,0.0)
    for i = 1:length(particles)
        q = particles[i]
        if (i == j) 
            continue
        end
        diff = q.pos - p.pos
        dist2 = diff.x*diff.x+diff.y*diff.y+eps
        invDist = 1 / sqrt(dist2)
        force += diff * G * p.mass * q.mass * invDist * invDist * invDist
    end
    return force
end    

function check_mass(particles,total_mass)
    s = 0.0
    for j = 1:length(particles)
        s = s+particles[j].mass
    end
    if abs(s-total_mass) > 0.000001
        println(particles)
        sleep(10)
    end
end

function total_mass(particles)
    t = 0.0
    for j = 1:length(particles)
        t = t+particles[j].mass
    end
    return t
end

function collide(p::Particle, q::Particle)
    mass = p.mass+q.mass
    vel = (p.vel*p.mass+q.vel*q.mass) * (1 / (p.mass+q.mass))
    pos = (p.pos*p.mass+q.pos*q.mass) * (1 / (p.mass+q.mass))
    Estart = (p.vel ⋅ p.vel * 0.5p.mass) + (q.vel ⋅ q.vel * 0.5q.mass)
    Eend =  (vel ⋅ vel * 0.5mass)
    ΔE = Estart - Eend
    if ΔE < 0
        println("ending energy: $Eend > $Estart, starting energy")
    end
    E = p.E + q.E + ΔE
    return Particle(pos,vel,mass,E)
end

function collision_particles(particles)
    r_scale = 0.01
    t_mass = total_mass(particles)
    # println("begin collision_particles ",t_mass)
    # sleep(0.5)
    i = 1
    while i <= length(particles)
        p = particles[i]
        j = i+1
        while j <= length(particles)
            q = particles[j]
            diff = q.pos-p.pos
            dist = sqrt(diff.x*diff.x+diff.y*diff.y)
            if dist < r_scale * (cbrt(p.mass) + cbrt(q.mass))
                # println("two particles collapse into one",i,j)
                r = collide(p,q)
                particles[i] = r
                p = r
                if j < length(particles)
                    particles[j] = particles[length(particles)]
                else
                    j = j+1
                end
                particles = particles[1:length(particles)-1]
            else
                j = j+1
            end
        end
        
        i = i+1
        check_mass(particles,t_mass)
    end
    t_mass = total_mass(particles)
    # println("end collision_particles ",t_mass)
    return particles
end

function new_position(particles,h)
    particles = collision_particles(particles)
    dt = 0.0001
    # println("point 4 :",total_mass(particles))
    particles_new = copy.(particles)
    # println("point 5 :",total_mass(particles_new))
    i = 0
    for j = 1:length(particles)
        p = particles[j]
        k1 = compute_force(p,particles,j)
        p_new = copy(p)
        p_new.vel += k1 * dt * h * 0.5 * (1/p.mass)
        p_new.pos += p.vel*dt * h * 0.5
        k2 = compute_force(p_new,particles,j)
        p_new = copy(p)
        p_new.vel += k2 * dt * h * 0.5 * (1/p.mass)
        p_new.pos += p.vel*dt * h * 0.5
        k3 = compute_force(p_new,particles,j)
        p_new = copy(p)
        p_new.vel += k3 * dt * h * (1/p.mass)
        p_new.pos += p.vel*dt * h
        k4 = compute_force(p_new,particles,j)
        k_average = (k1+2*k2+2*k3+k4)*(1/6.0)
        d1 = k2 - k1
        d2 = k3 - k2
        d3 = k4 - k3
        # println(k1)
        # println(d1)
        # println(d2)
        # println(d3)
        # sleep(1)
        i = i+1
        particles_new[i].vel += k1 * dt * h * (1/p.mass)
        particles_new[i].pos += p.vel*dt * h 
        # if abs(particles_new[i].pos.x) > +1
        #     particles_new[i].vel = Point(-particles_new[i].vel.x,particles_new[i].vel.y)
        # end
        # if abs(particles_new[i].pos.y) > +1
        #     particles_new[i].vel = Point(particles_new[i].vel.x,-particles_new[i].pos.y)
        # end
    end
    return particles_new
end


function center_particle_1!(particles)
    pos1 = particles[1].pos
    for i = 2:length(particles)
        particles[i].pos = particles[i].pos - pos1
    end
    particles[1].pos = particles[1].pos - pos1
end

end