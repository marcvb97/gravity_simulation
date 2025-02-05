using SimpleCanvas

w,h = 1080,1080
c = canvas(zeros(h,w))
m = zeros(h,w);
c = canvas(m); # c can be used as if it were a Matrix
close(c)

for i = 1:size(m,1)
    j = (sin(i/100.0)+1.1)*350;
    j = ceil(Int,j);
    c[i,j-5:j+5] .= 1.0;
end
# Drawing random rectangles
function draw_rectangles(n, rectangle_size)
    rs = rectangle_size-1
    for i = 1:n
        i%100 == 0 && println("drawing $i")
        x,y = rand(1:h-rs), rand(1:w-rs)
        c[x:x+rs, y:y+rs].= 1 .-c[x:x+rs, y:y+rs]
    end
end
draw_rectangles(10000, 50)

function draw_rand(n)
    for i = 1:n
        println(i)
        c.=rand(size(c)...)
    end
end
draw_rand(100)

function cmap(value::Float64)
    if value < 0.0; r = g = b = 0.0;
    elseif value <= 0.25; r = 0.0; g = 4 * value; b = 1.0; 
    elseif value <= 0.5; r = 0.0; g = 1.0; b = 1.0 - 4 * (value - 0.25); 
    elseif value <= 0.75; r = 4 * (value - 0.5); g = 1.0; b = 0.0
    elseif value <= 1.0; r = 1.0; g = 1.0 - 4 * (value - 0.75); b = 0.0
    else; r = g = b = 0.0; end
    round.(UInt8, 255 .*(r,g,b))
end
colormap!(c, cmap) # c::Canvas



using Revise
# using GLMakie

import Base: +, -, *

# Number of particles
const N = 100


struct Point{T}
    x::T 
    y::T 
end
+(p1::Point{T}, p2::Point{T}) where T = Point(p1.x + p2.x, p1.y + p2.y)
-(p1::Point{T}, p2::Point{T}) where T = Point(p1.x - p2.x, p1.y - p2.y)
-(p1::Point{T}) where T = Point(-p1.x, -p1.y)
*(p1::Point{T}, c) where T = Point(p1.x*c, p1.y*c) 
*(c,p1::Point{T}) where T = Point(c*p1.x, c*p1.y) 

import LinearAlgebra: norm
function norm(p::Point{T} where T)
    n = sqrt(p.x*p.x+p.y*p.y)
    return n
end

# Create particle data in Julia
mutable struct Particle
    pos::Point{Float64}
    vel::Point{Float64}
    mass::Float64
    temp::Float64
end

# Initialize particles with random positions
function gen_particle()
    speed = 0.0
    T = Float64
    pos = Point(rand(T) * 2 - 1, rand(T) * 2 - 1) * 0.50
    vel = Point(rand(T) * 2 - 1, rand(T) * 2 - 1) * speed
    temp = 1.0;
    mass = 0.01;
    return Particle(pos,vel,mass,temp)
end

# Initialize particles in spiralform
function gen_particle02()
    speed = 2.0
    T = Float64
    pos = Point(rand(T) * 2 - 1, rand(T) * 2 - 1) * 0.1
    vel = Point(-pos.y,pos.x)
    vel = vel * (1/norm(vel)^1.0) * speed
    temp = 1.0;
    mass = 0.01;
    return Particle(pos,vel,mass,temp)
end

particles = [gen_particle02() for _ in 1:N]

# Assumes all particles have the same mass
function set_zero_momentum!(particles::Array{Particle})
    momentum = sum(p.vel*p.mass for p in particles)
    mass = sum(p.mass for p in particles)
    momentum_per_mass = momentum * (1/mass)
    for p in particles 
        p.vel -= momentum_per_mass * p.mass
    end
end
set_zero_momentum!(particles)

# Visualization using GLMakie
# fig, ax, plt = scatter([p.pos.x for p in particles], [p.pos.y for p in particles], color=:blue, markersize=10)
# display(fig)
# xlims!(-1.5,+1.5)
# ylims!(-1.5,+1.5)
import Base.copy 
copy(p::Particle) = Particle(p.pos, p.vel, p.mass, p.temp)
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

function collision_particles(particles)
    r = 0.001
    t_mass = total_mass(particles)
    println("begin collision_particles ",t_mass)
    # sleep(0.5)
    i = 1
    while i <= length(particles)
        p = particles[i]
        j = i+1
        while j <= length(particles)
            q = particles[j]
            diff = q.pos-p.pos
            dist = sqrt(diff.x*diff.x+diff.y*diff.y)
            if dist < r
                println("two particles collapse into one",i,j)
                p.mass = p.mass+q.mass
                p.vel = (p.vel*p.mass+q.vel*q.mass) * (1 / (p.mass+q.mass))
                p.pos = (p.pos*p.mass+q.pos*q.mass) * (1 / (p.mass+q.mass))
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
        particles[i] = p
        i = i+1
        check_mass(particles,t_mass)
    end
    t_mass = total_mass(particles)
    println("end collision_particles ",t_mass)
    return particles
end

function new_position(particles,h)
    particles = collision_particles(particles)
    dt = 0.0001
    println("point 4 :",total_mass(particles))
    particles_new = copy.(particles)
    println("point 5 :",total_mass(particles_new))
    i = 0
    for j = 1:length(particles)
        p = particles[j]
        k1 = compute_force(p,particles,j)
        p_new = copy(p)
        p_new.vel += k1 * dt * h * 0.5
        p_new.pos += p.vel*dt * h * 0.5
        k2 = compute_force(p_new,particles,j)
        p_new = copy(p)
        p_new.vel += k2 * dt * h * 0.5
        p_new.pos += p.vel*dt * h * 0.5
        k3 = compute_force(p_new,particles,j)
        p_new = copy(p)
        p_new.vel += k3 * dt * h 
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
        if abs(particles_new[i].pos.x) > +1
            particles_new[i].vel = Point(-particles_new[i].vel.x,particles_new[i].vel.y)
        end
        if abs(particles_new[i].pos.y) > +1
            particles_new[i].vel = Point(particles_new[i].vel.x,-particles_new[i].pos.y)
        end
    end
    return particles_new
end

function plot_point(cc,p,v,a,b,c,d,k)
    m,n = size(cc)
    j = 1 + (p.x-a)*(m-1)/(b-a)
    i = n + (p.y-c)*(1-n)/(d-c)
    i = round(Int,i)
    j = round(Int,j)
    cc[i-k:i+k,j-k:j+k] .= v

end


c[1:h,1:w] .= 0.0
# Run the simulation
particles_plot = particles
for step in 1:100000
    println(step)
    # update position of particles
    println("point 3 :",total_mass(particles))
    particles_new = new_position(particles,1.0)
    println("point 1 :",total_mass(particles_new))

    # Update visualization
    # plt[1] = [p.pos.x for p in particles]
    # plt[2] = [p.pos.y for p in particles]
    if mod(step,10) == 0
        for p in particles_plot
            plot_point(c,p.pos,0.0,-1.25,+1.25,-1.25,+1.25,round(Int,sqrt(p.mass*100)))
        end
        
        particles_plot = particles
        for p in particles
            plot_point(c,p.pos,1.0,-1.25,+1.25,-1.25,+1.25,round(Int,sqrt(p.mass*100)))
        end
    end
    particles = particles_new
    println("point 2 :",total_mass(particles_new))
    # sleep(0.000001)
end

close(c)