using SimpleCanvas

w,h = 1080,1080
c = canvas(zeros(h,w))

for j = 1:size(c,2)
    i = (sin(j/100.0)+1.1)*350;
    i = round(Int,i);
    c[i-5:i+5,j] .= 1.0;
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

name!(c,"Test SimpleCanvas")