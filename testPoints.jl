# test the module Points.jl

include("Points.jl")
using .Points

using SimpleCanvas

p = Point(1.0,2.0)

norm(p)

q = p + 2*p

w,h = 1080,1080
cc = canvas(zeros(h,w))

for x = 0.0:0.1:7.0
    y = sin(x)
    p = Point(x,y)
    v = 1.0     # white
    a = 0.0
    b = 7.0
    c = -1.2
    d = +1.2
    k = 2       # size of plotted point
    plot_point(cc,p,v,a,b,c,d,k)
end