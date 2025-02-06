module Points

import Base: +, -, *
import LinearAlgebra: norm, ⋅

export +,-,*,norm, plot_point, Point, ⋅

struct Point{T}
    x::T 
    y::T 
end
+(p1::Point{T}, p2::Point{T}) where T = Point(p1.x + p2.x, p1.y + p2.y)
-(p1::Point{T}, p2::Point{T}) where T = Point(p1.x - p2.x, p1.y - p2.y)
-(p1::Point{T}) where T = Point(-p1.x, -p1.y)
*(p1::Point{T}, c) where T = Point(p1.x*c, p1.y*c) 
*(c,p1::Point{T}) where T = Point(c*p1.x, c*p1.y) 

function generate_point(x,y)
    return Point(x,y)
end

function norm(p::Point{T} where T)
    n = sqrt(p ⋅ p)
    return n
end

function ⋅(p::Point, q::Point)
    return p.x * q.x + p.y * q.y
end


function plot_point(cc,p,v,a,b,c,d,k)
    # plot point p in matrix cc
    # cc corresponds to the rectangle [a,b;c,d]
    # size of the plotted point is k
    m,n = size(cc)
    j = 1 + (p.x-a)*(m-1)/(b-a)
    i = n + (p.y-c)*(1-n)/(d-c)
    i = round(Int,i)
    j = round(Int,j)
    if (i-k >= 1) & (i+k <= n) & (j-k >= 1) & (j+k <=m)
        cc[i-k:i+k,j-k:j+k] .= v
    end
end

end