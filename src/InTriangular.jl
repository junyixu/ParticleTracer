module InTriangular
export is_inside
using OffsetArrays
mycross(a, b)= a[1]*b[2]- a[2]*b[1]
function sign_cross(x, a, b)
	mycross(a-x, a-b) |> sign
end

"""
 x is inside triangle (a,b,c) ?
"""
function is_inside(x, a::Vector, b::Vector, c::Vector)
	s1 = sign_cross(x, a, b)
	s2 = sign_cross(x, b, c)
	s3 = sign_cross(x, c, a)

    # 在三角形内部
	s1 == s2 == s3 && return true

    # 在三角形边界上
	s1 == 0.0 && s2 == s3 && return true
    s2 == 0.0 && s3 == s1 && return true
    s3 == 0.0 && s1 == s2 && return true

	return false 
end

is_inside(x,m::Matrix)=is_inside(x, m[:,1], m[:,2], m[:,3] )
end
