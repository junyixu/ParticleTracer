module Constants

# R0 = 1995.0
# x0 = [2112.0, 0, 0]
# p0 = [5.0, 1, 0]

import PhysicalConstants.CODATA2018 as C

export Unit

constants=(:e, :m_p, :m_e)
for constant in constants
    eval(:($constant=C.$constant.val))
end

c = C.c_0.val
m_α = 4m_p


struct Unit
    B::Float64
    E::Float64 # 电场
    Ω::Float64
    t::Float64
    x::Float64 # 位置
    p::Float64
    v::Float64
    ε::Float64 # 能量
    q::Float64
    m::Float64
    function Unit(B::Float64, q::Float64, m::Float64)
        E = B*c;
        Ω=q*B/m;
        t= 1/Ω;
        x = m*c/(q*B);
        p=m*c;
        v=c;
        ε=m*c^2;
        new(B, E, Ω, t, x, p, v, ε, q, m)
    end
end

function Unit(B::Float64, type::String)
	if type == "alpha"
		return Unit(B, 2*e, m_α) # 有时候需要改单位电量和单位质量
	end
	return Unit(B, e, me)
end
Unit(B::Float64)=Unit(B, e, m_e)


end
