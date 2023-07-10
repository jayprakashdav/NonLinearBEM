mutable struct PortVolt{T,P,P1,P2,F1,F2} <: TDFunctional{T}
    mesh::P
    port1::P1
    port2::P2
    pv1::F1
    pv2::F2
    ids
    speedoflight::T
end

BEAST.scalartype(p::TDFunctional) = eltype(p.speedoflight)

function portvoltage(mesh,port1,port2,v1,v2,ids)
    PortVolt(mesh,port1,port2,v1,v2,ids,1.0)
end

struct Sinusoidal{T}
    amplitude::T
    angfrequency::T
    delay::T
end

(g::Sinusoidal)(s::Real) = g.amplitude*sin(g.angfrequency*s)

function createsinusoidal(amplitude, angfrequency)
    Sinusoidal(amplitude,angfrequency,0.0)
end