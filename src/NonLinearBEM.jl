module NonLinearBEM

using BEAST, CompScienceMeshes, LinearAlgebra
using StaticArrays

import BEAST.LocalOperator, BEAST.Functional
import BEAST.SpaceTimeBasis
import BEAST.Polynomial
import BEAST.LinearRefSpaceTriangle
import BEAST.SingleNumQStrat
import BEAST.RefSpace
import BEAST.DiscreteEquation
import BEAST.TDFunctional
import BEAST.RTBasis
import BEAST.td_assemble

export IdentityTime
export NLConductivity
export ConductivityTD
export raviartthomaswithport

include("maxwell/timedomain/finite_conductivity.jl")
include("maxwell/timedomain/conductivityTD.jl")
include("maxwell/timedomain/finite_conductivityNL.jl")
include("bases/port_basis.jl")
include("bases/portvoltage.jl")
include("bases/test_gmsh.jl")

end
