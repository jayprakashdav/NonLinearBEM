using DrWatson
quickactivate("nonlinearmodeling")
using CompScienceMeshes, BEAST, LinearAlgebra, NonLinearBEM
r = 1e-8
h = 0.25e-8
include(srcdir("assets/genmesh.jl"))
Γ = meshsphere2(r, h)
X = raviartthomas(Γ)
X1 = BEAST.nedelec(Γ)
Δt0 = 0.07376
reduce_tstep = 1
Δt = Δt0/reduce_tstep
Nt = Int(floor(100*reduce_tstep))
T = timebasisshiftedlagrange(Δt, Nt, 3)
U = timebasisdelta(Δt, Nt)

V = X ⊗ T
W = X ⊗ U
V1 = X1 ⊗ T
W1 = X1 ⊗ U

duration = 40 * Δt0 * 2
delay = 1.0 * duration
amplitude = 3.5
gaussian = creategaussian(duration, delay, amplitude)
direction, polarisation = ẑ, x̂
E = planewave(polarisation, direction, derive(gaussian), 1.0)
chr = BEAST.Polynomial(0.6810,-1.68,1.2000)
#chr = BEAST.Polynomial(1.0,0.0,0.0)
SL = TDMaxwell3D.singlelayer(speedoflight=1.0, numdiffs=1)
idST = Identity()⊗Identity()
#σ = NLConductivity(chr,numdiffs=1)⊗IdentityTime()
σ = NonLinearBEM.tdconductivityb(chr,numdiffs=1)
#σ = BEAST.tdconductivityb(chr, numdiffs=1, init=true, init_j=xj, init_e=xe)
@hilbertspace j
@hilbertspace k
@hilbertspace l
@hilbertspace m

#=tdefie = @discretise(
         1.0SL[k,j] - 1.0id[k,l]  
        +1.0id[m,j] - 1.0 σ[m,l] == 0.0Identity()[m] + 1.0Ex[k], 
                                    j∈V, k∈W, l∈V1, m∈W1)=#
tdefie = @discretise 1.0SL[k,j] == -1.0E[k] j∈V k∈W
ohmslaw = @discretise 1.0idST[m,l] == 1.0σ[m] l∈V1 m∈W1
eq1 = tdefie
eq2 = ohmslaw
xj, xe, xj_all, xe_all = NonLinearBEM.tdsolve(eq1,eq2)

include(srcdir("save_solution.jl"))
save_sol([eq1, eq2], [xj, xe, xj_all, xe_all])