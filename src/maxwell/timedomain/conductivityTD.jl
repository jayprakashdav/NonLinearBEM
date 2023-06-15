abstract type ConductivityTD_model <: BEAST.Functional end

abstract type ConductivityTDop_model <: BEAST.LocalOperator end

mutable struct ConductivityTD <: ConductivityTD_model
    chr::BEAST.Polynomial
    dchr::BEAST.Polynomial
    numdiffs::Int
    efield::Array{SVector{3, Float64},2}
    jflux::Array{SVector{3, Float64},2}
end

#= mutable struct ConductivityTD{L, L2, M, N} <: ConductivityTD_model
    chr::L
    dchr::L2
    numdiffs::M
    efield::N
    jflux::N
end =#

mutable struct ConductivityTD_b <: ConductivityTD_model
    chr::BEAST.Polynomial
    dchr::BEAST.Polynomial
    numdiffs::Int
    efield::Array{SVector{3, Float64},2}
    jflux::Array{SVector{3, Float64},2}
    initialcondition::Bool
    initialxj::Array{Float64, 2}
    initialxe::Array{Float64, 2}
end

#= mutable struct ConductivityTD_b{L, L2, M, N, O, P} <: ConductivityTD_model
    chr::L
    dchr::L2
    numdiffs::M
    efield::N
    jflux::N
    initialcondition::O
    initialxj::P
    initialxe::P
end =#

mutable struct ConductivityTDop <: ConductivityTDop_model
    op::ConductivityTD
end

mutable struct ConductivityTDop_b <: ConductivityTDop_model
    op::ConductivityTD_b
end

BEAST.scalartype(p::NonLinearBEM.ConductivityTD_model) = eltype(p.efield)
BEAST.scalartype(op::NonLinearBEM.ConductivityTDop_model) = Float64

function tdconductivity(chr::BEAST.Polynomial;numdiffs=0)
    dchr = derive(chr)
    efield = fill(SVector(0.0,0.0,0.0), (2,2))
    ConductivityTD(chr, dchr, numdiffs, efield, efield)
end

function tdconductivityb(chr::BEAST.Polynomial;numdiffs=0, init=false, init_j=fill(0.0,(2,2)), init_e=fill(0.0,(2,2)))
    dchr = derive(chr)
    efield = fill(SVector(0.0,0.0,0.0), (2,2))
    ConductivityTD_b(chr, dchr, numdiffs, efield, efield, init, init_j, init_e)
end

function (f::ConductivityTD)(cell, cqdpt)
    ei = f.efield[cell, cqdpt]
    if norm(ei)==0
        ei = 1e-9.+ei
    end
    dsigma = f.dchr(norm(ei))*kron(ei, ei')/norm(ei)+f.chr(norm(ei))*I(3)
    fn = f.chr(norm(ei))*ei
    return inv(dsigma)*fn
end

function (f::ConductivityTD_b)(cell, cqdpt, mp)
    ei = f.efield[cell, cqdpt]
    if norm(ei)==0
        ei = 1e-9.+ei
    end
    #dsigma = f.dchr(norm(ei))*kron(ei, ei')/norm(ei)+f.chr(norm(ei))*I(3)
    #= if cell ==1 && cqdpt==1
        println(ei)
    end =#
    λ = 1.0
    margin = 0.2
    marginh = 0.5-margin
    xmp = mp.cart[1]-0.5
    ymp = mp.cart[2]-0.5
    #nearedge = false
    r = sqrt(xmp^2+ymp^2)
    λ = 1.0
    if r > marginh
        if r > 0.5
            λ = 0.0
        #= else
            λ = (r-marginh)/(margin) =#
        end
    end
    fn = (λ*f.chr(norm(ei))+(1-λ)*10.0)*ei
    return fn
end

BEAST.integrand(::NonLinearBEM.ConductivityTD_model, gx, ϕx) = gx[1] ⋅ ϕx
BEAST.defaultquadstrat(::NonLinearBEM.ConductivityTD_model, ::BEAST.LinearRefSpaceTriangle, ::BEAST.LinearRefSpaceTriangle) = BEAST.SingleNumQStrat(8)

BEAST.quaddata(exc::NonLinearBEM.ConductivityTD_model, g::BEAST.LinearRefSpaceTriangle, f::BEAST.LinearRefSpaceTriangle, tels, bels,
qs::BEAST.SingleNumQStrat) = 
NonLinearBEM.quaddata(exc::NonLinearBEM.ConductivityTD_model, g::BEAST.LinearRefSpaceTriangle, f::BEAST.LinearRefSpaceTriangle, tels, bels,
qs::BEAST.SingleNumQStrat)

function quaddata(exc::NonLinearBEM.ConductivityTD_model, g::BEAST.LinearRefSpaceTriangle, f::BEAST.LinearRefSpaceTriangle, tels, bels,
    qs::BEAST.SingleNumQStrat)

    u, w = BEAST.trgauss(qs.quad_rule)
    qd = [(w[i],SVector(u[1,i],u[2,i])) for i in 1:length(w)]
    A = BEAST._alloc_workspace(qd, g, f, tels, bels)

    return qd, A
end

BEAST.quadrule(exc::NonLinearBEM.ConductivityTD_model, ψ::BEAST.RefSpace, ϕ::BEAST.RefSpace, τ, (qd,A), qs::BEAST.SingleNumQStrat) = 
NonLinearBEM.quadrule(exc::NonLinearBEM.ConductivityTD_model, ψ::BEAST.RefSpace, ϕ::BEAST.RefSpace, τ, (qd,A), qs::BEAST.SingleNumQStrat)

function quadrule(exc::NonLinearBEM.ConductivityTD_model, ψ::BEAST.RefSpace, ϕ::BEAST.RefSpace, τ, (qd,A), qs::BEAST.SingleNumQStrat)
    for i in eachindex(qd)
        q = qd[i]
        w, p = q[1], BEAST.neighborhood(τ,q[2])
        A[i] = (w, p, ψ(p), ϕ(p))
    end
    return A
end

BEAST.kernelvals(f::NonLinearBEM.ConductivityTDop, mp, cell, cqdpt) = 
NonLinearBEM.kernelvals(f::ConductivityTDop, mp, cell, cqdpt)

function kernelvals(f::ConductivityTDop,mp, cell, cqdpt)
    ei = f.op.efield[cell, cqdpt]
    if norm(ei)==0
        ei = 1e-9.+ei
    end
    dsigma = f.op.dchr(norm(ei))*kron(ei, ei')/norm(ei)+f.op.chr(norm(ei))*I(3)
    if norm(inv(dsigma))>100
        println(norm(inv(dsigma)))
    end
    return inv(dsigma)
end

BEAST.kernelvals(f::NonLinearBEM.ConductivityTDop_b, mp, cell, cqdpt) = 
NonLinearBEM.kernelvals(f::ConductivityTDop_b, mp, cell, cqdpt)

function kernelvals(f::ConductivityTDop_b, mp, cell, cqdpt)
    ei = f.op.efield[cell, cqdpt]
    margin = 0.2
    marginh = 0.5-margin
    xmp = mp.cart[1]-0.5
    ymp = mp.cart[2]-0.5
    #nearedge = false
    r = sqrt(xmp^2+ymp^2)
    λ = 1.0
    if norm(ei)==0
        ei = 1e-9.+ei
    end
    if r > marginh
        if r > 0.5
            λ = 0.0
        #= else
            λ = (r-marginh)/(margin) =#
        end
    end
    #= println("not near edge")
    println(mp.cart) =#
    dsigma = λ*f.op.dchr(norm(ei))*kron(ei, ei')/norm(ei)+(λ*f.op.chr(norm(ei))+10.0*(1-λ))*I(3)
    return dsigma
end

BEAST.integrand(op::NonLinearBEM.ConductivityTDop_model, kernel, x, g, f) = dot(g[1],kernel*f[1])

BEAST.assemble_local_matched!(biop::NonLinearBEM.ConductivityTDop_b, tfs::BEAST.Space, bfs::BEAST.Space, store;
quadstrat=BEAST.defaultquadstrat(biop, tfs, bfs))=
NonLinearBEM.assemble_local_matched!(biop::NonLinearBEM.ConductivityTDop_b, tfs::BEAST.Space, bfs::BEAST.Space, store;
quadstrat=BEAST.defaultquadstrat(biop, tfs, bfs))

function assemble_local_matched!(biop::NonLinearBEM.ConductivityTDop_b, tfs::BEAST.Space, bfs::BEAST.Space, store;
    quadstrat=BEAST.defaultquadstrat(biop, tfs, bfs))

    tels, tad, ta2g = BEAST.assemblydata(tfs)
    bels, bad, ba2g = BEAST.assemblydata(bfs)

    bg2a = zeros(Int, length(BEAST.geometry(bfs)))
    for (i,j) in enumerate(ba2g) bg2a[j] = i end

    trefs = BEAST.refspace(tfs)
    brefs = BEAST.refspace(bfs)

    qd = BEAST.quaddata(biop, trefs, brefs, tels, bels, quadstrat)

    verbose = length(tels) > 10_000
    verbose && print("dots out of 20: ")
    todo, done, pctg = length(tels), 0, 0
    locmat = zeros(BEAST.scalartype(biop, trefs, brefs), BEAST.numfunctions(trefs), numfunctions(brefs))
    for (p,cell) in enumerate(tels)
        P = ta2g[p]
        q = bg2a[P]
        q == 0 && continue

        qr = BEAST.quadrule(biop, trefs, brefs, cell, qd, quadstrat)
        fill!(locmat, 0)
        BEAST.cellinteractions_matched!(locmat, biop, trefs, brefs, cell, qr,p)

        for i in 1 : size(locmat, 1), j in 1 : size(locmat, 2)
            for (m,a) in tad[p,i], (n,b) in bad[q,j]
                store(a * locmat[i,j] * b, m, n)
        
        end end

        new_pctg = round(Int, (done += 1) / todo * 100)
        verbose && new_pctg > pctg + 4 && (print("."); pctg = new_pctg)
    end
end

BEAST.cellinteractions_matched!(zlocal, biop::NonLinearBEM.ConductivityTDop_b, trefs, brefs, cell, qr, p) = 
NonLinearBEM.cellinteractions_matched!(zlocal, biop::NonLinearBEM.ConductivityTDop_b, trefs, brefs, cell, qr, p)

function cellinteractions(biop::ConductivityTDop_model, trefs::U, brefs::V, cell, qr, p) where {U<:RefSpace{T},V<:RefSpace{T}} where {T}

    num_tshs = length(qr[1][3])
    num_bshs = length(qr[1][4])

    zlocal = zeros(T, num_tshs, num_bshs)
    for (i,q) in enumerate(qr)

        w, mp, tvals, bvals = q[1], q[2], q[3], q[4]
        j = w * BEAST.jacobian(mp)
        kernel = BEAST.kernelvals(biop, mp, p, i)

        for m in 1 : num_tshs
            tval = tvals[m]

            for n in 1 : num_bshs
                bval = bvals[n]

                igd = BEAST.integrand(biop, kernel, mp, tval, bval)
                zlocal[m,n] += j * igd

            end
        end
    end

    return zlocal
end

BEAST.cellinteractions_matched!(zlocal, biop::ConductivityTDop_b, trefs, brefs, cell, qr, p) = 
NonLinearBEM.cellinteractions_matched!(zlocal, biop::ConductivityTDop_b, trefs, brefs, cell, qr, p)


function cellinteractions_matched!(zlocal, biop::ConductivityTDop_b, trefs, brefs, cell, qr, p)

    num_tshs = length(qr[1][3])
    num_bshs = length(qr[1][4])

    # zlocal = zeros(Float64, num_tshs, num_bshs)
    for (i,q) in enumerate(qr)

        w, mp, tvals, bvals = q[1], q[2], q[3], q[4]
        j = w * BEAST.jacobian(mp)
        kernel = BEAST.kernelvals(biop, mp, p,i)
        
        for n in 1 : num_bshs
            bval = bvals[n]
            for m in 1 : num_tshs
                tval = tvals[m]

                igd = BEAST.integrand(biop, kernel, mp, tval, bval)
                zlocal[m,n] += j * igd
            end
        end
    end

    return zlocal
end

BEAST.assemble!(field::NonLinearBEM.ConductivityTD_b, tfs::BEAST.Space, store;
quadstrat=BEAST.defaultquadstrat(field, tfs)) =
NonLinearBEM.assemble!(field::NonLinearBEM.ConductivityTD_b, tfs::BEAST.Space, store;
quadstrat=BEAST.defaultquadstrat(field, tfs))

function assemble!(field::ConductivityTD_b, tfs::BEAST.Space, store;
    quadstrat=BEAST.defaultquadstrat(field, tfs))

    tels, tad = BEAST.assemblydata(tfs)

    trefs = BEAST.refspace(tfs)
    qd = BEAST.quaddata(field, trefs, tels, quadstrat)

    for (t, tcell) in enumerate(tels)

        # compute the testing with the reference elements
        qr = BEAST.quadrule(field, trefs, t, tcell, qd, quadstrat)
        blocal = BEAST.celltestvalues(trefs, t, tcell, field, qr)

        for i in 1 : BEAST.numfunctions(trefs)
            for (m,a) in tad[t,i]
                store(a*blocal[i], m)
            end
        end

    end

end

BEAST.celltestvalues(tshs::BEAST.RefSpace{T, NF}, t, tcell, field::NonLinearBEM.ConductivityTD_b, qr) where {T, NF} = 
NonLinearBEM.celltestvalues(tshs::BEAST.RefSpace{T, NF}, t, tcell, field::NonLinearBEM.ConductivityTD_b, qr)

function celltestvalues(tshs::BEAST.RefSpace{T, NF}, t, tcell, field::NonLinearBEM.ConductivityTD_b, qr) where {T, NF}

    num_tshs = numfunctions(tshs)
    interactions = zeros(Complex{T}, num_tshs)

    num_oqp = length(qr)

    for p in 1 : num_oqp
        mp = qr[p].point

        dx =qr[p].weight

        fval = field(t,p,mp)
        tvals = qr[p].value

        for m in 1 : num_tshs
            tval = tvals[m]

            igd = BEAST.integrand(field, tval, fval)
            interactions[m] += igd * dx
        end
    end

    return interactions
end
#= function momintegrals!(z, exc::ConductivityTD, testrefs, timerefs, ncell, τ, ntcell, ρ, qr)

    for (ns,p) in enumerate(qr.quad_points)
        x = p.point
        w = p.weight
        f = p.value
        dx = w

        # try
        #     @assert ρ.vertices[1][1] <= cartesian(x)[1] <= ρ.vertices[2][1]
        # catch
        #     @show ρ.vertices[1][1]
        #     @show cartesian(x)[1]
        #     @show ρ.vertices[2][1]
        #     error("")
        # end

        tqr = timequadrule(qr,p)
        timeintegrals!(z, exc, testrefs, timerefs, x, ncell, ns, ntcell, ρ, dx, tqr, f)

    end

end

function timeintegrals!(z, exc::ConductivityTD, testrefs, timerefs, testpoint, ncell, ns, ntcell, timeelement, dx, qr, f)

    for (nt,p) in enumerate(qr.quad_points)
        t = p.point
        w = p.weight
        U = p.value
        dt = w #* jacobian(t) # * volume(timeelement)

        for i in 1 : numfunctions(testrefs)
            for k in 1 : numfunctions(timerefs)
                z[i,k] += dot(f[i][1]*U[k], exc(τ,ns, timeelement, nt)) * dt * dx
            end
        end
    end
end

function timeintegrals!(z, exc::ConductivityTD,
    spacerefs, timerefs::DiracBoundary,
    testpoint, ncell, ns, ntcell, timeelement,
    dx, qr, testvals)

    # since timeelement uses barycentric coordinates,
    # the first/left vertex has coords u = 1.0!
    testtime = neighborhood(timeelement, point(0.0))
    @assert cartesian(testtime)[1] ≈ timeelement.vertices[2][1]

    for i in 1 : numfunctions(spacerefs)
        z[i,1] += dot(testvals[i][1], exc(ncell, ns, ntcell)) * dx
    end
end =#
  
  #= planewave(;signature, polarization, direction, speedoflight) =
      PlaneWaveMWTD(direction, polarization, speedoflight, signature)
   =#

  #= 
  *(a, pw::PlaneWaveMWTD) = PlaneWaveMWTD(
      pw.direction,
      a * pw.polarisation,
      pw.speedoflight,
      pw.amplitude
  ) =#
  
  #= cross(k, pw::PlaneWaveMWTD) = PlaneWaveMWTD(
      pw.direction,
      k × pw.polarisation,
      pw.speedoflight,
      pw.amplitude
  )
   =#

  #= function integrate(f::BEAST.PlaneWaveMWTD)
      planewave(
          signature = integrate(f.amplitude),
          direction = f.direction,
          polarization = f.polarisation,
          speedoflight = f.speedoflight)
  end
  
  function differentiate(f::BEAST.PlaneWaveMWTD)
      planewave(
          signature = derive(f.amplitude),
          direction = f.direction,
          polarization = f.polarisation,
          speedoflight = f.speedoflight)
  end =#
  