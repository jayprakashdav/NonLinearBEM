#abstract type TDLocalOperator <: LocalOperator end
#= mutable struct σ2DMaterial{A,B,C,D,E,F,G,H,I,J} <: LocalOperator
	chr::BEAST.Polynomial
	numdiffs::Int
    timeindex::Int
	coeffs::Array{Float64, 2}
    numrefs::Int
    tnumrefs::Int
    field::SVector{3, Float64}
    ptlocs::Vector{Int64}
    lastloc::Int
    charts::Vector{A}
    tcells::Vector{B}
    pts::Vector{C}
    space_basis::D
    time_basis::E
    refs::F
    trefs::G
    ad::H
    tad::I
    chart_tree::J
    #= space_basis::BEAST.AbstractSpace{A}
    time_basis::BEAST.AbstractTimeBasisFunction{B}
    refs::BEAST.RefSpace{Float64}
    trefs::BEAST.RefSpace{Float64}
    ad::BEAST.AssemblyData{Float64}
    tad::BEAST.AssemblyData{Float64}
    chart_tree::CollisionDetection.Octree{F} =#
end =#
mutable struct σ2DMaterial <: BEAST.LocalOperator
    chr::BEAST.Polynomial
    numdiffs::Int
    field::Array{Float64,2}
end

function NLConductivity(cond;numdiffs=0)
    σ2DMaterial(cond,numdiffs,zeros(1,1))
end

#BEAST.kernelvals(op::σ2DMaterial,element,qdpt) = op.chr(norm(op.field[element,qdpt]))
BEAST.kernelvals(op::σ2DMaterial, mp, cell, qdpt) = op.chr(op.field[cell,qdpt])
#BEAST.kernelvals(op::σ2DMaterial, mp) = 50.0
BEAST.integrand(op::σ2DMaterial, kernel, x, g, f) = kernel*dot(g[1],f[1])
BEAST.scalartype(op::σ2DMaterial) = Float64

function update!(op::σ2DMaterial,coeffs,k,eq)
    #= testST = eq.test_space_dict[1]
    trialST = eq.trial_space_dict[1]
	strat = defaultquadstrat(op, testST.space, trialST.space) =#
	quadpoint_refines!(op,coeffs,k, eq.test_space_dict[1], eq.trial_space_dict[1])
end

#= function Base.:*(a::Number, op::σ2DMaterial)
	@info "scalar product a * op (σ2DMaterial)"
	σ2DMaterial(
		a*op.chr,
		op.numdiffs,
		a.*op.field)
end =#

function mpeval(point, op::σ2DMaterial; type=nothing)
    tcell = op.tcells[op.timeindex]
    k = op.timeindex
    #i = CompScienceMeshes.findchart(op.charts, op.chart_tree, point)
    i=nothing
    if k==1
        i = CompScienceMeshes.findchart(op.charts, op.chart_tree, point)
        push!(op.ptlocs, i)
        push!(op.pts, cartesian(point))
    else
        op.lastloc +=1
        if isapprox(op.pts[op.lastloc],cartesian(point), atol=1e-12)
            i=op.ptlocs[op.lastloc]
        else
            println("incorrect point related from last iteration")
            i = CompScienceMeshes.findchart(op.charts, op.chart_tree, point)
        end
    end
    if i !== nothing
        # @show i
        chart = op.charts[i]
        u = carttobary(chart, point)
        @assert isapprox(cartesian(point), barytocart(chart, u), atol=1e-12)
        vals = op.refs(neighborhood(chart,u))
        tvals = op.trefs(neighborhood(tcell,1))
        op.field[1] = [0.0,0.0,0.0]
        for si in 1:op.numrefs
            for ti in 1:op.tnumrefs
                for (m,a) in op.ad[i,si]
                    for (n,b) in op.tad[k,ti]
                        fx = vals[si].value
                        tfx = tvals[ti]
                        op.field[1] += (op.coeffs[m,n] * tfx * b) * a * fx
                    end
                end
            end
        end
    end
    return op.field[1]
end

function quadpoint_refines!(biop::LocalOperator, coeffs, k, tfs, bfs;
        quadstrat=BEAST.defaultquadstrat(biop, tfs.space, bfs.space))
    
        tgeo = BEAST.geometry(tfs.space)
        bgeo = BEAST.geometry(bfs.space)
        @assert CompScienceMeshes.refines(tgeo, bgeo)
    
        trefs = BEAST.refspace(tfs.space)
        brefs = BEAST.refspace(bfs.space)
        btrefs = BEAST.refspace(bfs.time)

        numbrefs = BEAST.numfunctions(brefs)
        numbtrefs = BEAST.numfunctions(btrefs)
    
        tels, tad, ta2g = BEAST.assemblydata(tfs.space)
        bels, bad, ba2g = BEAST.assemblydata(bfs.space)
        btels, btad = BEAST.assemblydata(bfs.time)

        bg2a = zeros(Int, length(BEAST.geometry(bfs.space)))
        for (i,j) in enumerate(ba2g) bg2a[j] = i end
    
        qd = BEAST.quaddata(biop, trefs, brefs, tels, bels, quadstrat)
        if k==0
            biop.field = zeros(length(tels),length(qd[1]))
            k=1
        end
        btcell = btels[k]
        todo, done, pctg = length(tels), 0, 0
        for (p,tcell) in enumerate(tels)
    
            P = ta2g[p]
            Q = CompScienceMeshes.parent(tgeo, P)
            q = bg2a[Q]
    
            bcell = bels[q]
            @assert overlap(tcell, bcell)
    
            isct = intersection(tcell, bcell)
            for cell in isct    
                qr = BEAST.quadrule(biop, trefs, brefs, cell, qd, quadstrat)
                for (qi,qdpt) in enumerate(qr)
                    mp = carttobary(bcell, qdpt[2])
                    vals = brefs(neighborhood(bcell, mp))
                    tvals = btrefs(neighborhood(btcell,1))
                    temp = [0.0,0.0,0.0]
                    for si in 1:numbrefs
                        for ti in 1:numbtrefs
                            for (m,a) in bad[q,si]
                                for (n,b) in btad[k,ti]
                                    fx = vals[si].value
                                    tfx = tvals[ti]
                                    temp += (coeffs[m,n] * tfx * b) * a * fx
                                end
                            end
                        end
                    end
                    biop.field[p,qi] = norm(temp)
                end
            end # next cell in intersection
        end # next cell in the test geometry
end

function tdsolve(eq1::BEAST.DiscreteEquation, eq2::BEAST.DiscreteEquation)
    Z = BEAST.td_assemble(eq1.equation.lhs, eq1.test_space_dict, eq1.trial_space_dict)
    b = BEAST.td_assemble(eq1.equation.rhs, eq1.test_space_dict)
    h = eq2.trial_space_dict[eq2.equation.lhs.terms[1].trial_id]
    h1 = h.space⊗BEAST.derive(h.time)
    f1 = eq1.test_space_dict[1]
    f2 = eq2.test_space_dict[1]
    g = eq1.trial_space_dict[1]
    idST = BEAST.Identity()⊗BEAST.Identity()
    Ġ = BEAST.assemble(idST, f1, h1)
    G_j = BEAST.assemble(idST, f2, g)
    Nt = BEAST.numfunctions(g.time)
    if typeof(eq2.equation.rhs.terms[1].functional)==ConductivityTD
        G_nl = BEAST.assemble(idST, f2, h)
        return marchonintimenl(eq1, eq2, Z, b, Ġ, G_j, G_nl, Nt)
    elseif typeof(eq2.equation.rhs.terms[1].functional)==ConductivityTD_b
        G_nl = BEAST.assemble(idST, f2, h)
        if eq2.equation.rhs.terms[1].functional.initialcondition
            return marchonintimenl5(eq1, eq2, Z, b, Ġ, G_j, G_nl, Nt)
        else
            return marchonintimenl4(eq1, eq2, Z, b, Ġ, G_j, G_nl, Nt)
        end
    end
    marchonintime(eq1, eq2, Z, b, Ġ, G_j, Nt)
end

function marchonintime(eq1, eq2,  Z, inc, Ġ, G_j, Nt)
    T = eltype(Z)
    Z0 = zeros(T, size(Z)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Z0,Z,1)
    Ġ0 = Ġ.data[1,:,:]
    G_j0 = G_j.data[1,:,:]
    T = eltype(Z0)
    M,N = size(Z0)
    Me,Ne = size(Ġ0)
    @assert M == size(inc,1)
    xj = zeros(T,N,Nt)
	xe = zeros(T,Ne,Nt)
    yj = zeros(T,N)
	ye = zeros(T,N)
    jk_start = 2
    ek_start = 2
    jk_stop = Nt
    ek_stop = BEAST.numfunctions(eq1.trial_space_dict[1].time)+1
    csxj = zeros(T,N,Nt)
    csxe = zeros(T,Ne,Nt)
    σ = eq2.equation.lhs.terms[1].kernel
    Q = zeros(eltype(Z0), size(Z0))
    for i in 1:Nt
        R = inc[:,i]
		fill!(yj,0)
		fill!(ye,0)
        update!(σ.spatial_factor,xe,i-1,eq2)
        Q = BEAST.assemble(σ.spatial_factor, eq2.test_space_dict[1].space, eq2.trial_space_dict[1].space)
		BEAST.ConvolutionOperators.convolve!(yj,Z,xj,csxj,i,jk_start,jk_stop)
        BEAST.ConvolutionOperators.convolve!(ye,Ġ,xe,csxe,i,ek_start,ek_stop)
        b = R - yj - ye         #solve
        P = inv(Matrix(Q))      #|Z   -Ġ|J =|b|
        V0 = inv(Z0+Ġ0*P*G_j0)  #|G_j  Q|E =|0|
        xj[:,i] = V0*b
        xe[:,i] = P*G_j0*xj[:,i]
        #= for iter in 1:3
            update!(σ.spatial_factor,xe,i,eq2)
            Q = assemble(σ.spatial_factor, eq2.test_space_dict[1].space, eq2.trial_space_dict[1].space)
            b = R - yj + ye         #solve
            P = inv(Matrix(Q))      #|Z   -Ġ|J =|b|
            V0 = inv(Z0-Ġ0*P*G_j0)  #|G_j Q|E =|0|
            xj[:,i] = V0*b
            xe[:,i] = P*G_j0*xj[:,i]
        end =#
        if i > 1
            csxj[:,i] .= csxj[:,i-1] .+ xj[:,i]
            csxe[:,i] .= csxe[:,i-1] .+ xe[:,i]
        else
            csxj[:,i] .= xj[:,i]
            csxe[:,i] .= xe[:,i]
        end
        (i % 10 == 0) && print(i, "[", Nt, "] - ")
    end
    return xj, xe
end
## Modified facecurrents to make it able to compute face charges##
function facecharges(coeffs, basis::SpaceTimeBasis)

    space_basis = basis.space
    time_basis = basis.time

    Nt = BEAST.numfunctions(time_basis)
    Δt = BEAST.timestep(time_basis)

    refs = BEAST.refspace(space_basis)
    trefs = BEAST.refspace(time_basis)
    numrefs = BEAST.numfunctions(refs)
    tnumrefs = BEAST.numfunctions(trefs)

    cells, ad = BEAST.assemblydata(space_basis)
    tcells, tad = BEAST.assemblydata(time_basis)

    mesh = BEAST.geometry(space_basis)
    T = eltype(coeffs)
    D = dimension(mesh)
    U = 1

    # TODO: express relative to input types
    PT = SVector{U, T}
    fcr = zeros(T, numcells(mesh), Nt)

    for (k,tcell) in enumerate(tcells)
        tmps = neighborhood(tcell,1)
        tvals = trefs(tmps)
        for (p,cell) in enumerate(cells)
            mps = center(cell)
            vals = refs(mps)

            # assemble
            for i in 1:numrefs
                fx = vals[i][1]
                for j in 1:tnumrefs
                    tfx = tvals[j]
                    for (m,a) in ad[p,i]
                        for (n,b) in tad[k,j]
                            fcr[p,k] += (coeffs[m,n] * tfx * b) * a * fx
                        end
                    end
                end
            end
        end
    end

    return fcr, mesh
end

function cellinteractions(biop::σ2DMaterial, trefs::U, brefs::V, cell, qr, p) where {U<:RefSpace{T},V<:RefSpace{T}} where {T}

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
#= #have to duplicate the function definition intended for assembly of Identity operator also for IdentityTime operator.
#The scalartype of Identiy operator is not defined probably
#because it is meant to be used with other kind of operator whose scalartype will be inherited when needed
function assemble(op::IdentityTime,
    testnfs::AbstractTimeBasisFunction,
    trialfns::AbstractTimeBasisFunction)

    tbf = convolve(testnfs, trialfns)
    has_zero_tail = all(tbf.polys[end].data .== 0)
    # @show has_zero_tail

    T = scalartype(tbf)
    if has_zero_tail
        z = zeros(T, numintervals(tbf)-1)
    else
        z = zeros(T, numfunctions(tbf))
    end

    Δt = timestep(tbf)
    #for i in eachindex(z)
    for i in 1:numintervals(tbf)-1
        p = tbf.polys[i]
        t = (i-1)*Δt
        z[i] = evaluate(p,t)
    end

    for i in numintervals(tbf):length(z)
        p = tbf.polys[end]
        t = (i-1)*Δt
        z[i] = evaluate(p,t)
    end

    return z
end =#