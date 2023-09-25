function update!(op::ConductivityTD_model, jcoeffs, ecoeffs, k, eq1, eq2)
    tgeo = BEAST.geometry(eq2.test_space_dict[1].space)
    bgeo = BEAST.geometry(eq2.trial_space_dict[1].space)
    if CompScienceMeshes.refines(tgeo, bgeo)
        op.efield = quadpoint_field_refines(op, ecoeffs, k, eq2.test_space_dict[1], eq2.trial_space_dict[1])
        op.jflux = quadpoint_field_refines(op, jcoeffs, k, eq2.test_space_dict[1], eq1.trial_space_dict[1])
    else
        op.efield = quadpoint_field(op, ecoeffs, k, eq2.test_space_dict[1], eq2.trial_space_dict[1])
        op.jflux = quadpoint_field(op, jcoeffs, k, eq2.test_space_dict[1], eq1.trial_space_dict[1])
    end	
end

function mpeval(point, op::ConductivityTD_model; type=nothing)
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

function quadpoint_field_refines(biop::ConductivityTD_model, coeffs, k, tfs, bfs;
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
        T = eltype(coeffs)
        D = dimension(BEAST.geometry(bfs.space))
        U = D+1
        PT = SVector{U, T}
        field = zeros(PT, (length(tels),length(qd[1])))
        btcell = btels[k]
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
                    for si in 1:numbrefs
                        fx = vals[si].value
                        for ti in 1:numbtrefs
                            tfx = tvals[ti]
                            for (m,a) in bad[q,si]
                                for (n,b) in btad[k,ti]
                                    field[p,qi] += (coeffs[m,n] * tfx * b) * a * fx
                                end
                            end
                        end
                    end
                end
            end # next cell in intersection
        end # next cell in the parent geometry
        return field
end

function quadpoint_field(biop::ConductivityTD_model, coeffs, k, tfs, bfs;
    quadstrat=BEAST.defaultquadstrat(biop, tfs.space, bfs.space))

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
    T = eltype(coeffs)
    D = dimension(BEAST.geometry(bfs.space))
    U = D+1
    PT = SVector{U, T}
    field = zeros(PT, (length(tels),length(qd[1])))
    btcell = btels[k]
    for (p,bcell) in enumerate(bels)   
        qr = BEAST.quadrule(biop, trefs, brefs, bcell, qd, quadstrat)
        #println(length(qr))
        for (qi,qdpt) in enumerate(qr)
            mp = carttobary(bcell, qdpt[2])
            vals = brefs(neighborhood(bcell, mp))
            tvals = btrefs(neighborhood(btcell,1))
            for si in 1:numbrefs
                fx = vals[si].value
                for ti in 1:numbtrefs
                    tfx = tvals[ti]
                    for (m,a) in bad[p,si]
                        for (n,b) in btad[k,ti]
                            field[p,qi] += (coeffs[m,n] * tfx * b) * a * fx
                        end
                    end
                end
            end
        end
    end
    return field
end

function marchonintimenl(eq1, eq2,  Z, inc, Ġ, G_j, G_nl, Nt)
    T = eltype(Z)
    Z0 = zeros(T, size(Z)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Z0,Z,1)
    Ġ0 = Ġ.data[1,:,:]
    G_j0 = G_j.data[1,:,:]
    G_nl0 = G_nl.data[1,:,:]
    T = eltype(Z0)
    M,N = size(Z0)
    Me,Ne = size(Ġ0)
    @assert M == size(inc,1)
    xj = zeros(T,N,Nt)
	xe = zeros(T,Ne,Nt)
	xeprev = zeros(T,N)
    yj = zeros(T,N)
	ye = zeros(T,N)
    jk_start = 2
    ek_start = 2
    jk_stop = Nt
    ek_stop = BEAST.numfunctions(eq1.trial_space_dict[1].time)+1
    csxj = zeros(T,N,Nt)
    csxe = zeros(T,Ne,Nt)
    σ = eq2.equation.rhs.terms[1].functional
    σop = BEAST.ConductivityTDop(σ)
    bσ = zeros(T, N)
    G_nl_inv = inv(Matrix(G_nl0))
    for i in 1:Nt
        R = inc[:,i]
		fill!(yj,0)
		fill!(ye,0)
		BEAST.ConvolutionOperators.convolve!(yj,Z,xj,csxj,i,jk_start,jk_stop)
        BEAST.ConvolutionOperators.convolve!(ye,Ġ,xe,csxe,i,ek_start,ek_stop)
        iter_max = 10
        for l in 1:iter_max
            if l==1
                if i==1
                    i=2
                end
                xeprev = xe[:,i-1]
                update!(σ, xj, xe, i-1, eq1, eq2)
            else
                xeprev = xe[:,i]
                update!(σ, xj, xe, i, eq1, eq2)
            end
            bσ = BEAST.assemble(σ, eq2.test_space_dict[1].space)
            Q = BEAST.assemble(σop, eq2.test_space_dict[1].space, eq1.trial_space_dict[1].space)
            #= if i>40
                println(i)
                println(norm(Q))
                #println(norm(bσ))
            end =#
            V0 = inv(Z0 - Ġ0*G_nl_inv*Q)
            b = R - yj + ye + Ġ0*xeprev - Ġ0*G_nl_inv*bσ         #solve |Z   -Ġ|J =|b|                     
            xj[:,i] = V0*b                            #|G_j  Q|E =|0|
            #xe[:,i] = xeprev - G_e0_inv*bσ
            xe[:,i] = xeprev - G_nl_inv*bσ + G_nl_inv*Q*xj[:,i]
            if norm(xe[:,i]-xeprev) < 1e-6
                break
            end
        end
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

function marchonintimenl2(eq1, eq2,  Z, inc, Ġ, G_j, G_nl, Nt)
    T = eltype(Z)
    Z0 = zeros(T, size(Z)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Z0,Z,1)
    Ġ0 = Ġ.data[1,:,:]
    G_j0 = G_j.data[1,:,:]
    G_nl0 = G_nl.data[1,:,:]
    T = eltype(Z0)
    M,N = size(Z0)
    Me,Ne = size(Ġ0)
    @assert M == size(inc,1)
    xj = zeros(T,N,Nt)
	xe = zeros(T,Ne,Nt)
	xeprev = zeros(T,N)
    yj = zeros(T,N)
	ye = zeros(T,N)
    jk_start = 2
    ek_start = 2
    jk_stop = Nt
    ek_stop = BEAST.numfunctions(eq1.trial_space_dict[1].time)+1
    csxj = zeros(T,N,Nt)
    csxe = zeros(T,Ne,Nt)
    σ = eq2.equation.rhs.terms[1].functional
    σop = BEAST.ConductivityTDop_b(σ)
    bσ = zeros(T, N)
    G_nl_inv = inv(Matrix(G_nl0))
    for i in 1:Nt
        println(i)
        R = inc[:,i]
		fill!(yj,0)
		fill!(ye,0)
		BEAST.ConvolutionOperators.convolve!(yj,Z,xj,csxj,i,jk_start,jk_stop)
        BEAST.ConvolutionOperators.convolve!(ye,Ġ,xe,csxe,i,ek_start,ek_stop)
        iter_max = 2
        for l in 1:iter_max
            if l==1
                if i==1
                    i=2
                end
                xeprev = xe[:,i-1]
                update!(σ, xj, xe, i-1, eq1, eq2)
            else
                xeprev = xe[:,i]
                update!(σ, xj, xe, i, eq1, eq2)
            end
            bσ = BEAST.assemble(σ, eq2.test_space_dict[1].space)
            Q = BEAST.assemble(σop, eq2.test_space_dict[1].space, eq2.trial_space_dict[1].space)
            Q_inv = inv(Matrix(Q))
            #= println("normQ ", norm(Q))
            println("normbsigma ", norm(bσ)) =#
            V0 = inv(Z0 - Ġ0*Q_inv*G_j0)
            rhs1 = R - yj + ye
            rhs2 = -Ġ0*Q_inv*bσ + Ġ0*xeprev
            b = rhs1 + rhs2         #solve |Z   -Ġ|J =|b|                     
            xj[:,i] = V0*b                            #|G_j  Q|E =|0|
            xe[:,i] = Q_inv*G_j0*xj[:,i]+ xeprev - Q_inv*bσ
            println("norm xe ", norm(xe[:,i]-xeprev))
            if norm(xe[:,i]-xeprev) < 1e-6
                break
            end
        end
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

function marchonintimenl3(eq1, eq2,  Z, inc, Ġ, G_j, G_nl, Nt)
    T = eltype(Z)
    Z0 = zeros(T, size(Z)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Z0,Z,1)
    Ġ0 = Ġ.data[1,:,:]
    G_j0 = G_j.data[1,:,:]
    G_nl0 = G_nl.data[1,:,:]
    T = eltype(Z0)
    M,N = size(Z0)
    Me,Ne = size(Ġ0)
    @assert M == size(inc,1)
    xj = zeros(T,N,Nt)
	xe = zeros(T,Ne,Nt)
	xeprev = zeros(T,N)
    yj = zeros(T,N)
	ye = zeros(T,N)
    jk_start = 2
    ek_start = 2
    jk_stop = Nt
    ek_stop = BEAST.numfunctions(eq1.trial_space_dict[1].time)+1
    csxj = zeros(T,N,Nt)
    csxe = zeros(T,Ne,Nt)
    σ = eq2.equation.rhs.terms[1].functional
    σop = BEAST.ConductivityTDop_b(σ)
    bσ = zeros(T, N)
    G_nl_inv = inv(Matrix(G_nl0))
    G_j0_inv = inv(Matrix(G_j0))
    println(Nt)
    for i in 1:Nt
        println(i)
        R = inc[:,i]
		fill!(yj,0)
		fill!(ye,0)
		BEAST.ConvolutionOperators.convolve!(yj,Z,xj,csxj,i,jk_start,jk_stop)
        BEAST.ConvolutionOperators.convolve!(ye,Ġ,xe,csxe,i,ek_start,ek_stop)
        iter_max = 100
        for l in 1:iter_max
            if l==1
                if i==1
                    i=2
                end
                xeprev = xe[:,i-1]
                update!(σ, xj, xe, i-1, eq1, eq2)
            else
                xeprev = xe[:,i]
                update!(σ, xj, xe, i, eq1, eq2)
            end
            bσ = BEAST.assemble(σ, eq2.test_space_dict[1].space)
            Q = BEAST.assemble(σop, eq2.test_space_dict[1].space, eq2.trial_space_dict[1].space)
            #Q_inv = inv(Matrix(Q))
            #= println("normQ ", norm(Q))
            println("normbsigma ", norm(bσ)) =#
            V0 = inv(Z0*G_j0_inv*Q - Ġ0)
            rhs1 = R - yj + ye
            rhs2 = -Z0*G_j0_inv*bσ + Z0*G_j0_inv*Q*xeprev
            b = rhs1 + rhs2         #solve |Z   -Ġ|J =|b|                     
            xe[:,i] = V0*b                            #|G_j  Q|E =|0|
            xj[:,i] = G_j0_inv*Q*xe[:,i] - rhs2
            println("norm xe ", norm(xe[:,i]-xeprev))
            if norm(xe[:,i]-xeprev) < 1e-6
                break
            end
        end
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

function marchonintimenl4(eq1, eq2,  Z, inc, Ġ, G_j, G_nl, Nt)
    T = eltype(Z)
    Z0 = zeros(T, size(Z)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Z0,Z,1)
    Ġ0 = Ġ.data[1,:,:]
    G_j0 = G_j.data[1,:,:]
    G_nl0 = G_nl.data[1,:,:]
    T = eltype(Z0)
    M,N = size(Z0)
    Me,Ne = size(Ġ0)
    V0 = zeros(N+Ne,N+Ne)
    V0[1:N, 1:N] = Z0
    V0[1:N, N+1:N+Ne] = -Ġ0
    V0[N+1:N+Ne, 1:N] = G_j0
    sol = zeros(2*N)
    xj_all = zeros(N)
    xe_all = zeros(N)
    @assert M == size(inc,1)
    xj = zeros(T,N,Nt)
	xe = zeros(T,Ne,Nt)
	xeprev = zeros(T,N)
    yj = zeros(T,N)
	ye = zeros(T,N)
    jk_start = 2
    ek_start = 2
    jk_stop = Nt
    ek_stop = BEAST.numfunctions(eq1.trial_space_dict[1].time)+1
    csxj = zeros(T,N,Nt)
    csxe = zeros(T,Ne,Nt)
    σ = eq2.equation.rhs.terms[1].functional
    σop = NonLinearBEM.ConductivityTDop_b(σ)
    bσ = zeros(T, N)
    #G_nl_inv = inv(Matrix(G_nl0))
    #G_j0_inv = inv(Matrix(G_j0))
    println(Nt)
    for i in 1:Nt
        println(i)
        R = inc[:,i]
        ##do a clean implementation
        #= g = BEAST.creategaussian(2.9504, 2.2128, 0.8)
        R[end] = g(0.01844*i) =#
		fill!(yj,0)
		fill!(ye,0)
		BEAST.ConvolutionOperators.convolve!(yj,Z,xj,csxj,i,jk_start,jk_stop)
        BEAST.ConvolutionOperators.convolve!(ye,Ġ,xe,csxe,i,ek_start,ek_stop)
        iter_max = 200
        for l in 1:iter_max
            if l==1
                if i==1
                    i=2
                end
                xeprev = xe[:,i-1]
                update!(σ, xj, xe, i-1, eq1, eq2)
            else
                xeprev = xe[:,i]
                update!(σ, xj, xe, i, eq1, eq2)
            end
            bσ = BEAST.assemble(σ, eq2.test_space_dict[1].space)
            Q = BEAST.assemble(σop, eq2.test_space_dict[1].space, eq2.trial_space_dict[1].space)
            #Q_inv = inv(Matrix(Q))
            #= println("normQ ", norm(Q))
            println("normbsigma ", norm(bσ)) =#
            #V0 = inv(Z0*G_j0_inv*Q - Ġ0)
            V0[N+1:N+Ne,N+1:N+Ne] = -Q
            rhs1 = R - yj + ye
            rhs2 = bσ - Q*xeprev
            #b = rhs1 + rhs2         #solve |Z   -Ġ|J =|b|                     
            b = [rhs1; rhs2]
            sol = inv(Matrix(V0))*b
            #= invV0 = GMRESSolver(V0, maxiter=1000, restart=0, tol=1e-6)
            mul!(sol, invV0, b) =#
            xj[:,i] = sol[1:N]                            #|G_j  Q|E =|0|
            xe[:,i] = sol[N+1:end]
            #xe_all = hcat(xe_all, xe[:,i])
            #xj_all = hcat(xj_all, xj[:,i])
            println("norm xe ", norm(xe[:,i]-xeprev))
            if norm(xe[:,i]-xeprev) < 1e-4
                break
            end
        end
        if i > 1
            csxj[:,i] .= csxj[:,i-1] .+ xj[:,i]
            csxe[:,i] .= csxe[:,i-1] .+ xe[:,i]
        else
            csxj[:,i] .= xj[:,i]
            csxe[:,i] .= xe[:,i]
        end
        (i % 10 == 0) && print(i, "[", Nt, "] - ")
    end
    return xj, xe, xj_all, xe_all
end

#initial condition implementation
function marchonintimenl5(eq1, eq2,  Z, inc, Ġ, G_j, G_nl, Nt)
    T = eltype(Z)
    Z0 = zeros(T, size(Z)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Z0,Z,1)
    Ġ0 = Ġ.data[1,:,:]
    G_j0 = G_j.data[1,:,:]
    G_nl0 = G_nl.data[1,:,:]
    T = eltype(Z0)
    M,N = size(Z0)
    Me,Ne = size(Ġ0)
    @assert M == size(inc,1)
    xj = zeros(T,N,Nt)
	xe = zeros(T,Ne,Nt)
	xeprev = zeros(T,N)
    yj = zeros(T,N)
	ye = zeros(T,N)
    jk_start = 2
    ek_start = 2
    jk_stop = Nt
    ek_stop = BEAST.numfunctions(eq1.trial_space_dict[1].time)+1
    csxj = zeros(T,N,Nt)
    csxe = zeros(T,Ne,Nt)
    σ = eq2.equation.rhs.terms[1].functional
    σop = BEAST.ConductivityTDop_b(σ)
    bσ = zeros(T, N)
    init_steps = size(σ.init_e,2)
    @assert N = size(σ.init_e,1)
    @assert N = size(σ.init_j,1)
    println(Nt)
    xj[1:init_steps] = init_j
    xe[1:init_steps] = init_e
    for i in init_steps+1:Nt
        println(i)
        R = inc[:,i]
		fill!(yj,0)
		fill!(ye,0)
		BEAST.ConvolutionOperators.convolve!(yj,Z,xj,csxj,i,jk_start,jk_stop)
        BEAST.ConvolutionOperators.convolve!(ye,Ġ,xe,csxe,i,ek_start,ek_stop)
        iter_max = 100
        for l in 1:iter_max
            if l==1
                if i==1
                    i=2
                end
                xeprev = xe[:,i-1]
                update!(σ, xj, xe, i-1, eq1, eq2)
            else
                xeprev = xe[:,i]
                update!(σ, xj, xe, i, eq1, eq2)
            end
            bσ = BEAST.assemble(σ, eq2.test_space_dict[1].space)
            Q = BEAST.assemble(σop, eq2.test_space_dict[1].space, eq2.trial_space_dict[1].space)
            #Q_inv = inv(Matrix(Q))
            #= println("normQ ", norm(Q))
            println("normbsigma ", norm(bσ)) =#
            #V0 = inv(Z0*G_j0_inv*Q - Ġ0)
            V0 = [Z0 -Ġ0; G_j0 -Q]
            rhs1 = R - yj + ye
            rhs2 = bσ - Q*xeprev
            #b = rhs1 + rhs2         #solve |Z   -Ġ|J =|b|                     
            b = [rhs1; rhs2]
            sol = inv(Matrix(V0))*b
            xj[:,i] = sol[1:N]                            #|G_j  Q|E =|0|
            xe[:,i] = sol[N+1:end]
            println("norm xe ", norm(xe[:,i]-xeprev))
            if norm(xe[:,i]-xeprev) < 1e-5
                break
            end
        end
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

#= function marchonintimenl6(eq1, eq2,  Z, inc, Ġ, G_j, G_nl, Nt)
    T = eltype(Z)
    Z0 = zeros(T, size(Z)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Z0,Z,1)
    Ġ0 = Ġ.data[1,:,:]
    G_j0 = G_j.data[1,:,:]
    G_nl0 = G_nl.data[1,:,:]
    T = eltype(Z0)
    M,N = size(Z0)
    Me,Ne = size(Ġ0)
    @assert M == size(inc,1)
    xj = zeros(T,N,Nt)
	xe = zeros(T,Ne,Nt)
	xeprev = zeros(T,N)
    yj = zeros(T,N)
	ye = zeros(T,N)
    jk_start = 2
    ek_start = 2
    jk_stop = Nt
    ek_stop = BEAST.numfunctions(eq1.trial_space_dict[1].time)+1
    csxj = zeros(T,N,Nt)
    csxe = zeros(T,Ne,Nt)
    σ = eq2.equation.rhs.terms[1].functional
    σop = BEAST.ConductivityTDop_b(σ)
    bσ = zeros(T, N)
    padrow1 = zeros(T, N)
    padrow2 = zeros(T, N)
    padrow2[end] = -1
    padmat = [1.0 0.0; -10.0 -1.0]
    Z01 = [Z0 padrow1 padrow2; [[padrow2'; padrow1'] padmat]]
    Ġ01 = [Ġ0; padrow1'; padrow1']
    G_j01 = [G_j0 padrow1 padrow1]
    init_steps = size(σ.init_e,2)
    @assert N = size(σ.init_e,1)
    @assert N = size(σ.init_j,1)
    println(Nt)
    xj[1:init_steps] = init_j
    xe[1:init_steps] = init_e
    for i in init_steps+1:Nt
        println(i)
        R = inc[:,i]
        vi = R[end]
        R[end] = 0.0
		fill!(yj,0)
		fill!(ye,0)
		BEAST.ConvolutionOperators.convolve!(yj,Z,xj,csxj,i,jk_start,jk_stop)
        BEAST.ConvolutionOperators.convolve!(ye,Ġ,xe,csxe,i,ek_start,ek_stop)
        iter_max = 100
        for l in 1:iter_max
            if l==1
                if i==1
                    i=2
                end
                xeprev = xe[:,i-1]
                update!(σ, xj, xe, i-1, eq1, eq2)
            else
                xeprev = xe[:,i]
                update!(σ, xj, xe, i, eq1, eq2)
            end
            bσ = BEAST.assemble(σ, eq2.test_space_dict[1].space)
            Q = BEAST.assemble(σop, eq2.test_space_dict[1].space, eq2.trial_space_dict[1].space)
            #Q_inv = inv(Matrix(Q))
            #= println("normQ ", norm(Q))
            println("normbsigma ", norm(bσ)) =#
            #V0 = inv(Z0*G_j0_inv*Q - Ġ0)
            V0 = [Z01 -Ġ01; G_j01 -Q]
            rhs1 = R - yj + ye
            rhs2 = bσ - Q*xeprev
            #b = rhs1 + rhs2         #solve |Z   -Ġ|J =|b|                     
            b = [rhs1; rhs2; 0 ; -vi]
            sol = inv(Matrix(V0))*b
            xj[:,i] = sol[1:N]                            #|G_j  Q|E =|0|
            xe[:,i] = sol[N+1:end]
            println("norm xe ", norm(xe[:,i]-xeprev))
            if norm(xe[:,i]-xeprev) < 1e-5
                break
            end
        end
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
end =#

#marchonintimenl6 integrates the scattering problem with a circuit involving a DC source and a 
#resistor in series
function marchonintimenl6(eq1, eq2,  Z, inc, Ġ, G_j, G_nl, Nt)
    T = eltype(Z)
    Z0 = zeros(T, size(Z)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Z0,Z,1)
    Ġ0 = Ġ.data[1,:,:]
    G_j0 = G_j.data[1,:,:]
    G_nl0 = G_nl.data[1,:,:]
    println(size(G_j0))
    println(size(Ġ0))
    T = eltype(Z0)
    M,N = size(Z0)
    Me,Ne = size(Ġ0)
    Δt = eq1.trial_space_dict[1].time.timestep
    V0 = zeros(N+Ne+2,N+Ne+2)
    padcol1 = zeros(T, N+Ne)
    padcol2 = zeros(T, N+Ne)
    padrow1e = zeros(T, N+Ne)'# this computes the adjoint and not just the transpose, practice caution when complex matrix elements are involved
    padrow2e = zeros(T, N+Ne)'
    padcol2[N] = -1.0
    padrow1e[N] = 1.0
    res= 500.0
    padmat = [1.0 0.0; res/Δt 1.0]
    V0[1:N, 1:N] = Z0
    V0[1:N, N+1:N+Ne] = -Ġ0
    V0[N+1:N+Ne, 1:N] = G_j0
    V0[1:N+Ne, N+Ne+1] = padcol1
    V0[1:N+Ne, N+Ne+2] = padcol2
    V0[N+Ne+1, 1:N+Ne] = padrow1e
    V0[N+Ne+2, 1:N+Ne] = padrow2e
    V0[N+Ne+1:end, N+Ne+1:end] = padmat
    sol = zeros(2*N)
    xj_all = zeros(N)
    xe_all = zeros(N)
    @assert M == size(inc,1)
    xj = zeros(T,N,Nt)
	xe = zeros(T,Ne,Nt)
    xj_cir = zeros(T,2, Nt)
	xeprev = zeros(T,N)
    yj = zeros(T,N)
	ye = zeros(T,N)
    jk_start = 2
    ek_start = 2
    jk_stop = Nt
    ek_stop = BEAST.numfunctions(eq1.trial_space_dict[1].time)+1
    csxj = zeros(T,N,Nt)
    csxe = zeros(T,Ne,Nt)
    σ = eq2.equation.rhs.terms[1].functional
    σop = NonLinearBEM.ConductivityTDop_b(σ)
    bσ = zeros(T, N)
    #G_nl_inv = inv(Matrix(G_nl0))
    #G_j0_inv = inv(Matrix(G_j0))
    println(Nt)
    for i in 1:Nt
        println(i)
        R = inc[:,i]
        vi = eq1.equation.rhs.terms[1].functional.pv1(i*Δt)
        R[end] = 0.0
		fill!(yj,0)
		fill!(ye,0)
		BEAST.ConvolutionOperators.convolve!(yj,Z,xj,csxj,i,jk_start,jk_stop)
        BEAST.ConvolutionOperators.convolve!(ye,Ġ,xe,csxe,i,ek_start,ek_stop)
        iter_max = 200
        for l in 1:iter_max
            if l==1
                if i==1
                    i=2
                end
                xeprev = xe[:,i-1]
                update!(σ, xj, xe, i-1, eq1, eq2)
            else
                xeprev = xe[:,i]
                update!(σ, xj, xe, i, eq1, eq2)
            end
            bσ = BEAST.assemble(σ, eq2.test_space_dict[1].space)
            Q = BEAST.assemble(σop, eq2.test_space_dict[1].space, eq2.trial_space_dict[1].space)
            #Q_inv = inv(Matrix(Q))
            #= println("normQ ", norm(Q))
            println("normbsigma ", norm(bσ)) =#
            #V0 = inv(Z0*G_j0_inv*Q - Ġ0)
            V0[N+1:N+Ne,N+1:N+Ne] = -Q
            rhs1 = R - yj + ye
            rhs2 = bσ - Q*xeprev
            #b = rhs1 + rhs2         #solve |Z   -Ġ|J =|b|                     
            b = [rhs1; rhs2; 0 ; vi+xj_cir[1,i-1]*res/Δt]
            sol = inv(Matrix(V0))*b
            #= invV0 = GMRESSolver(V0, maxiter=1000, restart=0, tol=1e-6)
            mul!(sol, invV0, b) =#
            xj[:,i] = sol[1:N]                            #|G_j  Q|E =|0|
            xe[:,i] = sol[N+1:N+Ne]
            xj_cir[:,i] = sol[N+Ne+1:end]
            #xe_all = hcat(xe_all, xe[:,i])
            #xj_all = hcat(xj_all, xj[:,i])
            println("norm xe ", norm(xe[:,i]-xeprev))
            if norm(xe[:,i]-xeprev) < 1e-4
                break
            end
        end
        if i > 1
            csxj[:,i] .= csxj[:,i-1] .+ xj[:,i]
            csxe[:,i] .= csxe[:,i-1] .+ xe[:,i]
        else
            csxj[:,i] .= xj[:,i]
            csxe[:,i] .= xe[:,i]
        end
        (i % 10 == 0) && print(i, "[", Nt, "] - ")
    end
    return xj, xe, xj_cir, xj_all, xe_all
end

#in marchonintime7, the circuit is complicated to include an inductor and a capacitor to drive the 
#circuit in oscillations
function marchonintimenl7(eq1, eq2,  Z, inc, Ġ, G_j, G_nl, Nt)
    println("tank circuit added")
    T = eltype(Z)
    Z0 = zeros(T, size(Z)[1:2])
    BEAST.ConvolutionOperators.timeslice!(Z0,Z,1)
    Ġ0 = Ġ.data[1,:,:]
    G_j0 = G_j.data[1,:,:]
    G_nl0 = G_nl.data[1,:,:]
    println(size(G_j0))
    println(size(Ġ0))
    T = eltype(Z0)
    M,N = size(Z0)
    Me,Ne = size(Ġ0)
    Δt = eq1.trial_space_dict[1].time.timestep
    V0 = zeros(N+Ne+2,N+Ne+2)
    padcol1 = zeros(T, N+Ne)
    padcol2 = zeros(T, N+Ne)
    padrow1e = zeros(T, N+Ne)'# this computes the adjoint and not just the transpose, practice caution when complex matrix elements are involved
    padrow2e = zeros(T, N+Ne)'
    padcol2[N] = 1.0
    padrow1e[N] = -1.0
    res= 0.001
    L = 0
    C = 0
    padmat = [1.0 -1.0*C; ((L/(Δt)^2)+res/Δt) 1.0]
    V0[1:N, 1:N] = Z0
    V0[1:N, N+1:N+Ne] = -Ġ0
    V0[N+1:N+Ne, 1:N] = G_j0
    V0[1:N+Ne, N+Ne+1] = padcol1
    V0[1:N+Ne, N+Ne+2] = padcol2
    V0[N+Ne+1, 1:N+Ne] = padrow1e
    V0[N+Ne+2, 1:N+Ne] = padrow2e
    V0[N+Ne+1:end, N+Ne+1:end] = padmat
    sol = zeros(2*N)
    xj_all = zeros(N)
    xe_all = zeros(N)
    @assert M == size(inc,1)
    xj = zeros(T,N,Nt)
	xe = zeros(T,Ne,Nt)
    xj_cir = zeros(T,2, Nt)
	xeprev = zeros(T,N)
    yj = zeros(T,N)
	ye = zeros(T,N)
    jk_start = 2
    ek_start = 2
    jk_stop = Nt
    ek_stop = BEAST.numfunctions(eq1.trial_space_dict[1].time)+1
    csxj = zeros(T,N,Nt)
    csxe = zeros(T,Ne,Nt)
    σ = eq2.equation.rhs.terms[1].functional
    σop = NonLinearBEM.ConductivityTDop_b(σ)
    bσ = zeros(T, N)
    #G_nl_inv = inv(Matrix(G_nl0))
    #G_j0_inv = inv(Matrix(G_j0))
    println(Nt)
    for i in 1:Nt
        println(i)
        R = inc[:,i]
        vi = eq1.equation.rhs.terms[1].functional.pv1(i*Δt)
        R[end] = 0.0
        ##do a clean implementation
        #= g = BEAST.creategaussian(2.9504, 2.2128, 0.8)
        R[end] = g(0.01844*i) =#
		fill!(yj,0)
		fill!(ye,0)
		BEAST.ConvolutionOperators.convolve!(yj,Z,xj,csxj,i,jk_start,jk_stop)
        BEAST.ConvolutionOperators.convolve!(ye,Ġ,xe,csxe,i,ek_start,ek_stop)
        iter_max = 200
        for l in 1:iter_max
            if l==1
                if i==1
                    i=2
                end
                xeprev = xe[:,i-1]
                update!(σ, xj, xe, i-1, eq1, eq2)
            else
                xeprev = xe[:,i]
                update!(σ, xj, xe, i, eq1, eq2)
            end
            bσ = BEAST.assemble(σ, eq2.test_space_dict[1].space)
            Q = BEAST.assemble(σop, eq2.test_space_dict[1].space, eq2.trial_space_dict[1].space)
            #Q_inv = inv(Matrix(Q))
            #= println("normQ ", norm(Q))
            println("normbsigma ", norm(bσ)) =#
            #V0 = inv(Z0*G_j0_inv*Q - Ġ0)
            V0[N+1:N+Ne,N+1:N+Ne] = -Q
            rhs1 = R - yj + ye
            rhs2 = bσ - Q*xeprev
            #b = rhs1 + rhs2         #solve |Z   -Ġ|J =|b|
            if i==2
                i_2 = xj_cir[1, i-1]
            else
                i_2 = xj_cir[1, i-2]
            end
            b = [rhs1; rhs2; 0 ; vi+xj_cir[1,i-1]*(res/Δt+(2*L/(Δt)^2)) - i_2*(L/(Δt^2))]
            sol = inv(Matrix(V0))*b
            #= invV0 = GMRESSolver(V0, maxiter=1000, restart=0, tol=1e-6)
            mul!(sol, invV0, b) =#
            xj[:,i] = sol[1:N]                            #|G_j  Q|E =|0|
            xe[:,i] = sol[N+1:N+Ne]
            xj_cir[:,i] = sol[N+Ne+1:end]
            #xe_all = hcat(xe_all, xe[:,i])
            #xj_all = hcat(xj_all, xj[:,i])
            println("norm xe ", norm(xe[:,i]-xeprev))
            if norm(xe[:,i]-xeprev) < 1e-4
                break
            end
        end
        if i > 1
            csxj[:,i] .= csxj[:,i-1] .+ xj[:,i]
            csxe[:,i] .= csxe[:,i-1] .+ xe[:,i]
        else
            csxj[:,i] .= xj[:,i]
            csxe[:,i] .= xe[:,i]
        end
        (i % 10 == 0) && print(i, "[", Nt, "] - ")
    end
    return xj, xe, xj_cir, xj_all, xe_all
end