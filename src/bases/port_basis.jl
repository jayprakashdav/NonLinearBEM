function raviartthomaswithport(mesh::CompScienceMeshes.AbstractMesh{U,D1,T}, portmesh::CompScienceMeshes.AbstractMesh{U,D2,T},direction) where {U,D1,D2,T}

    #determine the number of varying RWGs for port edges
    eps = nlmodelling.edgepairs(portmesh,skeleton(portmesh,0))
    cps = CompScienceMeshes.cellpairs(mesh, portmesh) 
    numpairs = size(eps,2)
    numportedges = size(portmesh.faces,1)  
    functions = Vector{Vector{BEAST.Shape{T}}}(undef,numpairs+1)
    positions = Vector{vertextype(mesh)}(undef,numpairs+1)
    Cells = cells(mesh)
    for i in axes(eps,2)
        e1 = portmesh.faces[eps[1,i]]
        e2 = portmesh.faces[eps[2,i]]
        c1 = cps[1,eps[1,i]]
        c2 = cps[1,eps[2,i]]
        oppv1 = 1
        oppv2 = 1
        for vert = 1:3
            if !((e1[1]==Cells[c1][vert]) | (e1[2]==Cells[c1][vert]))
                oppv1 = vert
                break
            end
        end
        for vert = 1:3
            if !((e2[1]==Cells[c2][vert]) | (e2[2]==Cells[c2][vert]))
                oppv2 = vert
                break
            end
        end
        functions[i] = [
            BEAST.Shape{T}(c1, abs(oppv1), T(+1.0)),
            BEAST.Shape{T}(c2, abs(oppv2), T(-1.0))]
        ctr1 = cartesian(center(chart(mesh, c1)))
        ctr2 = cartesian(center(chart(mesh, c2)))
        positions[i] = (ctr1 + ctr2) / 2
    end

    #define the global constant RWG for the port edges
    functions[numpairs+1] = []
    for i=numpairs+1
        temppos = []
        for j in axes(cps,2)
            e1 = portmesh.faces[j]
            c1 = cps[1,j]
            oppv1 = 1
            for vert = 1:3
                if !((e1[1]==Cells[c1][vert]) | (e1[2]==Cells[c1][vert]))
                    oppv1 = vert
                    break
                end
            end
            push!(functions[i],BEAST.Shape{T}(c1, abs(oppv1), T(direction*(0.05))))
            if j==1
                temppos = cartesian(center(chart(mesh,c1)))
            else
                temppos = temppos+cartesian(center(chart(mesh,c1)))
            end
        end
        positions[i] = temppos / size(cps,2)
        #positions[i] = sum(globalports.vertices)/length(globalports.vertices)
    end

    geo = mesh
    RTBasis(geo, functions, positions)
end

function raviartthomaswithport(dev_mesh::CompScienceMeshes.AbstractMesh{U,D1,T}, port1mesh::CompScienceMeshes.AbstractMesh{U,D2,T}, port2mesh::CompScienceMeshes.AbstractMesh{U,D3,T}) where{U,D1,D2,D3,T}
    outward = 1.0
    inward = -1.0
    X = raviartthomas(dev_mesh)
    Y = raviartthomaswithport(dev_mesh, port1mesh, outward)
    Z = raviartthomaswithport(dev_mesh, port2mesh, inward)
    
    totalfns = Vector{Vector{BEAST.Shape{T}}}(undef,0)
    append!(totalfns,X.fns)
    append!(totalfns,Y.fns[1:end-1])
    append!(totalfns,Z.fns[1:end-1])
    globalfn = append!(Y.fns[end],Z.fns[end])
    append!(totalfns,[globalfn])


    totalpos = Vector{vertextype(dev_mesh)}(undef,0)
    append!(totalpos,X.pos)
    append!(totalpos,Y.pos[1:end-1])
    append!(totalpos,Z.pos[1:end-1])
    globalfnpos = (Y.pos[end]+Z.pos[end])/2
    append!(totalpos,[globalfnpos])

    ids_internals = fill(1, numfunctions(X))
    ids_emitters = fill(2, numfunctions(Y)-1)
    ids_collectors = fill(3, numfunctions(Z)-1)
    ids_global = fill(4,1)

    return RTBasis(X.geo,totalfns,totalpos), vcat(ids_internals,ids_emitters,ids_collectors,ids_global)
end