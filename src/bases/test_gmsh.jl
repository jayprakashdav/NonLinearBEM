using CompScienceMeshes


"""
find the cells over which the basis function should be defined for a given internal 
vertex in the port
"""

function findcell(m0,m1)
    cps = CompScienceMeshes.cellpairs(m0,m1)
    eps = edgepairs(m1,skeleton(m1,0))
    for i in axes(eps,2)
        edg1 = m1.faces[eps[1,i]]
        edg2 = m1.faces[eps[2,i]]
        cell1 = cps[1,eps[1,i]]
        cell2 = cps[1,eps[2,i]]
        println(edg1)
        println(edg2)
        println(m0.faces[cell1])
        println(m0.faces[cell2])
        println("           ")
    end
end 

"""
Obtained the edge pairs over which the varying RWG basis functions should be defined
the associated cells to those edges are also obtained in cps1 and cps2 for m1 and m2
respectively.
Thus, have all the geometric information ready to create the varying and global RWG basis
functions.
"""

"""
Given a mesh of an edge and set of vertices from that edge mesh, 'cellpairs' will
generate a 2 x K matrix, where K is the number of pairs and each column contains
a pair of indices in the cell array of 'mesh' that have one of the supplied vertices in
common 
"""

"""
position in the basis function struct definition and 
"""

function edgepairs(mesh, vertices; dropjunctionpair=false)

    ndrops = dropjunctionpair ? 1 : 0
    @assert dimension(vertices)+1 == dimension(mesh)

    numedges = numcells(vertices)

    # perform a dry run to determine the number of cellpairs
    v2e, nn = vertextocellmap(mesh)
    c = findall(nn.==2)
    k = length(c)
    edgpairs = zeros(Int, 2, k)

    for i = 1:k
        edgpairs[1,i] = v2e[c[i],1]
        edgpairs[2,i] = v2e[c[i],2]
    end
    edgpairs
end