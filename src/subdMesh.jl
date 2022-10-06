export SEdge,SElement,SVertex,subdMesh

mutable struct SElement
    RingNodes::Vector{Int64}
    Edges::Vector{Int64}
end

mutable struct SEdge
    EndPoints::Vector{Int64}
    RelatedPoints::Vector{Int64}
    TwoElems::Vector{Int64}
    orientation::Vector{Int64}
end

mutable struct SVertex{T}
    Elements::Vector{Int64}
    Valence::Int64
    Edges::Vector{Int64}
    Coords::SVector{3,T}
end

mutable struct subdMesh
    edges::Vector{SEdge}
    vertices::Vector{SVertex}
    elements::Vector{SElement}
    mesh::Mesh
end

function coordtype(m::subdMesh) coordtype(m.mesh) end

"""
    GSubdMesh(mesh::Mesh)

    Given a linear polygon mesh, construct a data structure for subdivision surfaces
"""
function GSubdMesh(mesh::Mesh)
    D1 = 3
    edges = skeleton(mesh,1)
    faces = skeleton(mesh,2)
    connectMat = connectivity(edges,faces)
    NV, NE, NF = numvertices(mesh), numcells(edges), numcells(faces)
    # nv = NV + NE
    Svertices = Vector{SVertex}(undef,NV)
    Sedges = Vector{SEdge}(undef,NE)
    Selements = Vector{SElement}(undef,NF)
    for V in 1 : numvertices(mesh)
        ElementIndices = Int[]
        EdgeIndices = Int[]
        for (F, Face) in enumerate(cells(faces))
            for iv = 1 : 3
                # if V is a vertex of the face add all the vertices to the mask
                if Face[iv] == V
                    append!(ElementIndices, F)
                end
            end
        end
        for (E, Edge) in enumerate(cells(edges))
            for ie = 1 : 2
                if Edge[ie] == V
                    append!(EdgeIndices,E)
                end
            end
        end
        coords = mesh.vertices[V]
        # delete duplicate vertices
        # RingverticesIndices = unique(RingverticesIndices)
        # calculate valence for vertex V
        valence = length(ElementIndices)
        SV=SVertex(ElementIndices,valence,EdgeIndices,coords)
        # append!(Svertices, SV)
        Svertices[V] = SV
    end
    for (E, Edge) in enumerate(cells(edges))
        EndVerticesIndices = Edge
        FacesIndex = []
        RelatedVerticesIndices = []
        orientation = []
        for (F, Face) in enumerate(cells(faces))
            if connectMat[F,E] != 0
                    # Find out two faces in the mask
                append!(FacesIndex, F)
                append!(orientation,connectMat[F,E])
                    ## Loop over vertices in the faces find out related vertices
                for iv = 1 : length(Face)
                    if Face[iv] != Edge[1] && Face[iv] != Edge[2]
                            append!(RelatedVerticesIndices, Face[iv])
                    end
                end
            end
        end
        SE=SEdge(EndVerticesIndices,RelatedVerticesIndices,FacesIndex,orientation)
        # append!(Sedges,SE)
        Sedges[E] = SE
    end
    for (F, Face) in enumerate(cells(faces))
        VertexIndices = Face
        irreg_count = 0
        irreg_index = 0
        irreg_valence = 0
        EdgesIndices,orientation = find_edges(F,Face,edges,connectMat)
        for iv = 1 : length(Face)
            ivalence = Svertices[VertexIndices[iv]].Valence
            if ivalence != 6
                irreg_count += 1
                irreg_index = iv
                irreg_valence = ivalence
            end
        end
        if irreg_count > 1
            error("An element has more than 1 extraordinary point, please refine mesh")
        elseif irreg_count == 0
            neighborFaces,Nn = find_neighbor(faces,edges,F,EdgesIndices,Sedges)
            neighborFaces2,Mn = find_neighbor2(neighborFaces,faces,edges,VertexIndices,EdgesIndices,Sedges,connectMat)
            # reorder nodes
            # Ringnodes=[Mn[1],Mn[6],Nn[1],VertexIndices[1],Nn[3],Mn[2],VertexIndices[2],VertexIndices[3],Mn[5],Mn[3],Nn[2],Mn[4]]
            Ringnodes=[VertexIndices[1],VertexIndices[2],Nn[1],Mn[1],Mn[6],Nn[3],VertexIndices[3],Nn[2],Mn[3],Mn[2],Mn[4],Mn[5]]
            SF=SElement(Ringnodes,EdgesIndices)
        else
            # print("element $(F) has 1 extraordinary point!")
            Ringnodes = Vector{Int64}(undef,irreg_valence+6)
            Ringnodes[1] = VertexIndices[irreg_index]
            next = [2,3,1]
            next2 = [3,1,2]
            Ringnodes[2] = VertexIndices[next[irreg_index]]
            Ringnodes[irreg_valence+1] = VertexIndices[next2[irreg_index]]
            # if irreg_index == 1
                # Ringnodes[2] = VertexIndices[irreg_index+1]
                # Ringnodes[irreg_valence+1] = VertexIndices[irreg_index+2]
            # elseif irreg_index == 2
            #     Ringnodes[2] = VertexIndices[irreg_index+1]
            #     Ringnodes[irreg_valence+1] = VertexIndices[irreg_index-1]
            # elseif irreg_index == 3
            #     Ringnodes[2] = VertexIndices[irreg_index-2]
            #     Ringnodes[irreg_valence+1] = VertexIndices[irreg_index-1]
            # end
            #start a loop to find the ring of nodes 1:valence-1
            firstEdge = Sedges[EdgesIndices[irreg_index]]
            firstF = F
            for ir = 1:irreg_valence - 2
                twoelement =firstEdge.TwoElems
                tworelatedpoints = firstEdge.RelatedPoints
                neighborF = 0
                neighborV = 0
                if twoelement[1] == firstF
                    neighborF = twoelement[2]
                    neighborV = tworelatedpoints[2]
                elseif twoelement[2] == firstF
                    neighborF = twoelement[1]
                    neighborV = tworelatedpoints[1]
                end
                Ringnodes[ir+2] = neighborV
                neighborface = cells(faces)[neighborF]
                firstF = neighborF
                edgesIndices_n,orientation_n = find_edges(neighborF,neighborface,edges,connectMat)
                for E in edgesIndices_n
                    Edge = cells(edges)[E]
                    if (Edge[1] == Ringnodes[1] && Edge[2] == neighborV) || (Edge[2] == Ringnodes[1] && Edge[1] == neighborV)
                        firstEdge = Sedges[E]
                    end
                end
            end
            secondEdge = Sedges[EdgesIndices[next[irreg_index]]]
            firstF = F
            for ir = 1:3
                twoelement =secondEdge.TwoElems
                tworelatedpoints = secondEdge.RelatedPoints
                neighborF = 0
                neighborV = 0
                if twoelement[1] == firstF
                    neighborF = twoelement[2]
                    neighborV = tworelatedpoints[2]
                elseif twoelement[2] == firstF
                    neighborF = twoelement[1]
                    neighborV = tworelatedpoints[1]
                end
                Ringnodes[irreg_valence + 1 + ir] = neighborV
                neighborface = cells(faces)[neighborF]
                firstF = neighborF
                edgesIndices_n,orientation_n = find_edges(neighborF,neighborface,edges,connectMat)
                for E in edgesIndices_n
                    Edge = cells(edges)[E]
                    if (Edge[1] == Ringnodes[2] && Edge[2] == neighborV)||(Edge[2] == Ringnodes[2] && Edge[1] == neighborV)
                        secondEdge = Sedges[E]
                    end
                end
            end

            thirdEdge = Sedges[EdgesIndices[next2[irreg_index]]]
            firstF = F
            for ir = 1:3
                twoelement =thirdEdge.TwoElems
                tworelatedpoints = thirdEdge.RelatedPoints
                neighborF = 0
                neighborV = 0
                if twoelement[1] == firstF
                    neighborF = twoelement[2]
                    neighborV = tworelatedpoints[2]
                elseif twoelement[2] == firstF
                    neighborF = twoelement[1]
                    neighborV = tworelatedpoints[1]
                end
                if ir == 2
                    Ringnodes[irreg_valence + 6] = neighborV
                elseif ir ==3
                    Ringnodes[irreg_valence + 5] = neighborV
                end

                neighborface = cells(faces)[neighborF]
                firstF = neighborF
                edgesIndices_n,orientation_n = find_edges(neighborF,neighborface,edges,connectMat)
                for E in edgesIndices_n
                    Edge = cells(edges)[E]
                    if (Edge[1] == Ringnodes[irreg_valence+1] && Edge[2] == neighborV)||(Edge[2] == Ringnodes[irreg_valence+1] && Edge[1] == neighborV)
                        thirdEdge = Sedges[E]
                    end
                end
            end
            SF=SElement(Ringnodes,EdgesIndices)
        end
        # append!(Selements,SF)
        Selements[F] = SF
    end
    subdmesh=subdMesh(Sedges,Svertices,Selements,mesh)
end
"""
    find_neighbor(faces,edges,F,EdgesIndices,orientation,Sedges)

    Given a face find out the neighbor elements (share edge) and vertices.
"""
function find_neighbor(faces,edges,F,EdgesIndices,Sedges)
    neighborFaces = Vector{Int64}(undef,3)
    neighborVertices = Vector{Int64}(undef,3)
    for iE = 1:3
        E = EdgesIndices[iE]
        twoelement = Sedges[E].TwoElems
        tworelatedpoints = Sedges[E].RelatedPoints
        if twoelement[1] == F
            neighborFaces[iE] = twoelement[2]
            neighborVertices[iE] = tworelatedpoints[2]
        elseif twoelement[2] == F
            neighborFaces[iE] = twoelement[1]
            neighborVertices[iE] = tworelatedpoints[1]
        end
    end
    neighborFaces,neighborVertices
end

function find_neighbor2(neighborFaces,faces,edges,VertexIndices,EdgesIndices,Sedges,connectMat)
    neighborFaces2 = Vector{Int64}(undef,6)
    neighborVertices2 = Vector{Int64}(undef,6)
    C=[1,2,2,3,3,1]
    for i = 1 : 3
        F = neighborFaces[i]
        Face = cells(faces)[F]
        edgesIndices_n,orientation = find_edges(F,Face,edges,connectMat)
        for E in edgesIndices_n
            Edge = cells(edges)[E]
            if E != EdgesIndices[i] && (Edge[1] == VertexIndices[C[2*i-1]] || Edge[2] == VertexIndices[C[2*i-1]])
                twoelement = Sedges[E].TwoElems
                tworelatedpoints = Sedges[E].RelatedPoints
                if twoelement[1] == F
                    neighborFaces2[2*i-1] = twoelement[2]
                    neighborVertices2[2*i-1] = tworelatedpoints[2]
                elseif twoelement[2] == F
                    neighborFaces2[2*i-1] = twoelement[1]
                    neighborVertices2[2*i-1] = tworelatedpoints[1]
                end
            elseif E != EdgesIndices[i] && (Edge[1] == VertexIndices[C[2*i]] || Edge[2] == VertexIndices[C[2*i]])
                twoelement = Sedges[E].TwoElems
                tworelatedpoints = Sedges[E].RelatedPoints
                if twoelement[1] == F
                    neighborFaces2[2*i] = twoelement[2]
                    neighborVertices2[2*i] = tworelatedpoints[2]
                elseif twoelement[2] == F
                    neighborFaces2[2*i] = twoelement[1]
                    neighborVertices2[2*i] = tworelatedpoints[1]
                end
            end
        end
    end
    neighborFaces2,neighborVertices2
end

function find_edges(F,Face,edges,connectMat)
    EdgesIndices = Vector{Int64}(undef,3)
    orientation = Vector{Int64}(undef,3)
    for (E, Edge) in enumerate(cells(edges))
        if connectMat[F,E] == 1
            if Edge[1] == Face[1]||Edge[2] == Face[2]
                EdgesIndices[1] = E
                orientation[1] = 1
            elseif Edge[1] == Face[2]||Edge[2] == Face[3]
                EdgesIndices[2] = E
                orientation[1] = 1
            elseif Edge[1] == Face[3]||Edge[2] == Face[1]
                EdgesIndices[3] = E
                orientation[1] = 1
            end
        elseif connectMat[F,E] == -1
            if Edge[1] == Face[2]||Edge[2] == Face[1]
                EdgesIndices[1] = E
                orientation[1] = -1
            elseif Edge[1] == Face[3]||Edge[2] == Face[2]
                EdgesIndices[2] = E
                orientation[1] = -1
            elseif Edge[1] == Face[1]||Edge[2] == Face[3]
                EdgesIndices[3] = E
                orientation[1] = -1
            end
        end
    end
    EdgesIndices,orientation
end
