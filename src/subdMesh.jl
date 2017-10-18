export SEdge,SElement,SVertex,subdMesh,shapefun_reg,shapefun_reg_der,shapefun_reg_der2

type SElement
    # RingElements::Vector{Int64}
    RingNodes::Vector{Int64}
    Edges::Vector{Int64}
end

type SEdge
    EndPoints::Vector{Int64}
    RelatedPoints::Vector{Int64}
    TwoElems::Vector{Int64}
    orientation::Vector{Int64}
end

type SVertex
    Elements::Vector{Int64}
    Valence::Int64
    Edges::Vector{Int64}
end

type subdMesh
    edges::Vector{SEdge}
    vertices::Vector{SVertex}
    elements::Vector{SElement}
end
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
    Svertices = Vector{SVertex}(NV)
    Sedges = Vector{SEdge}(NE)
    Selements = Vector{SElement}(NF)
    for V in 1 : numvertices(mesh)
        ElementIndices = []
        EdgeIndices = []
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
        # delete duplicate vertices
        # RingverticesIndices = unique(RingverticesIndices)
        # calculate valence for vertex V
        valence = length(ElementIndices)
        SV=SVertex(ElementIndices,valence,EdgeIndices)
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
            Ringnodes = Vector{Int64}(irreg_valence+6)
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
        print(Ringnodes)
        print("\n")
    end
    mesh=subdMesh(Sedges,Svertices,Selements)
end
"""
    find_neighbor(faces,edges,F,EdgesIndices,orientation,Sedges)

    Given a face find out the neighbor elements (share edge) and vertices.
"""
function find_neighbor(faces,edges,F,EdgesIndices,Sedges)
    neighborFaces = Vector{Int64}(3)
    neighborVertices = Vector{Int64}(3)
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
    neighborFaces2 = Vector{Int64}(6)
    neighborVertices2 = Vector{Int64}(6)
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
    EdgesIndices = Vector{Int64}(3)
    orientation = Vector{Int64}(3)
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

function shapefun_reg(v,w)
    u=1.0-v-w;

    x = zeros(1,12)

    x[4] = (u*u*u*u + 2.0*u*u*u*v)/12.0

    x[5] = (u*u*u*u + 2.0*u*u*u*w)/12.0

    x[3] = (u*u*u*u + 2.0*u*u*u*w + 6.0*u*u*u*v + 6.0*u*u*v*w + 12.0*u*u*v*v + 6.0*u*v*v*w + 6.0*u*v*v*v + 2.0*v*v*v*w + v*v*v*v) / 12.0

    x[1] = (6.0*u*u*u*u + 24.0*u*u*u*w + 24.0*u*u*w*w + 8.0*u*w*w*w + w*w*w*w + 24.0*u*u*u*v + 60.0*u*u*v*w + 36.0*u*v*w*w + 6.0*v*w*w*w + 24.0*u*u*v*v + 36.0*u*v*v*w + 12.0*v*v*w*w + 8.0*u*v*v*v + 6.0*v*v*v*w + v*v*v*v) / 12.0

    x[6] = (u*u*u*u + 6.0*u*u*u*w + 12.0*u*u*w*w + 6.0*u*w*w*w + w*w*w*w + 2.0*u*u*u*v + 6.0*u*u*v*w + 6.0*u*v*w*w + 2.0*v*w*w*w) / 12.0

    x[10] = (2.0*u*v*v*v + v*v*v*v) / 12.0

    x[2] = (u*u*u*u + 6.0*u*u*u*w + 12.0*u*u*w*w + 6.0*u*w*w*w + w*w*w*w + 8.0*u*u*u*v + 36.0*u*u*v*w + 36.0*u*v*w*w + 8.0*v*w*w*w + 24.0*u*u*v*v + 60.0*u*v*v*w + 24.0*v*v*w*w + 24.0*u*v*v*v + 24.0*v*v*v*w + 6.0*v*v*v*v) / 12.0

    x[7] = (u*u*u*u + 8.0*u*u*u*w + 24.0*u*u*w*w + 24.0*u*w*w*w + 6.0*w*w*w*w + 6.0*u*u*u*v + 36.0*u*u*v*w + 60.0*u*v*w*w + 24.0*v*w*w*w + 12.0*u*u*v*v + 36.0*u*v*v*w + 24.0*v*v*w*w + 6.0*u*v*v*v + 8.0*v*v*v*w + v*v*v*v)/12.0

    x[12] = (2.0*u*w*w*w + w*w*w*w) / 12.0

    x[9] = (2.0*v*v*v*w + v*v*v*v) / 12.0

    x[8] = (2.0*u*w*w*w + w*w*w*w + 6.0*u*v*w*w + 6.0*v*w*w*w + 6.0*u*v*v*w + 12.0*v*v*w*w + 2.0*u*v*v*v + 6.0*v*v*v*w + v*v*v*v) / 12.0

    x[11] = (w*w*w*w + 2.0*v*w*w*w) / 12.0

    return x
end

function shapefun_reg_der(v,w)
    u=1.0-v-w
    der1=zeros(2,12)
    der1[1,4] = (-6.0*v*pow(u,2.0) - 2.0*pow(u,3.0))/12.0
    der1[2,4] = (-6.0*v*pow(u,2.0) - 4.0*pow(u,3.0))/12.0

    der1[1,5] = (-4.0*pow(u,3.0)-6.0*pow(u,2.0)*w)/12.0
    der1[2,5] = (-2.0*pow(u,3.0)-6.0*pow(u,2.0)*w)/12.0

    der1[1,3] = (-2.0*pow(v,3.0)-6.0*pow(v,2.0)*u + 6.0*v*pow(u,2.0)+2.0*pow(u,3.0))/12.0
    der1[2,3] = (-4.0*pow(v,3.0)-18.0*pow(v,2.0)*u - 12.0*v*pow(u,2.0)-2.0*pow(u,3.0) - 6.0*pow(v,2.0)*w-12.0*v*u*w - 6.0*pow(u,2.0)*w)/12.0

    der1[1,1] = (-4.0*pow(v,3.0)-24.0*pow(v,2.0)*u - 24.0*v*pow(u,2.0)-18.0*pow(v,2.0)*w - 48.0*v*u*w-12.0*pow(u,2.0)*w - 12.0*v*pow(w,2.0) - 12.0*u*pow(w,2.0) - 2.0*pow(w,3.0))/12.0
    der1[2,1] = (-2.0*pow(v,3.0)-12.0*pow(v,2.0)*u - 12.0*v*pow(u,2.0)-12.0*pow(v,2.0)*w - 48.0*v*u*w-24.0*pow(u,2.0)*w - 18.0*v*pow(w,2.0)-24.0*u*pow(w,2.0) - 4.0*pow(w,3.0))/12.0

    der1[1,6] = (-6.0*v*pow(u,2.0)-2.0*pow(u,3.0) - 12.0*v*u*w-12.0*pow(u,2.0)*w - 6.0*v*pow(w,2.0)-18.0*u*pow(w,2.0) - 4.0*pow(w,3.0))/12.0

    der1[2,6] = (2.0*pow(u,3.0)+6.0*pow(u,2.0)*w - 6.0*u*pow(w,2.0)-2.0*pow(w,3.0))/12.0

    der1[1,10] = (2.0*pow(v,3.0)+6.0*pow(v,2.0)*u)/12.0
    der1[2,10] = -pow(v,3.0)/6.0

    der1[1,2] = (24.0*pow(v,2.0)*u+24.0*v*pow(u,2.0) + 4.0*pow(u,3.0)+12.0*pow(v,2.0)*w + 48.0*v*u*w+18.0*pow(u,2.0)*w + 12.0*v*pow(w,2.0)+12.0*u*pow(w,2.0) + 2.0*pow(w,3.0))/12.0
    der1[2,2] = (12.0*pow(v,2.0)*u+12.0*v*pow(u,2.0) + 2.0*pow(u,3.0)-12.0*pow(v,2.0)*w + 6.0*pow(u,2.0)*w-12.0*v*pow(w,2.0) - 6.0*u*pow(w,2.0)-2.0*pow(w,3.0))/12.0

    der1[1,7] = (-2.0*pow(v,3.0)-6.0*pow(v,2.0)*u + 6.0*v*pow(u,2.0)+2.0*pow(u,3.0) - 12.0*pow(v,2.0)*w+12.0*pow(u,2.0)*w - 12.0*v*pow(w,2.0)+12.0*(1.0-v-w)*pow(w,2.0))/12.0
    der1[2,7] = (2.0*pow(v,3.0)+12.0*pow(v,2.0)*u + 18.0*v*pow(u,2.0)+4.0*pow(u,3.0) + 12.0*pow(v,2.0)*w+48.0*v*u*w + 24.0*pow(u,2.0)*w+12.0*v*pow(w,2.0) + 24.0*u*pow(w,2.0))/12.0

    der1[1,12] = -pow(w,3.0)/6.0
    der1[2,12] = (6.0*u*pow(w,2.0)+2.0*pow(w,3.0))/12.0

    der1[1,9] = (4.0*pow(v,3.0)+6.0*pow(v,2.0)*w)/12.0
    der1[2,9] = pow(v,3.0)/6.0

    der1[1,8]= (2.0*pow(v,3.0)+6.0*pow(v,2.0)*u + 12.0*pow(v,2.0)*w+12.0*v*u*w + 18.0*v*pow(w,2.0)+6.0*u*pow(w,2.0) + 4.0*pow(w,3.0))/12.0

    der1[2,8]= (4.0*pow(v,3.0)+6.0*pow(v,2.0)*u + 18.0*pow(v,2.0)*w+12.0*v*u*w + 12.0*v*pow(w,2.0)+6.0*u*pow(w,2.0) + 2.0*pow(w,3.0))/12.0

    der1[1,11] = pow(w,3.0)/6.0
    der1[2,11] = (6.0*v*pow(w,2.0)+4.0*pow(w,3.0))/12.0

    return der1

end

function shapefun_reg_der2(v,w)
    u = 1.0-v-w
    ders=zeros(3,12)
    der2[1,4] = v*u
    der2[2,4] = v*u+pow(u,2.0)
    der2[3,4] = (12.0*v*u+6.0*pow(u,2.0))/12.0

    der2[1,5] = pow(u,2.0)+u*w
    der2[2,5] = u*w
    der2[3,5] = (6.0*pow(u,2.0)+12.0*u*w)/12.0

    der2[1,3] = -2.0*v*u
    der2[2,3] = pow(v,2.0)+v*u+v*w+u*w
    der2[3,3] = (6.0*pow(v,2.0)-12.0*v*u -6.0*pow(u,2.0))/12.0

    der2[1,1] = pow(v,2.0)-2.0*pow(u,2.0) + v*w-2.0*u*w
    der2[2,1] = -2.0*v*u-2.0*pow(u,2.0) + v*w+pow(w,2.0)
    der2[3,1] = (6.0*pow(v,2.0)-12.0*pow(u,2.0) + 24.0*v*w+6.0*pow(w,2.0))/12.0

    der2[1,6] = v*u + v*w + u*w + pow(w,2.0)
    der2[2,6] = - 2.0*u*w
    der2[3,6] = (- 6.0*pow(u,2.0) - 12.0*u*w + 6.0*pow(w,2.0))/12.0

    der2[1,10] = v*u
    der2[2,10] = 0.0
    der2[3,10] = -pow(v,2.0)/2.0

    der2[1,2] = (-24.0*pow(v,2.0)+12.0*pow(u,2.0)-24.0*v*w + 12.0*u*w)/12.0
    der2[2,2] = (-24.0*pow(v,2.0)-24.0*v*u-24.0*v*w - 24.0*u*w)/12.0
    der2[3,2] = (-12.0*pow(v,2.0)+6.0*pow(u,2.0)-24.0*v*w - 12.0*u*w-6.0*pow(w,2.0))/12.0

    der2[1,7] = -2.0*v*u-2.0*v*w-2.0*u*w - 2.0*pow(w,2.0)
    der2[2,7] = v*u+pow(u,2.0)-2.0*v*w - 2.0*pow(w,2.0)
    der2[3,7] = (-6.0*pow(v,2.0)-12.0*v*u+6.0*pow(u,2.0) - 24.0*v*w-12.0*pow(w,2.0))/12.0

    der2[1,12] = 0.0
    der2[2,12] = u*w
    der2[3,12] = -pow(w,2.0)/2.0

    der2[1,9] = (12.0*pow(v,2.0)+12.0*v*w)/12.0
    der2[2,9] = 0.0
    der2[3,9] = pow(v,2.0)/2.0

    der2[1,8]= (12.0*v*u+12.0*v*w+12.0*u*w + 12.0*pow(w,2.0))/12.0
    der2[2,8]= pow(v,2.0)+v*u+v*w+u*w
    der2[3,8]= (6.0*pow(v,2.0)+12.0*v*u+24.0*v*w + 12.0*u*w+6.0*pow(w,2.0))/12.0

    der2[1,11]= 0.0
    der2[2,11]= v*w+pow(w,2.0)
    der2[3,11]= pow(w,2.0)/2.0

  return der2;

end

function pow(x,y)
    return x^y
end
