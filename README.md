# CompScienceMeshes

Geometry types and algorithms for computational science

[![Build status](https://github.com/krcools/CompScienceMeshes.jl/workflows/CI/badge.svg)](https://github.com/krcools/CompScienceMeshes.jl/actions)
[![codecov](https://codecov.io/gh/krcools/CompScienceMeshes.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/krcools/CompScienceMeshes.jl)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://krcools.github.io/CompScienceMeshes.jl/dev/)

## Installation

In addition to the dependencies declared in REQUIRE, this package relies for some of its functionality on `gmsh`. Make sure `gmsh` is installed and on the system path if you require this functionality.

## New in Version 0.4.0

From 0.4.0, a `Mesh` is an iterable object whose values are opaque pointers to geometric elements. Before, they were tuples of indices into the vertex buffer. The advantage is that starting from such a pointer, it is possible to traverse non-trivial Mesh data structures. The motivating example is that of a refinement-parent pair of Meshes. Given a pointer into the parent mesh, the new approach allows to iterate over the children in the refinement. With the deprecated approach this was simply not possible; the for example triple of vertex indices defining a triangle gives no access to what elements in the refinement were created by refining the parent triangle.

Some important changes for users of CompScienceMeshes:

- The mesh itself is an iterator over pointers. If an iterator of tuples of vertex indices is required, use `cells(mesh)`.
- Charts are constructed from element pointers: `chart(mesh, pointer_to_element)`, as oppposed to `chart(mesh, tuple_of_vertex_indices)`.
- Predicates passeed to `submesh` now require to arguments: the mesh *and* a pointer to one of its elements. It is the users responsibility to make sure that the mesh and the pointer match up.


## Introduction

This package provides the geometric framework to facilitated the construction of finite element spaces and the assembly of matrices stemming from the discretisation of local (differential) and global (integral) operators on those finite element spaces.

The package roughly contains three components:

* Meshes: allowing for the (almost) linear construction of connectivity matrices. A default implementation is provided but the algorithms should be easily extendable to user defined mesh structures. It is very common, for example, that mesh data structures contain problem specific information (local elasticity, permittivity, boundary conditions). User can use those enriched structures if they extend a limited number of functions.

* Charts: a concept designed after the differential geometric concept of a chart on a manifold. It allows for the construction of points in Euclidian space from a set of parameters and the other way around.

* Neighborhoods: a concept designed after the derivative of a chart as a map from the parametrising vector space to the tangent space of a point of the manifold. It allows querying for tangents, normal, and the Jacobian determinant for use in integration routines.
