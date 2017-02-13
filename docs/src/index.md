# CompScienceMeshes.jl

Mesh querying and differential geometry for finite and boundary element methods.

## Introduction

This package provides the geometric framework to facilitated the construction of finite element spaces and the assembly of matrices stemming from the discretisation of local (differential) and global (integral) operators on those finite element spaces.

The package roughly contains three components:

* Meshing algorithm: allowing for the (almost) linear construction of connectivity matrices. A default implementation is provided but the algorithms should be easily extendable to user defined mesh structures. It is very common, for example, that mesh data structures contain problem specific information (local elasticity, permittivity, boundary conditions). User can use those enriched structures if they extend a limited number of functions.
* Charts: a concept designed after the differential geometric concept of a chart on a manifold. It allows for the construction of points in Euclidian space from a set of parameters and the other way around.
* Neighborhoods: a concept designed after the derivative of a chart as a map from the parametrising vector space to the tangent space of a point of the manifold. It allows querying for tangents, normal, and the Jacobian determinant for use in integration routines.
