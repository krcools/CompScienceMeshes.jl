# Neighborhoods

The Neighborhood concept represents the differential geometric notion of the derivative of
a chart in a point of its parameter domain. Not only can it be queried for the parametric and
Cartesian coordinates of the input and output point respectively, it also contains information
on tangents and Jacobian determinants. For the special case where the dimension of the range
manifold is one less than that of the surrounding space, access is provided to the unit normal.

This concept is an enriched point concept and should allow the construction of the majority of
kernels encountered in finite and boundary element methods.

```@docs
parametric
cartesian
jacobian
tangents
normal
```

In addition to this interface, a neighborhood has appropriate methods for `getindex` defined, so
that it can be used with any function that expects a tuple of coordinates. In other words, a
neighborhood is a model for Point.

## Numerical Quadrature

This package provides a number of routines that aim to facilitate numerical quadrature over
charts. In addition, it includes a default set of quadrature rules for segments (1D) and triangles
(2D).

```@docs
quadpoints
```
