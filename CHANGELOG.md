# Changelog

## Version 0.10

- New method for union that returns the input meshes as submeshes of the union mesh.
- `parent` back-edge for submeshes.
- embedding implemented for surface meshes.
- Minimum supported Julia version is 1.10

## Version 0.9.1

- `meshcuboid` returns outward oriented meshes for both the `gmsh` and `compsiencemeshes` generators.

## Version 0.9.3

- The `subcharts` API is introduced, which returns an iterable container of tuples comprising the chart of a subentity such as a face or an edge, together with the overlap map from the domain of the subchart to the domain of the original chart. The APi is introduced to enable efficient implementations of local interpolation routines for finite element spaces.
