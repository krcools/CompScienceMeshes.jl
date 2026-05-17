# Version 0.12 

- Fixed bug where submeshes were sometimes treated as refinements.
- [BREAKING]: PlotlyJS support is now implemented as extension. Functions `patch`, `wireframe`, and `normals` are retired in favour of standard PlolyJS language `scatter3d`, `mesh3d` and `cone`.

# Version 0.11

- [BREAKING]: `cells(mesh)` returns no longer tuples of integers that are indices into the mesh vertex buffer, but typed objects that should be interpreted as *local graphs*, and that can be querried for adjacency and orientation information. This is done to avoid ambiguity when considering higher order meshes and to allow core algorithms such as `connectivity` and `skeleton` to be written generically for both lowest order and higher order meshes.
- [BREAKING]: PlotlyJS support is now an implemented as extension and patch, wireframe, and normals are retired in favour of PlotlyJS standard names mesh3d, scatter3d, and cone.

## Version 0.10

- New method for union that returns the input meshes as submeshes of the union mesh.
- `parent` back-edge for submeshes.
- embedding implemented for surface meshes.
- Minimum supported Julia version is 1.10

## Version 0.9.1

- `meshcuboid` returns outward oriented meshes for both the `gmsh` and `compsiencemeshes` generators.

## Version 0.9.3

- The `subcharts` API is introduced, which returns an iterable container of tuples comprising the chart of a subentity such as a face or an edge, together with the overlap map from the domain of the subchart to the domain of the original chart. The APi is introduced to enable efficient implementations of local interpolation routines for finite element spaces.
