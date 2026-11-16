# Version 1.1

- New primitive `CompScienceMeshes.meshbunny(h; order=1)` (not exported): a closed mesh of the Stanford bunny, flat (`order=1`) or with cubic PN triangles (`order=3`). The scan is downloaded on first use.
- `load_gmsh_mesh` merges the elements of all entities per element type, so physical groups that span several geometric entities are now read completely.
- `connectivity` supports meshes of equal dimension: `D[m,k]` is `op(±1)` when cell `m` and face `k` share the same vertices (with the sign of their relative orientation) and zero otherwise.
- The Dunavant and Gauss triangle quadrature tables now record the polynomial degree each rule integrates exactly.
- `sort_sfc` throws an `ArgumentError` for non-finite point coordinates instead of recursing forever.

# Version 1.0

- Support for curvilinear (higher order) simplicial meshes: new `CompScienceMeshes.CurvilinearMesh` type (not exported) with curvilinear charts for lines and triangles.
  - `load_gmsh_mesh` returns a `CurvilinearMesh` for `element=:line` or `element=:triangle` with `order > 1`.
  - `meshgeo` accepts an `order` keyword and now reads its output through `load_gmsh_mesh`.
  - New primitive `CompScienceMeshes.gmshcircle(radius, delta; order=1)` (not exported).
  - `skeleton`, `boundary`, `connectivity`, and barycentric refinement work for curvilinear meshes.
- New exports `meshorder` and `nodes`.
- `scatter3d` (PlotlyBase extension) can plot meshes of points.

# Version 0.13

- Fixed various bugs to better support curves
- Better support for triangles in a 2D universe
- Tolerances in overlap are now relative to entity sizes
- Refining a curve returns a `BarycentricRefinement`.
- Fixed `isinclosure` for zero-dimensional cells.
- Fixed negative volumes for charts in a 2D universe (#63).
- Hexahedral charts.
- `==` is defined for the mesh types.
- Structured `meshrectangle` and `meshcuboid` pick the edge length closest to the requested one that fits the geometry exactly, instead of approximating the dimensions by multiples of the edge length.
- Restored type stability of `SubMesh` and made its `length` and iteration allocation-free (iteration was quadratic in time and memory since 0.10).

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
