## Charts

Charts represent maps from parameter domains to configuration space. Their main use is to
provide the geometric description of the cells that make up meshes.

The Chart concept is defined by the following API

```@docs
dimension(ch::CompScienceMeshes.Simplex)
universedimension
volume(::CompScienceMeshes.Simplex)
neighborhood
```

The most important example of a Chart is a Simplex. For Simplices, the following additional
functions are available

```@docs
simplex
center
```
