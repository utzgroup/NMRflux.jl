# Working with SpectData

`SpectData` is the central data representation in `NMRflux.jl`. It holds an N dimensional numerical array together with one coordinate vector per dimension, and it subtypes `AbstractArray{T,N}`.

```julia
struct SpectData{T,N} <: AbstractArray{T,N}
    dat::AbstractArray{T,N}
    coord::NTuple{N,AbstractVector}
end
```

Because it implements the Julia array interface, ordinary indexing, slicing and broadcasting all work on it directly, and it can be passed to functions that expect an array. The point of the type is bookkeeping: the values and the physical axes travel together, so no processing step has to be told what its own x axis is.

## 1. Loading a JEOL example

```@example JEOLEg
using NMRflux
using NMRflux.Examples
using Plots: plot, savefig

data_jeol = NMRflux.Examples.Data["Spheroid culture medium"]
jdf_file = joinpath(data_jeol["path"], "yp-5-fu-2.5-100.jdf")

params_jeol, data_td_jeol = NMRflux.load(jdf_file, :JEOL)

plot(coords(data_td_jeol, 1), real.(data_td_jeol.dat);
    xlabel = "time / s",
    ylabel = "signal (a.u.)",
    title = "JEOL FID (real part)"
    )

savefig("spectdata_jeol_fid.svg"); nothing # hide
```

![](spectdata_jeol_fid.svg)

```@example JEOLEg
typeof(data_td_jeol)
```

## 2. Reading the coordinates

`coords(S)` returns the full tuple of coordinate vectors and `coords(S, k)` returns the coordinate of dimension `k`. Both are exported, and both read the same `coord` field you can reach directly.

```@example JEOLEg
t_axis = coords(data_td_jeol, 1)
(first(t_axis), last(t_axis), length(t_axis))
```

Reading `S.coord[k]` gives the same vector. This documentation uses whichever form is clearer in context.

## 3. Indexing and slicing

A range returns a new `SpectData` carrying the matching slice of the coordinate:

```@example JEOLEg
data_td_jeol[1:5]
```

Integer indices are dropped from the coordinate tuple, exactly as they are dropped from the array shape. Slicing a two dimensional dataset at a fixed row therefore hands back a one dimensional `SpectData` with the correct remaining axis. The raw values stay available through `dat`:

```@example JEOLEg
data_td_jeol.dat[1:5]
```

Broadcasts work directly on the object and preserve the coordinates. `real.(S)` comes back as a `SpectData`, not a bare array, and `imag.` behaves the same way:

```@example JEOLEg
real.(data_td_jeol)[1:5]
```

## 4. Building a SpectData by hand

The two argument constructor takes the data and a tuple of coordinates, one per dimension:

```@example JEOLEg
n = 8
SpectData(collect(1.0:n), (range(0.0, step=0.1, length=n),))
```

The single argument constructor fills in unit ranges, which is what processors use internally when they are handed a plain array:

```@example JEOLEg
SpectData(collect(1.0:4.0))
```

## 5. Higher dimensions

For an N dimensional dataset, `dat` becomes an N dimensional array and `coord[k]` holds the physical coordinate for dimension `k`. `FourierTransformPlan` replaces the transformed axis with a frequency axis, `ZeroFill` extends it, and every other processor passes the coordinates through untouched.

```@example JEOLEg
S2 = SpectData(reshape(collect(1.0:12.0), 4, 3),
    (range(0.0, step=0.25, length=4), range(-1.0, step=1.0, length=3)))

coords(S2)
```

Slicing it at a fixed row drops the first coordinate and keeps the second:

```@example JEOLEg
row = S2[2, :]
(size(row.dat), coords(row, 1))
```

The processors that consume these objects are covered in [Classical Processing Pipeline](DataProcessing.md).
