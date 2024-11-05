```@index
```

```@docs
ReferenceStore
```

## Correction interface

Kerchunk files often need corrections to the metadata.  

For example, the CF-convention `add_offset` and `scale_factor` metadata fields are stored as separate variables in the source data, but should ideally be stored as a single Zarr `FixedScaleOffset` filter so you can get performance as close to native as possible.  Some CF datasets also encode an `_Unsigned` metadata field, which should simply be used to edit the `dtype` of the Zarr array.

Kerchunk also sometimes places the compressor as the last filter, which is technically compliant with Zarr v3 but is not compliant with Zarr v2.  This is corrected by moving the compressor to the `compressor` field of the metadata, but this has to be done before the Zarr is loaded.

This is the point of the correction interface.  As more idiosyncrasies are discovered, they can be added to it.

```@docs
do_correction!
add_scale_offset_filter_and_set_mask!
move_compressor_from_filters!
```

