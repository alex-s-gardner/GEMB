"""
    managelayers(T, d, dz, W, a, adiff, m, EI, EW, dzMin, zMax, zMin, re, gdn, gsp, zTop, zY, CI, LF, CtoK)

Adjusts the depth and number of vertical layers in the model to ensure that the thickness of any single layer 
does not exceed thresholds set for the minimum and maximum allowable layer thickness.

# Arguments
- `T`: Temperature [K]
- `d`: Density [kg m⁻³]
- `dz`: Grid cell depth [m]
- `W`: Water content [kg]
- `a`: Surface albedo
- `adiff`: Surface albedo for diffuse radiation
- `m`: Mass [kg]
- `EI`: Energy of snow/ice [J]
- `EW`: Energy of water [J]
- `dzMin`: Minimum grid cell depth [m]
- `zMax`: Maximum model depth [m]
- `zMin`: Minimum model depth [m]
- `re`: Effective grain radius [mm]
- `gdn`: Grain dendricity
- `gsp`: Grain sphericity
- `zTop`: Top of model [m]
- `zY`: Year-old snow depth [m]
- `CI`: Specific heat capacity [J kg⁻¹ K⁻¹]
- `LF`: Latent heat of fusion [J kg⁻¹]
- `CtoK`: Celsius to Kelvin conversion

# Returns
- `d`: Updated density [kg m⁻³]
- `T`: Updated temperature [K]
- `dz`: Updated grid cell depth [m]
- `W`: Updated water content [kg]
- `mAdd`: Mass added/removed to/from base [kg]
- `dz_add`: Thickness of layer added/removed [m]
- `addE`: Energy added/removed [J]
- `a`: Updated surface albedo
- `adiff`: Updated diffuse albedo
- `m`: Updated mass [kg]
- `EI`: Updated energy of snow/ice [J]
- `EW`: Updated energy of water [J]
- `re`: Updated effective grain radius [mm]
- `gdn`: Updated grain dendricity
- `gsp`: Updated grain sphericity
"""
function managelayers(T, d, dz, W, a, adiff, m, EI, EW, dzMin, zMax, zMin, re, gdn, gsp, zTop, zY, CI, LF, CtoK)
    Dtol = 1e-11

    n = length(T)
    
    Zcum = cumsum(dz)
    
    # Logical mask for top layers
    top_layers = Zcum .<= (zTop + Dtol)
    
    # Initialize dzMin2 array
    dzMin2 = fill(dzMin, n)
    
    # Update bottom layers
    bottom_mask = .!top_layers
    num_bottom = sum(bottom_mask)
    if num_bottom > 0
        dzMin2[bottom_mask] = cumprod(fill(zY, num_bottom)) .* dzMin
    end
    
    # Initialize and update dzMax2 array
    dzMax2 = fill(2 * dzMin, n)
    if num_bottom > 0
        dzMax2[bottom_mask] = max.(zY .* dzMin2[bottom_mask], 2 * dzMin)
    end

    # Track cells to delete
    delete_cell = falses(n)

    # Check for cells that are too small
    for i in 1:n
        if dz[i] < (dzMin2[i] - Dtol)
            delete_cell[i] = true
            
            # Determine target location
            i_target = i == n ? findlast(.!delete_cell) : i + 1
            
            # Merge quantities
            m_new = m[i] + m[i_target]
            T[i_target] = (T[i]*m[i] + T[i_target]*m[i_target]) / m_new
            a[i_target] = (a[i]*m[i] + a[i_target]*m[i_target]) / m_new
            adiff[i_target] = (adiff[i]*m[i] + adiff[i_target]*m[i_target]) / m_new
            
            # Use grain properties from lower cell
            re[i_target] = re[i]
            gdn[i_target] = gdn[i]
            gsp[i_target] = gsp[i]
            
            # Merge remaining properties
            dz[i_target] = dz[i] + dz[i_target]
            d[i_target] = m_new / dz[i_target]
            W[i_target] = W[i] + W[i_target]
            m[i_target] = m_new
        end
    end

    # Remove deleted cells
    # TODO: we can do this without "allocating" new arrays,
    # by defining a function `filter_trunc!(array, filter_idxs)`
    # that first assigns the values in `filter_idxs` to `array[1:sum(filter_idxs)]`,
    # then truncates the array to the new length.
    filter_idx = .!delete_cell
    m = m[filter_idx]
    W = W[filter_idx]
    dz = dz[filter_idx]
    d = d[filter_idx]
    T = T[filter_idx]
    a = a[filter_idx]
    re = re[filter_idx]
    gdn = gdn[filter_idx]
    gsp = gsp[filter_idx]
    adiff = adiff[filter_idx]
    EI = EI[filter_idx]
    EW = EW[filter_idx]
    dzMax2 = dzMax2[filter_idx]
    
    n = length(T)

    # Split cells that are too large
    f = findall(dz .> dzMax2 .+ Dtol)
    
    # Adjust quantities for cells to be split
    dz[f] ./= 2
    W[f] ./= 2
    m[f] ./= 2
    EI[f] ./= 2
    EW[f] ./= 2
    
    # Create new indices including duplicates
    fs = sort(vcat(1:n, f))
    
    # Recreate variables with split cells
    dz = dz[fs]
    W = W[fs]
    m = m[fs]
    T = T[fs]
    d = d[fs]
    a = a[fs]
    adiff = adiff[fs]
    EI = EI[fs]
    EW = EW[fs]
    re = re[fs]
    gdn = gdn[fs]
    gsp = gsp[fs]

    # Correct for total model depth
    Ztot = sum(dz)
    
    if Ztot < zMin - Dtol
        # Add layer at bottom
        mAdd = m[end] + W[end]
        addE = T[end] * m[end] * CI + W[end] * (LF + CtoK*CI)
        dz_add = dz[end]
        
        push!.((dz, T, W, m, d, a, adiff, EI, EW, re, gdn, gsp), 
               getfield.((dz, T, W, m, d, a, adiff, EI, EW, re, gdn, gsp), fill(:, 12) #= TODO: Iterators.repeated(: 12) if this fill does not compile out=#)...)
               
    elseif Ztot > zMax + Dtol
        # Remove bottom layer
        mAdd = -(m[end] + W[end])
        addE = -(T[end] * m[end] * CI) - W[end] * (LF + CtoK*CI)
        dz_add = -dz[end]
        
        pop!.((dz, T, W, m, d, a, re, gdn, gsp, adiff, EI, EW))
        
    else
        mAdd = 0.0
        addE = 0.0
        dz_add = 0.0
    end

    return d, T, dz, W, mAdd, dz_add, addE, a, adiff, m, EI, EW, re, gdn, gsp
end
