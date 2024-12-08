function gridInitialize(z_top::Float64, dz_top::Float64, z_max::Float64, beta::Float64)
    # gridInitialize sets up the initial layer thickness and total grid depth.  
    # 
    ## Syntax 
    # 
    #  dz = gridInitialize(z_top, dz_top, z_max, beta)
    #  [dz,z_center] = gridInitialize(z_top, dz_top, z_max, beta)
    # 
    ## Description
    # 
    # dz = gridInitialize(z_top, dz_top, z_max, beta) creates a 1D grid structure
    # containing the depth dz of each cell in the column, where inputs match 
    # Fig. 1 of Gardner et al., 2023 (https://doi.org/10.5194/gmd-16-2277-2023)
    # and all inputs are scalars as follows: 
    # 
    #  * z_top (m): Thickness of the upper portion of the model grid, in which grid spacing is constant.
    #  * dz_top (m): Spacing of the upper portion of the model grid. 
    #  * z_max (m): Maximum thickness of the total column. 
    #  * beta (unitless): Grid cell stretching parameter for the lower portion of the model grid, in which grid length increases linearly with depth. 
    # 
    # [dz,z_center] = gridInitialize(z_top, dz_top, z_max, beta) also returns a
    # 1D array z_center containing the center depths of each grid cell in
    # meters.
    #
    ## Documentation
    # 
    # For complete documentation, see: https://github.com/alex-s-gardner/GEMB 
    # 
    ## References 
    # If you use GEMB, please cite the following: 
    # 
    # Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass 
    # Balance (GEMB): a model of firn processes for cryosphere research, Geosci. 
    # Model Dev., 16, 2277–2302, https://doi.org/10.5194/gmd-16-2277-2023, 2023.

    ## Error checks:

    # Calculate number of top grid points:
    n_top = z_top/dz_top

    Dtol = 1e-11  # Depth tolerance

    # Check to see if the top grid cell structure length (dz_top) goes evenly 
    # into specified top structure depth (zTop)
    @assert mod(n_top,1) == 0 "Top grid cell structure length does not go evenly into specified top structure depth, adjust dz_top or zTop."

    # Make sure top grid cell structure length (dz_top) is greater than 5 cm
    if dz_top < 0.05-Dtol
        @warn "Initial top grid cell length (dz_top) is < 0.05 m."
    end

    ## Generate grid:

    # Initialize top grid depth vector:
    dzT = fill(dz_top, Int(n_top))

    # Bottom grid cell depth = x*beta^(cells from to structure)

    # Initialize bottom vectors
    dzB = zeros(Int((z_max - z_top)/dz_top))
    gp0 = dz_top
    z0 = z_top
    k = 1

    while z_max > z0+Dtol
        dzB[k] = gp0 * beta
        gp0 = dzB[k]
        z0 = z0 + gp0
        k += 1
    end

    # Delete excess cells from bottom vector:
    filter!(!iszero, dzB)

    # Combine top and bottom dz vectors
    dz = [dzT; dzB]

    # Optional output:
    z_center = -cumsum(dz) .+ dz/2

    return dz, z_center

    ## ---------NEED TO IMPLEMENT A PROPER GRID STRECHING ALGORITHM------------
    # See https://github.com/alex-s-gardner/GEMB/issues/8.
end

