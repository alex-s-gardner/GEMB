"""
    accumulation(aIdx, dsnowIdx, Tmean, T_air, T, dz, d, P, W, dz_min, C, V, Vmean, a, adiff, a_SNOW, re, gdn, gsp, dIce)

Adds precipitation and deposition to the model grid.

# Arguments
- `aIdx`: Albedo index
- `dsnowIdx`: Snow density index
- `Tmean`: Mean temperature [K]
- `T_air`: Air temperature [K]
- `T`: Temperature [K]
- `dz`: Grid cell depth [m]
- `d`: Density [kg m⁻³]
- `P`: Precipitation [kg m⁻²]
- `W`: Water content [kg]
- `dz_min`: Minimum grid cell depth [m]
- `C`: Accumulation rate [mm w.e. yr⁻¹]
- `V`: Wind speed [m s⁻¹]
- `Vmean`: Mean wind speed [m s⁻¹]
- `a`: Surface albedo
- `adiff`: Surface albedo for diffuse radiation
- `a_SNOW`: Fresh snow albedo
- `re`: Effective grain radius [mm]
- `gdn`: Grain dendricity
- `gsp`: Grain sphericity
- `dIce`: Ice density [kg m⁻³]

# Returns
- `T`: Updated temperature [K]
- `dz`: Updated grid cell depth [m]
- `d`: Updated density [kg m⁻³]
- `Ra`: Rain amount [kg]
- `W`: Updated water content [kg]
- `a`: Updated surface albedo
- `adiff`: Updated diffuse albedo
- `re`: Updated effective grain radius [mm]
- `gdn`: Updated grain dendricity
- `gsp`: Updated grain sphericity

# Documentation
For complete documentation, see: https://github.com/alex-s-gardner/GEMB

# References
If you use GEMB, please cite the following:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass
Balance (GEMB): a model of firn processes for cryosphere research, Geosci.
Model Dev., 16, 2277–2302, https://doi.org/10.5194/gmd-16-2277-2023, 2023.
"""
function accumulation(aIdx, dsnowIdx, Tmean, T_air, T, dz, d, P, W, dz_min, C, V, Vmean, a, adiff, a_SNOW, re, gdn, gsp, dIce)
    # Constants
    Ttol = 1e-10
    Dtol = 1e-11
    Gdntol = 1e-10
    Ptol = 1e-6

    # Specify constants
    CtoK = 273.15    # Kelvin to Celsius conversion
    dSnow = 150.0    # density of snow [kg m⁻³]
    reNew = 0.05     # new snow grain size [mm]
    gdnNew = 1.0     # new snow dendricity
    gspNew = 0.5     # new snow sphericity
    Ra = 0.0         # rainfall [mm w.e. or kg m⁻³]

    # Density of fresh snow [kg m⁻³]
    if dsnowIdx == 0     # Default value defined above
        # Use default dSnow
    elseif dsnowIdx == 1 # Density of Antarctica snow
        dSnow = 350.0
        #dSnow = 360.0 #FirnMICE Lundin et al., 2017
    elseif dsnowIdx == 2 # Density of Greenland snow, Fausto et al., 2018
        dSnow = 315.0

        #From Vionnet et al., 2012 (Crocus)
        gdnNew = min(max(1.29 - 0.17*V, 0.20), 1.0)
        gspNew = min(max(0.08*V + 0.38, 0.5), 0.9)
        reNew = max(1e-1*(gdnNew/0.99+(1.0-1.0*gdnNew/0.99)*(gspNew/0.99*3.0+(1.0-gspNew/0.99)*4.0))/2.0, Gdntol)
    elseif dsnowIdx == 3 #Surface snow accumulation density from Kaspers et al., 2004, Antarctica
        #dSnow = alpha1 + beta1*T + delta1*C + epsilon1*W
        #     7.36x10-2  1.06x10-3  6.69x10-2  4.77x10-3
        dSnow = (7.36e-2 + 1.06e-3*min(Tmean,CtoK-Ttol) + 6.69e-2*C/1000.0 + 4.77e-3*Vmean)*1000.0
    elseif dsnowIdx == 4 # Kuipers Munneke and others (2015), Greenland
        dSnow = 481.0 + 4.834*(Tmean-CtoK)
    end

    mInit = d .* dz

    if P > 0 + Ptol
        # determine initial mass
        
        # if snow
        if T_air <= CtoK + Ttol
            z_snow = P/dSnow               # depth of snow
            dfall = gdnNew
            sfall = gspNew
            refall = reNew
            
            # if snow depth is greater than specified min dz, new cell created
            if z_snow > dz_min + Dtol
                T = vcat(T_air, T)         # new cell T
                dz = vcat(z_snow, dz)      # new cell dz  
                d = vcat(dSnow, d)         # new cell d
                W = vcat(0.0, W)           # new cell W
                a = vcat(a_SNOW, a)        # new cell a
                adiff = vcat(a_SNOW, adiff)# new cell adiff
                re = vcat(refall, re)      # new cell grain size
                gdn = vcat(dfall, gdn)     # new cell grain dendricity
                gsp = vcat(sfall, gsp)     # new cell grain sphericity
                
            # if snow depth is less than specified minimum dz snow
            else
                mass = mInit[1] + P        # grid cell adjust mass
                dz[1] = dz[1] + P/dSnow    # adjust grid cell depth
                d[1] = mass / dz[1]        # adjust grid cell density
                
                # adjust variables as a linearly weighted function of mass
                # adjust temperature (assume P is same temp as air)
                T[1] = (T_air * P + T[1] * mInit[1])/mass
                
                # adjust a, re, gdn & gsp
                if aIdx > 0
                    a[1] = (a_SNOW * P + a[1] * mInit[1])/mass
                end

                gdn[1] = dfall
                gsp[1] = sfall
                re[1] = max(0.1*(gdn[1]/0.99+(1.0-1.0*gdn[1]/0.99)*(gsp[1]/0.99*3.0+(1.0-gsp[1]/0.99)*4.0))/2, Gdntol)
            end
            
        # if rain    
        else
            # rain is added by increasing the mass and temperature of the ice
            # of the top grid cell. Temperatures are set artifically high to
            # account for the latent heat of fusion. This is the same as
            # directly adding liquid water to the the snow pack surface but
            # makes the numerics easier.
            
            LF = 0.3345e6  # latent heat of fusion(J kg⁻¹)
            CI = 2102.0    # specific heat capacity of snow/ice (J kg⁻¹ k⁻¹)
            
            # grid cell adjust mass
            mass = mInit[1] + P
            
            # adjust temperature
            # liquid: must account for latent heat of fusion
            T[1] = (P *(T_air + LF/CI) + T[1] * mInit[1]) / mass
            
            # adjust grid cell density
            d[1] = mass / dz[1]
            
            # if d > the density of ice, d = dIce
            if d[1] > dIce - Dtol
                d[1] = dIce           # adjust d
                dz[1] = mass / d[1]   # dz is adjusted to conserve mass
            end
            
            Ra = P
        end
        
        # check for conservation of mass
        mass = sum(d .* dz)
        mass_diff = mass - sum(mInit) - P
        mass_diff = round(mass_diff * 100)/100
        
        if mass_diff > 0
            error("mass not conserved in accumulation function")
        end
    end

    return T, dz, d, Ra, W, a, adiff, re, gdn, gsp
end
