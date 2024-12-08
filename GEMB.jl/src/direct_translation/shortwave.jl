"""
    shortwave(swIdx, aIdx, dsw, dswdiff, as, asdiff, d, dz, re, dIce)

Distribute absorbed shortwave radiation vertically within snow/ice.

# Arguments
- `swIdx`: shortwave allowed to penetrate surface (0 = No, 1 = Yes)
- `aIdx`: method for calculating albedo (1-4)
- `dsw`: downward shortwave radiative flux [W m⁻²]
- `dswdiff`: downward shortwave diffuse radiative flux [W m⁻²] 
- `as`: surface albedo
- `asdiff`: surface albedo for diffuse radiation
- `d`: grid cell density [kg m⁻³]
- `dz`: grid cell depth [m]
- `re`: grid cell effective grain radius [mm]
- `dIce`: ice density [kg m⁻³]

# Keyword arguments
- `Dtol`: tolerance for determining if shortwave radiation is absorbed by the top grid cell [kg m⁻³]

# Returns
- `swf`: absorbed shortwave radiation [W m⁻²]
"""
function shortwave(swIdx, aIdx, dsw, dswdiff, as, asdiff, d, dz, re, dIce; Dtol=1e-11)

    # Initialize variables
    m = length(d)
    swf = zeros(m)

    if (swIdx == 0) || ((dIce - d[1]) < Dtol)  # all sw radiation absorbed by top grid cell
        # calculate surface shortwave radiation fluxes [W m⁻²]
        if aIdx == 1  # albedo_method = "gardner_2009"
            swf[1] = (1.0 - as) * max(0.0, (dsw - dswdiff)) + (1.0 - asdiff) * dswdiff
        else
            swf[1] = (1 - as) * dsw
        end
    else  # sw radiation absorbed at depth within glacier
        if aIdx == 2  # albedo_method = "brun_1992" function of effective radius (3 spectral bands)
            # convert effective radius [mm] to grain size [m]
            gsz = (re * 2) / 1000

            # Spectral fractions [0.3-0.8um 0.8-1.5um 1.5-2.8um]
            # (Lefebre et al., 2003)
            sF = [0.606, 0.301, 0.093]

            # initialize variables
            B1_cum = ones(m + 1)
            B2_cum = ones(m + 1)

            # spectral albedos:
            # 0.3 - 0.8um
            a1 = min(0.98, 0.95 - 1.58 * sqrt(gsz[1]))
            # 0.8 - 1.5um
            a2 = max(0.0, 0.95 - 15.4 * sqrt(gsz[1]))
            # 1.5 - 2.8um
            a3 = max(0.127, 0.88 + 346.3*gsz[1] - 32.31*sqrt(gsz[1]))

            # separate net shortwave radiative flux into spectral ranges
            swfS = (sF * dsw) .* (1 .- [a1, a2, a3])

            # absorption coefficient for spectral range:
            h = d ./ sqrt.(gsz)
            B1 = 0.0192 * h                 # 0.3 - 0.8um
            B2 = 0.1098 * h                 # 0.8 - 1.5um
            # B3 = +inf                     # 1.5 - 2.8um

            # cumulative extinction factors
            B1_cum[2:end] = cumprod(exp.(-B1 .* dz))
            B2_cum[2:end] = cumprod(exp.(-B2 .* dz))

            # flux across grid cell boundaries
            Qs1 = swfS[1] * B1_cum
            Qs2 = swfS[2] * B2_cum

            # net energy flux to each grid cell
            swf = (Qs1[1:m] - Qs1[2:m+1]) + (Qs2[1:m] - Qs2[2:m+1])

            # add flux absorbed at surface
            swf[1] = swf[1] + swfS[3]
        else  # function of grid cell density
            # fraction of sw radiation absorbed in top grid cell
            # (wavelength > 0.8um)
            SWs = 0.36

            # calculate surface shortwave radiation fluxes [W m⁻²]
            swf_s = SWs * (1 - as) * dsw
            swf_ss = (1-SWs) * (1 - as) * dsw

            # SW allowed to penetrate into snowpack
            Bs = 10.0    # snow SW extinction coefficient [m⁻¹] (Bassford,2006)
            Bi = 1.3     # ice SW extinction coefficient [m⁻¹] (Bassford,2006)

            # calculate extinction coefficient B [m⁻¹] vector
            B = Bs .+ (300 .- d) .* ((Bs - Bi)/(dIce - 300))

            # cumulative extinction factor
            B_cum = [1.0; cumprod(exp.(-B .* dz))]

            # flux across grid cell boundaries
            Qs = swf_ss * B_cum

            # net energy flux to each grid cell
            swf = (Qs[1:m] - Qs[2:m+1])

            # add flux absorbed at surface
            swf[1] = swf[1] + swf_s
        end
    end

    return swf
end
