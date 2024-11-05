"""
    melt(T, d, dz, W, Ra, a, adiff, dzMin, zMax, zMin, zTop, zY, re, gdn, gsp, dIce)

Computes the quantity of meltwater due to snow temperature in excess of 0 deg C, determines pore water content and adjusts grid spacing.

# Arguments
- `T`: Temperature [K]
- `d`: Density [kg m⁻³]
- `dz`: Grid cell depth [m]
- `W`: Water content [kg]
- `Ra`: Rain amount [kg]
- `a`: Surface albedo
- `adiff`: Surface albedo for diffuse radiation
- `dzMin`: Minimum grid cell depth [m]
- `zMax`: Maximum model depth [m]
- `zMin`: Minimum model depth [m]
- `zTop`: Top of model [m]
- `zY`: Year-old snow depth [m]
- `re`: Effective grain radius [mm]
- `gdn`: Grain dendricity
- `gsp`: Grain sphericity
- `dIce`: Ice density [kg m⁻³]

# Returns
- `sumM`: Total melt [kg]
- `Msurf`: Surface layer melt [kg]
- `Rsum`: Sum of runoff [kg]
- `Fsum`: Sum of refreeze [kg]
- `T`: Updated temperature [K]
- `d`: Updated density [kg m⁻³]
- `dz`: Updated grid cell depth [m]
- `W`: Updated water content [kg]
- `mAdd`: Mass added/removed to/from base [kg]
- `dz_add`: Thickness of layer added/removed [m]
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
function melt(T, d, dz, W, Ra, a, adiff, dzMin, zMax, zMin, zTop, zY, re, gdn, gsp, dIce)
    # Constants
    Ttol = 1e-10
    Dtol = 1e-11
    Wtol = 1e-13

    ER = 0.0
    sumM = 0.0
    sumER = 0.0
    addE = 0.0
    mSum0 = 0.0
    sumE0 = 0.0
    mSum1 = 0.0
    sumE1 = 0.0
    dE = 0.0
    dm = 0.0
    X = 0
    Wi = 0.0

    # Specify constants
    CtoK = 273.15    # Celsius to Kelvin conversion
    CI = 2102.0      # specific heat capacity of snow/ice (J kg⁻¹ K⁻¹)
    LF = 0.3345e6    # latent heat of fusion (J kg⁻¹)
    dPHC = 830.0     # pore hole close off density [kg m⁻³]

    n = length(T)
    M = zeros(n)
    maxF = zeros(n)
    dW = zeros(n)

    # Store initial mass [kg] and energy [J]
    m = dz .* d                     # grid cell mass [kg]
    EI = m .* T * CI               # initial energy of snow/ice
    EW = W .* (LF .+ CtoK * CI)    # initial energy of water

    mSum0 = sum(W) + sum(m)        # total mass [kg]
    sumE0 = sum(EI) + sum(EW)      # total energy [J]

    # Initialize melt and runoff scalars
    R = 0.0        # runoff [kg]
    Rsum = 0.0     # sum runoff [kg]
    Fsum = 0.0     # sum refreeze [kg]
    sumM = 0.0     # total melt [kg]
    mAdd = 0.0     # mass added/removed to/from base of model [kg]
    addE = 0.0     # energy added/removed to/from base of model [J]
    dz_add = 0.0   # thickness of layer added/removed [m]
    Msurf = 0.0    # surface layer melt

    # Output
    surplusE = 0.0

    # Calculate temperature excess above 0°C
    exsT = max.(0, T .- CtoK)      # [K] to [°C]

    # New grid point center temperature, T [K]
    T = min.(T, CtoK)

    # Specify irreducible water content saturation [fraction]
    Swi = 0.07                     # assumed constant after Colbeck, 1974

    # REFREEZE PORE WATER
    # Check if any pore water
    if sum(W) > 0 + Wtol
        # Calculate maximum freeze amount, maxF [kg]
        maxF = max.(0, .-((T .- CtoK) .* m * CI) / LF)

        # Freeze pore water and change snow/ice properties
        dW = min.(maxF, W)                               # freeze mass [kg]
        W = W .- dW                                      # pore water mass [kg]
        m = m .+ dW                                      # new mass [kg]
        d = m ./ dz                                      # density [kg m⁻³]
        T = T .+ (m .> Wtol) .* (dW .* (LF .+ (CtoK .- T)*CI) ./ (m .* CI))  # temperature [K]

        # If pore water froze in ice then adjust d and dz thickness
        d[d .> dIce-Dtol] .= dIce
        dz = m ./ d
    end

    # Squeeze water from snow pack
    Wi = (dIce .- d) .* Swi .* (m ./ d)    # irreducible water content [kg]
    exsW = max.(0, W .- Wi)                # water "squeezed" from snow [kg]

    # MELT, PERCOLATION AND REFREEZE
    F = zeros(n)

    # Add previous refreeze to F and reset dW
    F = F .+ dW
    dW .= 0

    # Run melt algorithm if there is melt water or excess pore water
    if (sum(exsT) > 0.0 + Ttol) || (sum(exsW) > 0.0 + Wtol)
        # Check if thermal energy exceeds energy to melt entire cell
        # If so redistribute temperature to lower cells (temperature surplus)
        # (Maximum T of snow before entire grid cell melts is a constant 
        # LF/CI = 159.1342)
        surpT = max.(0, exsT .- LF/CI)

        if sum(surpT) > 0.0 + Ttol
            # Calculate surplus energy
            surpE = surpT .* CI .* m
            i = 1

            while sum(surpE) > 0.0 + Ttol && i < n+1
                if i < n
                    # Use surplus energy to increase temperature of lower cell
                    T[i+1] = surpE[i] / m[i+1] / CI + T[i+1]

                    exsT[i+1] = max(0, T[i+1] - CtoK) + exsT[i+1]
                    T[i+1] = min(CtoK, T[i+1])

                    surpT[i+1] = max(0, exsT[i+1] - LF/CI)
                    surpE[i+1] = surpT[i+1] * CI * m[i+1]
                else
                    surplusE = surpE[i]
                    @warn "Surplus energy at the base of GEMB column"
                end

                # Adjust current cell properties (159.1342 is max T)
                exsT[i] = LF/CI
                surpE[i] = 0
                i += 1
            end
        end

        # Convert temperature excess to melt [kg]
        Mmax = exsT .* d .* dz * CI / LF
        M = min.(Mmax, m)              # melt
        Msurf = M[1]
        sumM = max(0, sum(M) - Ra)     # total melt [kg] minus liquid rain

        # Calculate maximum refreeze amount, maxF [kg]
        maxF = max.(0, -((T .- CtoK) .* d .* dz * CI) / LF)

        # Initialize refreeze, runoff, flxDn and dW vectors [kg]
        R = zeros(n)
        flxDn = [R; 0]

        # Determine deepest grid cell where melt/pore water is generated
        X_idx = findlast(@. M > 0.0 + Wtol | exsW > 0.0 + Wtol)
        X = isnothing(X_idx) ? 1 : X_idx

        Xi = 1
        n = length(T)

        # Meltwater percolation
        for i in 1:n
            # Calculate total melt water entering cell
            inM = M[i] + flxDn[i]

            depthice = 0.0
            if d[i] >= dPHC-Dtol
                for l in i:n
                    if d[l] >= dPHC-Dtol
                        depthice += dz[l]
                    else
                        break
                    end
                end
            end

            # Break loop if no meltwater and depth > mw_depth
            if abs(inM) < Wtol && i > X
                break

            # If reaches impermeable ice layer all liquid water runs off (R)
            elseif d[i] >= dIce-Dtol || (d[i] >= dPHC-Dtol && depthice > 0.1+Dtol)
                # No water freezes in this cell
                # No water percolates to lower cell
                # Cell ice temperature & density do not change
                m[i] = m[i] - M[i]                           # mass after melt
                Wi = (dIce-d[i]) * Swi * (m[i]/d[i])        # irreducible water
                dW[i] = max(min(inM, Wi - W[i]), -1*W[i])   # change in pore water
                R[i] = max(0.0, inM - dW[i])                # runoff

            # Check if no energy to refreeze meltwater
            elseif abs(maxF[i]) < Dtol
                m[i] = m[i] - M[i]                           # mass after melt
                Wi = (dIce-d[i]) * Swi * (m[i]/d[i])        # irreducible water
                dW[i] = max(min(inM, Wi - W[i]), -1*W[i])   # change in pore water
                flxDn[i+1] = max(0.0, inM - dW[i])          # meltwater out
                R[i] = 0

            # Some or all meltwater refreezes
            else
                # Change in density and temperature
                # Melt water
                m[i] = m[i] - M[i]
                dz_0 = m[i]/d[i]
                dMax = (dIce - d[i])*dz_0                # d max = dIce
                F1 = min(min(inM,dMax), maxF[i])        # maximum refreeze
                m[i] = m[i] + F1                        # mass after refreeze
                d[i] = m[i]/dz_0

                # Pore water
                Wi = (dIce-d[i])* Swi * dz_0                     # irreducible water
                dW[i] = max(min(inM - F1, Wi-W[i]), -1*W[i])    # change in pore water
                F2 = 0.0

                # THIS HAS NOT BEEN CHECKED
                if dW[i] < 0.0-Wtol                      # excess pore water
                    dMax = (dIce - d[i])*dz_0            # maximum refreeze
                    maxF2 = min(dMax, maxF[i]-F1)        # maximum refreeze
                    F2 = min(-1.0*dW[i], maxF2)          # pore water refreeze
                    m[i] = m[i] + F2                     # mass after refreeze
                    d[i] = m[i]/dz_0
                end

                F[i] = F[i] + F1 + F2

                flxDn[i+1] = max(0.0, inM - F1 - dW[i])  # meltwater out
                if m[i] > Wtol
                    T[i] = T[i] + ((F1+F2)*(LF+(CtoK - T[i])*CI)/(m[i]*CI))
                end

                # Check if an ice layer forms
                if abs(d[i] - dIce) < Dtol
                    # Excess water runs off
                    R[i] = flxDn[i+1]

                    # No water percolates to lower cell
                    flxDn[i+1] = 0
                end
            end

            Xi += 1
        end

        # GRID CELL SPACING AND MODEL DEPTH

        if any(W .< 0.0-Wtol)
            error("Negative pore water generated in melt equations.")
        end

        # Delete all cells with zero mass
        # Adjust pore water
        W = W .+ dW

        # Calculate Rsum
        Rsum = sum(R) + flxDn[Xi]

        # Delete all cells with zero mass
        D = (m .<= 0+Wtol)
        m = m[.!D]
        W = W[.!D]
        d = d[.!D]
        T = T[.!D]
        a = a[.!D]
        re = re[.!D]
        gdn = gdn[.!D]
        gsp = gsp[.!D]
        adiff = adiff[.!D]
        EI = EI[.!D]
        EW = EW[.!D]

        # Calculate new grid lengths
        dz = m ./ d
    end

    Fsum = sum(F)

    # Manage the layering to match user defined requirements
    d, T, dz, W, mAdd, dz_add, addE, a, adiff, m, _, _, re, gdn, gsp = 
        managelayers(T, d, dz, W, a, adiff, m, EI, EW, dzMin, zMax, zMin, re, gdn, gsp, zTop, zY, CI, LF, CtoK)

    # CHECK FOR MASS AND ENERGY CONSERVATION

    # Calculate final mass [kg] and energy [J]
    sumER = Rsum * (LF + CtoK * CI)
    EI = m .* T * CI
    EW = W .* (LF .+ CtoK * CI)

    mSum1 = sum(W) + sum(m) + Rsum
    sumE1 = sum(EI) + sum(EW)

    dm = round((mSum0 - mSum1 + mAdd)*100)/100
    dE = round(sumE0 - sumE1 - sumER + addE - surplusE)

    if dm != 0 || dE != 0
        error("Mass and energy are not conserved in melt equations: dm: $dm dE: $dE")
    end

    if any(W .< 0.0-Wtol)
        error("Negative pore water generated in melt equations.")
    end

    return sumM, Msurf, Rsum, Fsum, T, d, dz, W, mAdd, dz_add, a, adiff, re, gdn, gsp
end
