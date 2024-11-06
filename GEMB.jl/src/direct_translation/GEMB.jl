"""
    GEMB(P0, Ta0, V0, dateN, dlw0, dsw0, eAir0, pAir0, S, isrestart)

Run the Glacier Energy and Mass Balance (GEMB) model.

GEMB calculates a 1-D surface glacier mass balance with detailed subsurface processes including:

- Melt water percolation and refreeze
- Pore water retention
- Dynamic albedo with long-term memory  
- Subsurface temperature diffusion
- Subsurface penetration of shortwave radiation

# Arguments
- `P0`: Precipitation [kg m⁻²]
- `Ta0`: Screen level air temperature [K]
- `V0`: Wind speed [m s⁻¹]
- `dateN`: Dates for simulation
- `dlw0`: Downward longwave radiation flux [W m⁻²]
- `dsw0`: Downward shortwave radiation flux [W m⁻²]
- `eAir0`: Screen level vapor pressure [Pa]
- `pAir0`: Screen level air pressure [Pa]
- `S`: Named tuple of model settings
- `isrestart`: Whether this is a restart run

# Documentation
For complete documentation, see: https://github.com/alex-s-gardner/GEMB

# References
Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass
Balance (GEMB): a model of firn processes for cryosphere research, Geosci.
Model Dev., 16, 2277–2302, https://doi.org/10.5194/gmd-16-2277-2023, 2023.
"""
function GEMB(P0, Ta0, V0, dateN, dlw0, dsw0, eAir0, pAir0, S, isrestart)
    dt = (dateN[2]-dateN[1]) * (60*60*24)  # input time step in seconds

    # Constants
    dIce = 910.0  # density of ice [kg m⁻³]

    # Generate model grid
    dz, z_centers = gridInitialize(S.zTop, S.dzTop, S.zMax, S.zY)

    # Initialize model variables
    if isrestart
        m = S.Sizeini
        a = S.Aini        # albedo [fraction]
        adiff = S.Adiffini # albedo [fraction] 
        dz = S.Dzini      # layering
        d = S.Dini        # density [kg m⁻³]
        EC = S.ECini      # surface evaporation (-) condensation (+) [kg m⁻²]
        gdn = S.Gdnini    # grain dentricity
        gsp = S.Gspini    # grain sphericity
        re = S.Reini      # grain size [mm]
        T = S.Tini        # snow temperature [K]
        W = S.Wini        # water content [kg m⁻²]
    else
        m = length(dz)
        a = fill(S.aSnow, m)     # albedo equal to fresh snow [fraction]
        adiff = fill(S.aSnow, m) # albedo equal to fresh snow [fraction]
        d = fill(dIce, m)        # density to that of ice [kg m⁻³]
        EC = 0.0                 # surface evaporation (-) condensation (+) [kg m⁻²]
        gdn = zeros(m)           # grain dentricity to old snow
        gsp = zeros(m)           # grain sphericity to old snow
        re = fill(2.5, m)        # grain size to old snow [mm]
        T = fill(S.Tmean, m)     # initial grid cell temperature to annual mean [K]
        W = zeros(m)             # water content to zero [kg m⁻²]
    end

    F = zeros(m)     # refreeze to zero [kg m⁻²]
    M = zeros(m)     # melt water to zero [kg m⁻²]
    Msurf = 0.0      # initialize surface melt for albedo parameterization
    Ra = zeros(m)    # rain amount to zero [kg m⁻²]

    # Fixed lower temperature boundary condition - T is fixed
    T_bottom = T[end]

    # Determine save time steps
    dateV = vcat(dateN, dateN[end] + (dateN[end]-dateN[end-1]))
    outIdx = if S.outputFreq == "monthly"
        month.(dateV[1:end-1]) .!= month.(dateV[2:end])
    elseif S.outputFreq == "daily"
        day.(dateV[1:end-1]) .!= day.(dateV[2:end])
    elseif S.outputFreq == "3hourly"
        hour.(dateV[1:end-1]) .!= hour.(dateV[2:end])
    end

    # Initialize output structure
    # Single level time series
    varName_monolevel = ["time", "Ta", "P", "M", "R", "F", "EC", "netSW", 
                        "netLW", "shf", "lhf", "a1", "netQ", "re1", "d1", "m", "FAC"]

    Z = fill(NaN, 1, sum(outIdx))
    O = Dict{String,Any}(
        "time" => dateN[outIdx],
        "M" => Z[:], "R" => Z[:], "F" => Z[:], "netSW" => Z[:], "netLW" => Z[:],
        "shf" => Z[:], "lhf" => Z[:], "a1" => Z[:], "netQ" => Z[:], "re1" => Z[:],
        "d1" => Z[:], "Ta" => Z[:], "P" => Z[:], "comp1" => Z[:], "comp2" => Z[:],
        "ps" => Z[:], "m" => Z[:]
    )

    # Time averages/totals
    I = findall(outIdx)
    for i in 1:length(I)
        if i == 1
            O["Ta"][i] = mean(Ta0[1:I[i]]) - 273.15  # convert from K to deg C
            O["P"][i] = sum(P0[1:I[i]])
        else
            O["Ta"][i] = mean(Ta0[(I[i-1]+1):I[i]]) - 273.15
            O["P"][i] = sum(P0[(I[i-1]+1):I[i]])
        end
    end

    # Multi level time series
    varName_profile = ["d", "T", "W", "a", "dz", "re", "gdn", "gsp"]
    Z = fill(NaN, length(d) + S.addCells, length(I))
    for var in varName_profile
        O[var] = copy(Z)
    end
    O["ps"] = copy(Z)

    R = 0.0

    # Initialize cumulative output values
    OV_varNames = ["R", "M", "F", "P", "EC", "Ra", "mAdd", "netSW", "netLW", "shf",
                   "lhf", "a1", "re1", "ulw", "d1", "comp1", "comp2", "m", "netQ", "FAC"]

    for var in OV_varNames
        if !haskey(O, var)
            O[var] = zero(Z)
        end
    end

    # Initialize variables that will be accessed throughout function
    R = 0.0
    M = 0.0 
    F = 0.0
    P = 0.0
    EC = 0.0
    Ra = 0.0
    mAdd = 0.0
    netSW = 0.0
    netLW = 0.0
    shf = 0.0
    lhf = 0.0
    a1 = 0.0
    re1 = 0.0
    ulw = 0.0
    d1 = 0.0
    comp1 = 0.0
    comp2 = 0.0
    m = 0.0
    netQ = 0.0
    FAC = 0.0

    OV = Dict{String,Any}()
    for var in OV_varNames
        OV[var] = 0.0
    end
    OV["count"] = 0

    # Start year loop for model spin up
    @showprogress for yIdx in 1:(S.spinUp + 1)
        # Determine initial mass [kg]:
        initMass = sum(dz .* d) + sum(W)

        # Initialize cumulative variables:
        sumR = 0.0
        sumF = 0.0
        sumM = 0.0
        sumEC = 0.0
        sumP = 0.0
        sumMassAdd = 0.0
        sumMsurf = 0.0
        sumRa = 0.0

        # Start loop for data frequency
        for dIdx in 1:length(dateN)
            # Extract daily data:
            dlw = dlw0[dIdx]     # downward longwave radiation flux [W m⁻²]
            dsw = dsw0[dIdx]     # downward shortwave radiation flux [W m⁻²]
            Ta = Ta0[dIdx]       # screen level air temperature [K]
            P = P0[dIdx]         # precipitation [kg m⁻²]
            V = V0[dIdx]         # wind speed [m s⁻¹]
            eAir = eAir0[dIdx]   # screen level vapor pressure [Pa]
            pAir = pAir0[dIdx]   # screen level air pressure [Pa]

            # Albedo calculations contained in switch to minimize passing of variables
            # to albedo function
            swf = zeros(length(dz))
            
            # Calculate albedo based on scheme
            if S.aIdx in [1,2]
                ccsnowValue = isa(S.ccsnowValue, Vector) ? S.ccsnowValue[dIdx] : S.ccsnowValue
                cciceValue = isa(S.cciceValue, Vector) ? S.cciceValue[dIdx] : S.cciceValue
                cotValue = isa(S.cotValue, Vector) ? S.cotValue[dIdx] : S.cotValue
                szaValue = isa(S.szaValue, Vector) ? S.szaValue[dIdx] : S.szaValue
                dswdiffrf = isa(S.dswdiffrf, Vector) ? S.dswdiffrf[dIdx] : S.dswdiffrf

                # Snow grain metamorphism
                re, gdn, gsp = grainGrowth(T, dz, d, W, re, gdn, gsp, dt, S.aIdx)

                # Calculate snow, firn and ice albedo
                a, adiff = albedo(S.aIdx, re, dz, d, nothing, S.aIce, S.aSnow, S.aValue, S.adThresh,
                                a, adiff, T, W, P, EC, Msurf, ccsnowValue, cciceValue, szaValue,
                                cotValue, nothing, nothing, nothing, dt, dIce)

                # Determine distribution of absorbed sw radiation with depth
                swf = shortwave(S.swIdx, S.aIdx, dsw, dswdiffrf, a[1], adiff[1], d, dz, re, dIce)

            elseif S.aIdx == 3
                cldFrac = isa(S.cldFrac, Vector) ? S.cldFrac[dIdx] : S.cldFrac

                # Calculate snow, firn and ice albedo
                a, adiff = albedo(S.aIdx, re, dz, d, cldFrac, S.aIce, S.aSnow, S.aValue, S.adThresh,
                                a, adiff, nothing, nothing, nothing, nothing, nothing, nothing, nothing,
                                nothing, nothing, nothing, nothing, nothing, nothing, dIce)

                # Determine distribution of absorbed sw radiation with depth
                swf = shortwave(S.swIdx, S.aIdx, dsw, nothing, a[1], adiff[1], d, dz, re, dIce)

            elseif S.aIdx == 4
                # Calculate snow, firn and ice albedo
                a, adiff = albedo(S.aIdx, nothing, nothing, d, nothing, S.aIce, S.aSnow, S.aValue,
                                S.adThresh, a, adiff, T, W, P, EC, nothing, nothing, nothing, nothing,
                                nothing, S.t0wet, S.t0dry, S.K, dt, dIce)

                # Determine distribution of absorbed sw radiation with depth
                swf = shortwave(S.swIdx, S.aIdx, dsw, nothing, a[1], adiff[1], d, dz, re, dIce)
            end

            # Calculate net shortwave [W m⁻²]
            netSW = sum(swf)

            # Calculate new temperature-depth profile and turbulent heat fluxes [W m⁻²]
            shf, lhf, T, EC, ulw = thermo(T, re, dz, d, swf, dlw, Ta, V, eAir, pAir, S.tcIdx,
                                         S.eIdx, S.teValue, S.dulwrfValue, S.teThresh, W[1], dt,
                                         S.dzMin, S.Vz, S.Tz, S.ThermoDeltaTScaling, dIce,
                                         S.isdeltaLWup)

            # Change in thickness of top cell due to evaporation/condensation
            # assuming same density as top cell
            dz[1] = dz[1] + EC / d[1]

            # Add snow/rain to top grid cell adjusting cell depth, temperature and density
            T, dz, d, Ra, W, a, adiff, re, gdn, gsp = accumulation(S.aIdx, S.dsnowIdx, S.Tmean,
                                                                   Ta, T, dz, d, P, W, S.dzMin,
                                                                   S.C, V, S.Vmean, a, adiff,
                                                                   S.aSnow, re, gdn, gsp, dIce)

            # Calculate water production, melt, runoff and resulting changes
            comp2 = sum(dz)
            M, Msurf, R, F, T, d, dz, W, mAdd, _, a, adiff, re, gdn, gsp = melt(T, d, dz, W,
                                                                                 Ra, a, adiff,
                                                                                 S.dzMin, S.zMax,
                                                                                 S.zMin, S.zTop,
                                                                                 S.zY, re, gdn,
                                                                                 gsp, dIce)
            comp2 = comp2 - sum(dz)

            # Allow non-melt densification
            comp1 = sum(dz)
            d, dz = densification(S.denIdx, S.aIdx, S.swIdx, S.adThresh, d, T, dz, S.C, dt,
                                re, S.Tmean, dIce)
            comp1 = comp1 - sum(dz)

            # Calculate net longwave [W m⁻²]
            netLW = dlw - ulw

            # Sum component mass changes [kg m⁻²]
            sumMassAdd += mAdd
            sumM += M
            sumMsurf += Msurf
            sumR += R
            sumW = sum(W)
            sumP += P
            sumEC += EC
            sumRa += Ra
            sumF += F

            # Calculate total system mass
            sumMass = sum(dz .* d)
            FAC = sum(dz .* (dIce .- min.(d,dIce)))/1000

            # Check mass conservation
            dMass = sumMass + sumR + sumW - sumP - sumEC - initMass - sumMassAdd
            dMass = round(dMass * 100)/100

            if dMass != 0
                error("total system mass not conserved in MB function")
            end

            # Check bottom grid cell T is unchanged
            if abs(T[end] - T_bottom) > 1e-8
                @debug "T(end)≠T_bottom"
            end

            if yIdx == S.spinUp + 1
                # Initialize cumulative and average variables for output
                d1 = d[1]
                a1 = a[1]
                re1 = re[1]
                netQ = netSW + netLW + shf + lhf

                OV["R"] += R
                OV["M"] += M 
                OV["F"] += F
                OV["P"] += P
                OV["EC"] += EC
                OV["Ra"] += Ra
                OV["mAdd"] += mAdd
                OV["netSW"] += netSW
                OV["netLW"] += netLW
                OV["shf"] += shf
                OV["lhf"] += lhf
                OV["a1"] += a1
                OV["re1"] += re1
                OV["ulw"] += ulw
                OV["d1"] += d1
                OV["comp1"] += comp1
                OV["comp2"] += comp2
                OV["m"] += m
                OV["netQ"] += netQ
                OV["FAC"] += FAC
                OV["count"] += 1

                if outIdx[dIdx]
                    # Store model output
                    r = sum(outIdx[1:dIdx])

                    for v in OV_varNames
                        if v in ["M", "R", "F", "EC", "P", "Ra", "mAdd", "comp1", "comp2"]
                            O[v][r] = OV[v]
                        else
                            O[v][r] = OV[v] / OV["count"]
                        end
                    end

                    # Instantaneous level data
                    o = length(d) - 1
                    @show o size(O["re"])
                    O["re"][(end-o):end,r] .= re
                    O["d"][(end-o):end,r] .= d
                    O["T"][(end-o):end,r] .= T
                    O["W"][(end-o):end,r] .= W
                    O["dz"][(end-o):end,r] .= dz
                    O["gdn"][(end-o):end,r] .= gdn
                    O["gsp"][(end-o):end,r] .= gsp
                    O["ps"][(end-o):end,r] .= sum(dz) - sumMass/910

                    O["m"][r] = o + 1

                    # Reset cumulative values
                    for v in OV_varNames
                        OV[v] = 0.0
                    end
                    OV["count"] = 0
                end
            end
        end

        # Display cycle completed and time to screen
        # println("$(S.runID): cycle: $yIdx of $(S.spinUp + 1), cpu time: $(round(time() - start_time)) sec, " *
        #         "avg melt: $(round(sumM/(dateN[end]-dateN[1])*365.25)) kg/m2/yr")
    end

    # Save model output and settings
    return O, S
end

function _populate_defaults!(S::Dict{Symbol,Any})

    # Model settings
    get!(S, :spinUp, 2)
    get!(S, :aIdx, 1)
    get!(S, :swIdx, 1)
    get!(S, :denIdx, 2)
    get!(S, :dsnowIdx, 1)
    get!(S, :eIdx, 1)
    get!(S, :tcIdx, 1)
    get!(S, :runID, "test_run")
    
    # Grid parameters
    get!(S, :zTop, 10.0)
    get!(S, :dzTop, 0.05)
    get!(S, :dzMin, S[:dzTop] / 2)
    get!(S, :zMax, 250.0)
    get!(S, :zMin, ceil(S[:zMax]/2 /10)*10)
    get!(S, :zY, 1.10)
    
    # Vertical profile parameters
    get!(S, :Vz, 2.0)#[2.0, 10.0]  # Heights for wind speed profile [m]
    get!(S, :Tz, 0.5)#[0.5, 2.0]   # Heights for temperature profile [m]
    
    # Time parameters
    get!(S, :dt, 3600.0)  # Time step [s]
    
    # Physical parameters
    get!(S, :ThermoDeltaTScaling, 1/11)
    get!(S, :Vmean, 10.0)
    get!(S, :teValue, 1.0)
    get!(S, :teThresh, 0.5)
    get!(S, :isdeltaLWup, false)
    get!(S, :dulwrfValue, 0.0)
    get!(S, :isrestart, false)
    get!(S, :outputFreq, "monthly")
    get!(S, :addCells, 100)  # Maximum number of new cells that can be added
    get!(S, :Tmean, 273.15)  # Mean annual temperature [K]
    get!(S, :C, 500.0)      # Annual accumulation rate [mm w.e. yr⁻¹]
    get!(S, :CtoK, 273.15)  # Celsius to Kelvin conversion
    get!(S, :Ttol, 0.01)    # Temperature tolerance
    get!(S, :Ptol, 1e-11)   # Precipitation tolerance
    get!(S, :Dtol, 1e-11)   # Layer thickness tolerance
    get!(S, :adThresh, 0.1) # Albedo threshold
    
    # Albedo parameters
    get!(S, :aSnow, 0.85)
    get!(S, :aIce, 0.48)
    get!(S, :aValue, S[:aSnow])
    get!(S, :dswdiffrf, 0.0)
    get!(S, :szaValue, 0.0)
    get!(S, :cotValue, 0.0)
    get!(S, :ccsnowValue, 0.0)
    get!(S, :cciceValue, 0.0)
    get!(S, :cldFrac, 0.0)
    get!(S, :t0wet, 0.0)
    get!(S, :t0dry, 0.0)
    get!(S, :K, 0.0)
    
    # Fresh snow parameters
    get!(S, :gdnNew, 1.0)   # Fresh snow grain dendricity
    get!(S, :gspNew, 0.5)   # Fresh snow grain sphericity  
    get!(S, :reNew, 0.1)    # Fresh snow grain radius [mm]

    return S
end

# Example usage
function run_example()
    # Set up model parameters similar to MASTER_RUN.m
    S = Dict{Symbol,Any}()

    _populate_defaults!(S)
    
    # Load test data (this would need to be implemented)
    # For example purposes, creating dummy data
    n = 1000 # number of timesteps
    P0 = zeros(n)
    Ta0 = 273.15 .+ zeros(n) 
    V0 = 10.0 .+ zeros(n)
    dateN = collect(1:n) ./ 365.25 # daily timesteps for 1 year
    dlw0 = 300.0 .+ zeros(n)
    dsw0 = 200.0 .+ zeros(n) 
    eAir0 = 610.0 .+ zeros(n)
    pAir0 = 101325.0 .+ zeros(n)
    
    # Run model
    output = GEMB(P0, Ta0, V0, dateN, dlw0, dsw0, eAir0, pAir0, NamedTuple(S), false)
    
    return output
end
