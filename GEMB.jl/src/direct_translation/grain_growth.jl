function grainGrowth(T, dz, d, W, re, gdn, gsp, dt, aIdx)
# grainGrowth models the effective snow grain size. 
# 
## Syntax 
# 
# 
#
## Description
# 
# 
# 
## Inputs
# 
# * T: grid cell temperature [K]
# * dz: grid cell depth [m]
# * d: grid cell density [kg m-3]
# * W: water content [kg]
# * re: effective grain size [mm]
# * gdn: grain dentricity
# * gsp: grain sphericity
# * dt: time step of input data [s]
# 
## Outputs 
# 
# * re: effective grain size [mm]
# * gdn: grain dentricity
# * gsp: grain sphericity
#
## Documentation
# 
# For complete documentation, see: https://github.com/alex-s-gardner/GEMB 
# 
## References
# Formulations in this function are from the following: 
# 
# DENDRITIC SNOW METAMORPHISM:
# Brun, E., P. David, M. Sudul, and G. Brunot, 1992: A numerical model to
# simulate snow-cover stratigraphy for operational avalanche forecasting.
# Journal of Glaciology, 38, 13-22.
#
# NONDENDRITIC SNOW METAMORPHISM:
# Marbouty, D., 1980: An experimental study of temperature-gradient
# metamorphism. Journal of Glaciology, 26, 303-312.
#
# WET SNOW METAMORPHISM:
# Brun, E., 1989: Investigation on wet-snow metamorphism in respect of
# liquid-water content. Annals of Glaciology, 13, 22-26.
#
# If you use GEMB, please cite the following: 
# 
# Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass 
# Balance (GEMB): a model of firn processes for cryosphere research, Geosci. 
# Model Dev., 16, 2277–2302, https://doi.org/10.5194/gmd-16-2277-2023, 2023.

    Ttol = 1e-10
    Gdntol = 1e-10
    Wtol = 1e-13

    # only when aIdx = 1 or 2 do we run grainGrowth:
    if aIdx != 1 && aIdx != 2
        # come out as we came in:
        return re, gdn, gsp
    end

    ## Function
    gsz = re * 2

    # convert dt from seconds to days
    dt = dt/86400

    # convert T from k to deg C

    # determine liquied-water content in terms
    lwc = W ./ (d .* dz) * 100

    # set maximum water content by mass to 9 percent (Brun, 1980)
    lwc[lwc .> 9+Wtol] .= 9

    ## Calculate temperature gradiant across grid cells
    # returns the average gradinet across the upper and lower grid cell

    # initialize
    dT = zeros(size(T))
    Ti = T

    # depth of grid point center from surface
    zGPC = cumsum(dz) .- dz/2

    # Take forward differences on left and right edges
    m = length(zGPC)
    if m > 2
        dT[1] = (T[3] - T[1])/(zGPC[3]-zGPC[1])
        dT[m] = (T[m] - T[m-2])/(zGPC[m]-zGPC[m-2])
    elseif m > 1
        dT[1] = (T[2] - T[1])/(zGPC[2]-zGPC[1])
        dT[m] = (T[m] - T[m-1])/(zGPC[m]-zGPC[m-1])
    end

    # Take centered differences on interior points
    zGPC = zGPC[3:end] .- zGPC[1:end-2]
    dT[2:end-1] = (T[3:end].-T[1:end-2])./zGPC

    # take absolute value of temperature gradient
    dT = abs.(dT)

    # index for dentricity > 0 & == 0
    G = gdn .> 0+Gdntol
    J = .!G

    ## DENDRITIC SNOW METAMORPHISM
    # FOR SNOW DENTRICITY > 0

    # if there is snow dentricity > 0
    if sum(G) != 0
        # disp ('DENDRITIC DRY SNOW METAMORPHISM')
        # index for dentricity > 0 and T gradients < 5 degC m-1 and >= 5 degC m-1
        H = (abs.(dT) .<= 5+Ttol) .& G .& (W .<= 0+Wtol) # asg not wet accounted for on 19/08/29
        I = (abs.(dT) .> 5+Ttol) .& G .& (W .<= 0+Wtol)

        # determine coefficients
        A = -2E8 * exp.(-6E3 ./ T[H]) * dt
        B = 1E9 * exp.(-6E3 ./ T[H]) * dt
        C = (-2E8 * exp.(-6E3 ./ T[I]) * dt) .* (abs.(dT[I]) .^ 0.4)

        # new dentricity and sphericity for dT < 5 degC m-1
        gdn[H] = gdn[H] .+ A
        gsp[H] = gsp[H] .+ B

        # new dendricity and sphericity for dT >= 5 degC m-1
        gdn[I] = gdn[I] .+ C
        gsp[I] = gsp[I] .+ C

        # WET SNOW METAMORPHISM

        # index for dendritic wet snow
        L = (W .> 0+Wtol) .& G

        # check if snowpack is wet
        if sum(L) != 0
            # disp('DENDRITIC WET SNOW METAMORPHISM')
            # determine coefficient
            D = (1/16) * (lwc[L] .^ 3) * dt

            # new dendricity and sphericity for wet snow
            gdn[L] = gdn[L] .- D
            gsp[L] = gsp[L] .+ D
        end

        # dendricity and sphericity can not be > 1 or < 0
        gdn[gdn .<= 0+Gdntol] .= 0
        gsp[gsp .<= 0+Gdntol] .= 0
        gdn[gdn .>= 1-Gdntol] .= 1
        gsp[gsp .>= 1-Gdntol] .= 1

        # determine new grain size (mm)
        gsz[G] = max.(0.1 .* (gdn[G]./0.99 .+ (1.0.-1.0.*gdn[G]./0.99).*(gsp[G]./0.99*3.0.+(1.0.-gsp[G]./0.99)*4.0)), Gdntol*2)
    end

    # if there is snow dentricity == 0
    if sum(J) != 0
        # disp('NONDENDRITIC SNOW METAMORPHISM')

        # When wet-snow grains (class 6) are submitted to a
        # temperature gradient higher than 5 degC m-1, their sphericity
        # decreases according to Equations (4). When sphericity
        # reaches 0, their size increases according to the functions
        # determined by Marbouty. (Brun et al., 1992)
        P1 = (!).(iszero.(J)) .& (gsp.>0+Gdntol) .& (gsp.<1-Gdntol) .& (abs.(dT) .> 5+Ttol)
        P2 = (!).(iszero.(J)) .& (gsp.>0+Gdntol) .& (gsp.<1-Gdntol) .& ((abs.(dT) .<= 5+Ttol) .& (W .> 0+Wtol))
        P3 = (!).(iszero.(J)) .& (gsp.>0+Gdntol) .& (gsp.<1-Gdntol) .& .!P1 .& .!P2

        F1 = (-2e8 .* exp.(-6e3 ./ T[P1]) .* dt) .* abs.(dT[P1]).^(0.4)
        F2 = (1.0/16.0) * lwc[P2].^(3.0) * dt
        F3 = 1e9 * exp.(-6e3 ./ T[P3]) * dt

        gsp[P1] = gsp[P1] .+ F1
        gsp[P2] = gsp[P2] .+ F2
        gsp[P3] = gsp[P3] .+ F3

        # sphericity can not be > 1 or < 0
        gsp[gsp .<= 0+Gdntol] .= 0
        gsp[gsp .>= 1-Gdntol] .= 1

        # DRY SNOW METAMORPHISM (Marbouty, 1980)
        # grouped model coefficinets from Marbouty, 1980: Figure 9
        P = (!).(iszero.(J)) .& ((W .<= 0+Wtol) .| ((gsp .<=0+Gdntol) .& (abs.(dT) .> 5+Ttol))) # asg not wet accounted for on 19/08/29
        dTi = dT
        Q = Marbouty(Ti[P], d[P], dTi[P])

        # calculate grain growth
        gsz[P] = gsz[P] .+ Q * dt

        # WET SNOW METAMORPHISM (Brun, 1989)

        # index for nondendritic wet snow
        K = (!).(iszero.(J)) .& .!((W .<= 0+Wtol) .| ((gsp .<=0+Gdntol) .& (abs.(dT) .> 5+Ttol)))

        # check if snowpack is wet
        if sum(K) != 0
            # disp('NONDENDRITIC WET SNOW METAMORPHISM')
            # wet rate of change coefficient
            E = (1.28E-8 .+ 4.22E-10 .* (lwc[K].^3)) .* (dt * 86400)   # [mm^3 s^-1]

            # calculate change in grain volume and convert to grain size
            gsz[K] = 2 .* (3/(pi * 4).*((4 / 3)*pi.*(gsz[K]./2).^3 .+ E)).^(1/3)
        end

        # grains with sphericity == 1 can not have grain sizes > 2 mm (Brun, 1992)
        gsz[(abs.(gsp.-1).<Wtol) .& (gsz .> 2-Wtol)] .= 2

        # grains with sphericity == 0 can not have grain sizes > 5 mm (Brun, 1992)
        gsz[(abs.(gsp.-1).>=Wtol) .& (gsz .> 5-Wtol)] .= 5
    end

    # convert grain size back to effective grain radius
    re = gsz/2

    return re, gdn, gsp
end

function Marbouty(T, d, dT)
## calculates grain growth according to Fig. 9 of Marbouty, 1980
# ------NO GRAIN GROWTH FOR d > 400 kg m-3 because (H is set to zero)------
# ---------------this is a major limitation of the model-------------------

## Initialize

    Ttol = 1e-10
    Dtol = 1e-11

    F = zeros(size(T))
    H = F
    G = F

    E = 0.09       # model time growth constant [mm d-1]
    T = T .- 273.15 # converts T from K to ºC
    dT = dT/100.0   # convert dT from degC/m to degC/cm

## Temperature coefficient F
    I = T .> -6+Ttol
    F[I] .= 0.7 .+ ((T[I]/-6) .* 0.3)

    I = (T .<= -6+Ttol) .& (T .> -22+Ttol)
    F[I] .= 1 .- ((T[I].+6)/-16 .* 0.8)

    I = (T .<= -22+Ttol) .& (T .> -40+Ttol)
    F[I] .= 0.2 .- ((T[I].+22)/-18 .* 0.2)

## Density coefficient H

    H[d .< 150-Dtol] .= 1

    I = (d .>= 150-Dtol) .& (d .< 400-Dtol)
    H[I] .= 1 .- ((d[I].-150)/250)

## Temperature gradient coefficient G

    I = (dT .>= 0.16-Ttol) .& (dT .< 0.25-Ttol)
    G[I] .= ((dT[I] .- 0.16)/0.09) .* 0.1

    I = (dT .>= 0.25-Ttol) .& (dT .< 0.40-Ttol)
    G[I] .= 0.10 .+ (((dT[I] .- 0.25)/0.15) .* 0.57)

    I = (dT .>= 0.40-Ttol) .& (dT .< 0.50-Ttol)
    G[I] .= 0.67 .+ (((dT[I] .- 0.40)/0.10) .* 0.23)

    I = (dT .>= 0.50-Ttol) .& (dT .< 0.70-Ttol)
    G[I] .= 0.90 .+ (((dT[I] .- 0.50)/0.20) .* 0.1)

    G[dT .>= 0.7-Ttol] .= 1

## Grouped coefficient Q

    Q = F.*H.*G.*E

    return Q
end
