abstract type AbstractAlbedoModel end

"""
    Gardner2009 <: AbstractAlbedoModel

Albedo parameterization from Gardner & Sharp (2010) that accounts for grain size evolution, 
impurities, and solar zenith angle effects on snow and ice albedo.

Reference:
Gardner, A. S., & Sharp, M. J. (2010). A review of snow and ice albedo and the development 
of a new physically based broadband albedo parameterization. Journal of Geophysical Research: 
Earth Surface, 115(F1).
"""
struct Gardner2009 <: AbstractAlbedoModel 
    "Surface effective grain radius [mm]"
    re::Float64
    "Concentration of light absorbing carbon [ppm1]"
    clabSnow::Float64
    "Solar zenith angle of incident radiation [deg]"
    SZA::Float64
    "Cloud optical thickness"
    COT::Float64
    "Concentration of light absorbing carbon of first ice layer [ppm1]"
    clabIce::Float64
end

"""
    Brun1992 <: AbstractAlbedoModel

Temperature-dependent snow albedo parameterization from Brun et al. (1992) used in the 
CROCUS snow model. Accounts for snow aging and metamorphism effects on albedo.

Reference:
Brun, E., David, P., Sudul, M., & Brunot, G. (1992). A numerical model to simulate 
snow-cover stratigraphy for operational avalanche forecasting. Journal of Glaciology, 38(128), 13-22.
"""
struct Brun1992 <: AbstractAlbedoModel 
    "Surface effective grain radius [mm]"
    re::Float64
end

"""
    Greuell1994 <: AbstractAlbedoModel

Snow and ice albedo parameterization from Greuell & Konzelmann (1994) that includes effects 
of snow depth, temperature, and cloudiness on surface albedo.

Reference:
Greuell, W., & Konzelmann, T. (1994). Numerical modelling of the energy balance and the 
englacial temperature of the Greenland Ice Sheet. Calculations for the ETH-Camp location 
(West Greenland, 1155 m asl). Global and Planetary change, 9(1-2), 91-114.
"""
struct Greuell1994 <: AbstractAlbedoModel 
    "Snow surface density [kg m-3]"
    d::Float64
    "Cloud amount"
    n::Float64
    "Albedo of ice"
    aIce::Float64
    "Albedo of fresh snow"
    aSnow::Float64
end

"""
    Bougamont2005 <: AbstractAlbedoModel

Albedo scheme from Bougamont et al. (2005) that parameterizes snow aging and melt effects 
on surface albedo for ice sheet applications.

Reference:
Bougamont, M., Bamber, J. L., & Greuell, W. (2005). A surface mass balance model for the 
Greenland Ice Sheet. Journal of Geophysical Research: Earth Surface, 110(F4).
"""
struct Bougamont2005 <: AbstractAlbedoModel 
    aIce::Float64
    aSnow::Float64
    EC
    t0wet::Float64
    t0dry::Float64
    K::Float64
end

# function apply_albedo(model::AbstractAlbedoModel, dT, ...)
#     error("Not implemented")
#     return some_prescribed_values
# end