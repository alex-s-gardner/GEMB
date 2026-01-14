# `densification_lookup_M01` documentation
The `densification_lookup_M01` function returns calibrated coefficients of a densification model. 

# Model Description
See [Section 4.2 Firn air content of Gardner et al (2023)](https://gmd.copernicus.org/articles/16/2277/2023/#section4). The model is based
on Ligtenberg et al (2011)'s Equations 8 and 9, where the output of this
function matches the coefficients MO<sub>550</sub> and MO<sub>830</sub> .

```matlab
 MO_density = MO_density_offset - MO_density_slope * log(accumulation_rate)
```

# Syntax
```matlab
M01 = densification_lookup_M01(densification_coeffs_M01) 
```

# Description

`M01 = densification_lookup_M01(densification_coeffs_M01)` returns
coefficients in the form: 

```matlab
  M01 = [M0_550_offset M0_550_slope M0_830_offset M0_830_slope] 
```
or 

```matlab
  M01 = [M0_550_offset M0_550_slope M0_830_offset M0_830_slope;
         M1_550_offset M1_550_slope M1_830_offset M1_830_slope] 
```

# Example 
Get the coefficients used in Gardner et al., 2023:

```matlab
densification_lookup_M01("Gre_ERA5_GS_SW0")
ans =
      1.3566    0.1350    1.8705    0.2290
      1.4318    0.1055    2.0453    0.2137
```
# References 

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
[https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 

Ligtenberg, S. R. M., Helsen, M. M., and van den Broeke, M. R.: An improved 
semi-empirical model for the densification of Antarctic firn, The Cryosphere, 
5, 809–819, [https://doi.org/10.5194/tc-5-809-2011](https://doi.org/10.5194/tc-5-809-2011), 2011. 

# Author Information
The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
at [https://github.com/alex-s-gardner/GEMB](https://github.com/alex-s-gardner/GEMB). Please cite any use of GEMB as:

Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
[https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023. 