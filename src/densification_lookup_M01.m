function M01 = densification_lookup_M01(densification_coeffs_M01)
switch densification_coeffs_M01
    % ------------------------ Antarctic -------------------------
    % ERA5 v4 (Paolo et al., 2023)
    case "Ant_ERA5v4_Paolo23"
        M01 = [[1.5131, 0.1317, 0.1317, 0.2158]; [1.8422, 0.1688, 2.4979, 0.3225]];

    % ERA5 new albedo_method="BruneLeFebre", sw_absorption_method=1
    case "Ant_ERA5_BF_SW1"
        M01 = [2.2191, 0.2301, 2.2917, 0.2710]; 

    % RACMO callibration, default (Gardner et al., 2023)
    case "Ant_RACMO_GS_SW0"
        M01 = [1.6383, 0.1691, 1.9991, 0.2414];

    % ------------------------- Greenland ------------------------
    % ERA5 new albedo_method="GardnerSharp", sw_absorption_method=0, bare ice
    case "Gre_ERA5_GS_SW0"
        M01 = [1.3566, 0.1350, 1.8705, 0.2290];

    % RACMO callibration, default (Gardner et al., 2023)
    case "Gre_RACMO_GS_SW0"
        M01 = [[1.2691, 0.1184, 1.9983, 0.2511]; [1.4318, 0.1055, 2.0453, 0.2137]];

    %  ismember(albedo_method,["GreuellKonzelmann","BougamontBamber"]) && sw_absorption_method>0
    case "Gre_RACMO_GB_SW1"
        M01 = [2.2191, 0.2301, 2.2917, 0.2710];
   
    % Kuipers Munneke and others (2015) [semi-emperical], Greenland
    case "Gre_KuipersMunneke"
        M01 = [1.042, 0.0916, 1.734, 0.2039];
end