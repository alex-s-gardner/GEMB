function Profile_spunup = gemb_spinup(Profile_initial, ClimateForcing, ModelParam, spinup_cycles, display_options)
% gemb_spinup spins up a gemb Profile.
% 
%% Syntax
% 
%  Profile_spunup = gemb_spinup(Profile_initial, ClimateForcing, ModelParam, spinup_cycles)
%  Profile_spunup = gemb_spinup(..., display_waitbar=false)
% 
%% Description
% 
% Profile_spunup = gemb_spinup(Profile_initial, ClimateForcing, ModelParam, spinup_cycles)
% runs the gemb function an integer number of spinup_cycles. Default
% spinup_cycles = 1. 
% 
% Profile_spunup = gemb_spinup(..., display_waitbar=false) disables the
% graphical waitbar. 
% 
%% Example 
% 
%  % Generate sample data: 
%  time_step_hours = 3;
%  ClimateForcing = simulate_climate_forcing("test_1", time_step_hours);
%
%  % Initialize model parameters:
%  ModelParam = model_initialize_parameters(output_frequency="daily");
%    
%  % Initialize grid:
%  Profile_initial = model_initialize_profile(ModelParam, ClimateForcing);
%   
%  % Convert full forcing time series to climatology:
%  ClimateForcingSpinup = forcing_climatology(ClimateForcing); 
%
%  % Spinup a profile for 50 climatological average years: 
%  Profile_spunup = gemb_spinup(Profile_initial, ClimateForcing, ModelParam, 50);
%
%  % Run GEMB: 
%  OutData = gemb(Profile_spunup, ClimateForcing, ModelParam);
%
%% Author Information
% The Glacier Energy and Mass Balance (GEMB) was created by Alex Gardner, with contributions
% from Nicole-Jeanne Schlegel and Chad Greene. Complete code and documentation are available
% at https://github.com/alex-s-gardner/GEMB. Please cite any use of GEMB as:
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): 
% a model of firn processes for cryosphere research, Geosci. Model Dev., 16, 2277–2302, 
% https://doi.org/10.5194/gmd-16-2277-2023, 2023. 
%
% See also gemb_profile and model_initialize_profile. 

%% Input checks

arguments 
    Profile_initial   (:,10)    table {mustContainVariables(Profile_initial, ["z_center","temperature", "dz", "density", "water", "grain_radius", "grain_dendricity", "grain_sphericity", "albedo", "albedo_diffuse"])}
    ClimateForcing    (:,7) timetable {mustContainVariables(ClimateForcing, ["shortwave_downward", "longwave_downward", "temperature_air", "pressure_air", "vapor_pressure", "wind_speed", "precipitation"])}
    ModelParam        (1,1)    struct {mustHaveFields(ModelParam, ["run_prefix", "output_frequency","output_padding","black_carbon_snow","black_carbon_ice","cloud_optical_thickness","solar_zenith_angle","shortwave_downward_diffuse","cloud_fraction","density_ice"])}
    spinup_cycles     (1,1)           {mustBeInteger} = 1
    display_options.display_waitbar (1,1) logical = true
end

%% Spin up 

ModelParam.output_frequency = "last";
Profile = Profile_initial;

if display_options.display_waitbar

    % Create waitbar if running in a graphical environment
    if usejava('desktop')
        h_bar = waitbar(0, 'Spinning up the first GEMB cycle...', 'Name', 'Spinup Progress');
    else
        h_bar = [];
    end
end


for i = 1:spinup_cycles
    % Run GEMB: 
    OutData = gemb(Profile, ClimateForcing, ModelParam, display_waitbar=false);

    % Extract the final profile: 
    Profile = gemb_profile(OutData);

    if display_options.display_waitbar & ~isempty(h_bar)
     
        % Create message: 
        msg = ['Cycle ',num2str(i),' of ',num2str(spinup_cycles),' complete.'];
        
        % Update the bar
        waitbar(i/spinup_cycles, h_bar, msg);
        
    end

end

Profile_spunup = Profile; 

if display_options.display_waitbar
    % Close the progress bar
    if ~isempty(h_bar) && ishandle(h_bar)
        close(h_bar);
    end
end

end


% Custom Validation Function
function mustHaveFields(s, requiredFields)
    for f = requiredFields
        if ~isfield(s, f)
            error('InvalidInput:MissingField', ...
                'Structure is missing the required field: %s', f);
        end
    end
end

% Custom validation function
function mustContainVariables(tab, requiredVars)
    if ~all(ismember(requiredVars, tab.Properties.VariableNames))
        eid = 'Table:MissingVariables';
        msg = 'Input table must contain variables: ' + strjoin(requiredVars, ", ");
        throwAsCaller(MException(eid, msg));
    end
end