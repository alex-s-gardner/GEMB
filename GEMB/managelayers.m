function [d, T, dz, W, mAdd, dz_add, addE, a, adiff, m, EI, EW, re, gdn, gsp] = ...
		managelayers(T, d, dz, W, a, adiff, m, EI, EW, dzMin, zMax, zMin, re, gdn, gsp, zTop, zY, CI, LF, CtoK)
% managelayers adjusts the depth and number of vertical layers in the model 
% to ensure that the thickness of any single layer does not exceed thresholds 
% set for the minimum and maximum allowable layer thickness.
%
%% Syntax 
% 
% 
%
%% Description
% 
% 
% 
%% Inputs
% 
%  T       K            Grid cell temperature.
%  d       kg m^-3      Grid cell density.
%  dz      m            Grid cell thickness.
%  W       kg m^-2      Water content. 
%  a       fraction     Albedo. 
%  adiff   fraction     Diffuse albedo.
%  m       kg m^-2      Grid cell mass.
%  EI      J m^-2       Initial energy of snow/ice.
%  EW      J m^-2       Initial energy of water.
%  dzMin   m            Minimum allowable grid spacing.
%  zMax    m            Maximum depth of the total column. 
%  zMin    m            Minimum depth of the total column. 
%  re      mm           Grain size
%  gdn     unitless     Grain dendricity
%  gsp     unitless     Grain sphericity  
%  zTop    m            Thickness of the upper portion of the model grid, in which grid spacing is constant. 
%  zY      unitless     Grid cell stretching parameter for the lower portion of the model grid, in which grid length increases linearly with depth.
%  CI      J kg^-1 K^-1 Specific heat capacity of snow/ice. 
%  LF      J kg^-1      Latent heat of fusion.
%  CtoK    K            273.15 conversion from C to K.
%  
%% Outputs
% 
%  d       kg m^-3      Grid cell density.
%  T       K            Grid cell temperature.
%  W       kg m^-2      Water content. 
%  mAdd:   kg m^-2      Mass added to the column.
%  dz_add: m            Thickness added to the column.
%  addE:   J m^-2       Energy added to the column.
%  a       fraction     Albedo. 
%  adiff   fraction     Diffuse albedo.
%  m       kg m^-2      Grid cell mass.
%  EI      J m^-2       Initial energy of snow/ice.
%  EW      J m^-2       Initial energy of water.
%  re      mm           Grain size
%  gdn     unitless     Grain dendricity
%  gsp     unitless     Grain sphericity  
%  
%% Documentation
% 
% For complete documentation, see: https://github.com/alex-s-gardner/GEMB 
% 
%% References 
% If you use GEMB, please cite the following: 
% 
% Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass
% Balance (GEMB): a model of firn processes for cryosphere research, Geosci.
% Model Dev., 16, 2277-2302, https://doi.org/10.5194/gmd-16-2277-2023, 2023

Dtol = 1e-11; % tolerance for numerical comparison. 

n = numel(T);

Zcum = cumsum(dz);

% A logical "mask" that indicates which cells are in the top layers: 
top_layers = Zcum <= (zTop + Dtol); 

% Define dzMin2 array using the top-layers' dzMin value for the entire column:  
dzMin2 = dzMin * ones(n,1); 

% Overwrite the bottom layers as the cumulative product times the stretching factor: 
dzMin2(~top_layers) = cumprod(zY * ones(sum(~top_layers),1))*dzMin; 

% Define dzMax2 array using the top-layers' dzMin value for the entire column:  
dzMax2 = 2 * dzMin * ones(n,1); 

% In the bottom layers, dzMax2 is the larger of (zY * dzMin2) or (2 * dzMin)  
dzMax2(~top_layers) = max(zY * dzMin2(~top_layers), 2 * dzMin);

%%

% Preallocate a logical array that will be true for any cell to be deleted: 
delete_cell = false(n,1); 

% Check to see if any cells are too small and need to be merged
for i=1:n

	if dz(i) < (dzMin2(i) - Dtol)

		% dz has not met minimum thickness requirements, so we will delete it 
		% and merge its contents into another cell: 
		delete_cell(i) = true; 

		% Detemine the target location for the cell contents to go: 
		if i==n
			% If the very bottom cell (i==n) is too small, find the lowermost
			% cell that isn't going to be deleted:  
			i_target = find(~delete_cell,1,'last'); 
		else
			i_target = i + 1;
		end

		% Move the quantities to the target location. Quantities are
		% calculated as linearly weighted functions of mass: 
		m_new           = m(i) + m(i_target);
		T(i_target)     = (    T(i)*m(i) +     T(i_target)*m(i_target)) / m_new;
		a(i_target)     = (    a(i)*m(i) +     a(i_target)*m(i_target)) / m_new;
		adiff(i_target) = (adiff(i)*m(i) + adiff(i_target)*m(i_target)) / m_new;

		% Use grain properties from lower cell:
		re(i_target)  =  re(i);
		gdn(i_target) = gdn(i);
		gsp(i_target) = gsp(i);

		% Merge with underlying grid cell and delete old cell:
		dz(i_target) = dz(i) + dz(i_target);         % combine cell depths
		d(i_target)  = m_new / dz(i_target);         % combine top densities
		W(i_target)  = W(i) + W(i_target);           % combine liquid water
		m(i_target)  = m_new;                        % combine top masses

	end
end

% Delete combined cells:
m(delete_cell)      = []; 
W(delete_cell)      = []; 
dz(delete_cell)     = []; 
d(delete_cell)      = []; 
T(delete_cell)      = []; 
a(delete_cell)      = [];
re(delete_cell)     = []; 
gdn(delete_cell)    = []; 
gsp(delete_cell)    = []; 
adiff(delete_cell)  = []; 
EI(delete_cell)     = []; 
EW(delete_cell)     = [];
dzMax2(delete_cell) = []; 

% Calculate *new* length of cells:
n = numel(T);

%% Split cells
% * An early implementation of this code used a loop which is included in comments at the bottom of this function for posterity. 

% Find the cells that exceed tolerances: 
f = find(dz > dzMax2+Dtol);

% Conserve quantities among the cells that will be split: 
dz(f) = dz(f)/2; 
W(f)  =  W(f)/2; 
m(f)  =  m(f)/2; 
EI(f) = EI(f)/2; 
EW(f) = EW(f)/2; 

% Sort the indices of all the cells including the ones that will be duplicated:  
fs = sort([(1:n)';f]); 

% Recreate the variables with split cells: 
dz    =    dz(fs); 
W     =     W(fs); 
m     =     m(fs); 
T     =     T(fs); 
d     =     d(fs); 
a     =     a(fs); 
adiff = adiff(fs); 
EI    =    EI(fs); 
EW    =    EW(fs); 
re    =    re(fs); 
gdn   =   gdn(fs); 
gsp   =   gsp(fs); 

%% CORRECT FOR TOTAL MODEL DEPTH
% WORKS FINE BUT HAS BEEN DISABLED FOR CONVIENCE OF MODEL OUTPUT
% INTERPRETATION

% Calculate total model depth:
Ztot = sum(dz);

if Ztot < zMin-Dtol

	% Mass and energy to be added:
	mAdd   = m(end) + W(end);
	addE   = T(end) * m(end) * CI + W(end) * (LF+CtoK*CI);
	dz_add = dz(end);

	% Add a grid cell of the same size and temperature to the bottom:
	dz    = [   dz;    dz(end)];
	T     = [    T;     T(end)];
	W     = [    W;     W(end)];
	m     = [    m;     m(end)];
	d     = [    d;     d(end)];
	a     = [    a;     a(end)];
	adiff = [adiff; adiff(end)];
	EI    = [   EI;    EI(end)];
	EW    = [   EW;    EW(end)];
	re    = [   re;    re(end)];
	gdn   = [  gdn;   gdn(end)];
	gsp   = [  gsp;   gsp(end)];

elseif Ztot > zMax+Dtol

	% Mass and energy loss:
	mAdd   = -(m(end) + W(end));
	addE   = -(T(end) * m(end) * CI) - W(end) * (LF+CtoK*CI);
	dz_add = -(dz(end));

	% Remove a grid cell from the bottom:
	dz(end)    = []; 
	T(end)     = []; 
	W(end)     = []; 
	m(end)     = [];
	d(end)     = []; 
	a(end)     = [];  
	re(end)    = [];  
	gdn(end)   = [];
	gsp(end)   = []; 
	adiff(end) = []; 
	EI(end)    = []; 
	EW(end)    = [];

else
	% No mass or energy is added or removed: 
	mAdd   = 0;
	addE   = 0;
	dz_add = 0;

end

end

%% Old Split Cells loop
% % This loop was in Alex Gardner's original code. Chad Greene vectorized
% % it in July 2024: 
%
% for j=n:-1:1
%     if (dz(j) > dzMax2(j)+Dtol)
% 
%         % split in two
%         dz =    [   dz(1:j-1) ;    dz(j)/2 ;    dz(j)/2 ;    dz(j+1:end)];
%         W =     [    W(1:j-1) ;     W(j)/2 ;     W(j)/2 ;     W(j+1:end)];
%         m =     [    m(1:j-1) ;     m(j)/2 ;     m(j)/2 ;     m(j+1:end)];
%         T =     [    T(1:j-1) ;     T(j)   ;     T(j)   ;     T(j+1:end)];
%         d =     [    d(1:j-1) ;     d(j)   ;     d(j)   ;     d(j+1:end)];
%         a =     [    a(1:j-1) ;     a(j)   ;     a(j)   ;     a(j+1:end)];
%         adiff = [adiff(1:j-1) ; adiff(j)   ; adiff(j)   ; adiff(j+1:end)];
%         EI =    [   EI(1:j-1) ;    EI(j)/2 ;    EI(j)/2 ;    EI(j+1:end)];
%         EW =    [   EW(1:j-1) ;    EW(j)/2 ;    EW(j)/2 ;    EW(j+1:end)];
%         re =    [   re(1:j-1) ;    re(j)   ;    re(j)   ;    re(j+1:end)];
%         gdn =   [  gdn(1:j-1) ;   gdn(j)   ;   gdn(j)   ;   gdn(j+1:end)];
%         gsp =   [  gsp(1:j-1) ;   gsp(j)   ;   gsp(j)   ;   gsp(j+1:end)];
%     end
% end
