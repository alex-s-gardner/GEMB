GEMB output, value summed per output timestep
For values with layers, above the surface is padded with nan
e.g., to get the full layer of densities for time step 1: density_column=d(end-m+1:end,1);
and the surface density at time step 1 (or d1(1) ) should be equal to density_column(1)

		EC: surface evaporation (-) condensation (+) [kg m-2]
		F: refreeze of column [kg m-2]
		FAC: firn air content within column [m]
		M: melt of column [kg m-2]
		P: precipitation input [kg m-2]
		R: runoff of column [kg m-2]
		Ra: rain added to column [kg m-2]
		S: structure of all model input
		T: temperature of each layer [K]
		Ta: 2m air temperature input [K]
		W: liquid water content of each layer [kg m-2]
		a1: albedo of the top layer [fraction]
		comp1: dry snow ccompaction [m]
		comp2: melt compaction [m]
		d: density of each layer [kg m-3]
		d1: density of the top layer [kg m-3]
		dz: layer depths [m]
		gdn: grain dendricity within each layer
		gsp: grain sphericity within each layer
		lhf: latent heat flux [W m-2]
		m: number of snow layers
		mAdd: mass added/removed to/from base of model [kg m-2]
		netLW: net longwave radiation [W m-2]
		netQ: net energy (netSW + netLW + shf + lhf) [W m-2]
		netSW: net shortwave radiation [W m-2]
		ps: pore space within each layer [m^3/m^2]
		re: optically equivelant grain radius within each layer [mm]
		re1: optically equivelant grain radius of the top layer [mm]
		shf: sensible heat flux [W m-2]
		time: time since beginning of simulation [day]
		ulw: upward longwave radiation [W m-2]
