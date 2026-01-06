%% merge GEMB model output for each point and save to netcdf

%% USER INPUT

varMerge = {'T_air','P','M','R','EC','elev','a1', 'compaction_dens', ...
    'compaction_melt', 'd_50m','sw_net', 'lw_net', 'shf', 'lhf',};
S.run_prefix = 'S2A1D2';
S.inputDIR = '../input/CFSR/T62';

%% Combine output and place into a netcdf

load(fullfile(S.inputDIR, 'mask'))

% find all output data files
f = dir(fullfile('..','Output', [S.run_prefix '*.mat']));
fName = {f(:).name};

% read in time stamp from first file
load(fullfile('..','Output', fName{1}), 'time', 'S')

% read in level data to determine how many vertical levels
load(fullfile('..','Output', fName{1}), 'dz')
nLevel = size(dz,1);
clear dz

for v = 1:length(varMerge)
    var = varMerge{v};
    
    switch var
        case 'P'
            longName    = 'Total precipitaition';
            shortName   = 'precip';
            units       = 'kg/m2';
            monolevel   = true;
            depthAvg    = false;
        case 'elev'
            longName    = 'change in surface elevation';
            shortName   = 'dh';
            units       = 'm';
            monolevel   = true;
            depthAvg    = false;
        case 'R'
            longName    = 'runoff';
            shortName   = 'runoff';
            units       = 'kg/m2';
            monolevel   = true;
            depthAvg    = false;
        case 'M'
            longName    = 'melt';
            shortName   = 'melt';
            units       = 'kg/m2';
            monolevel   = true;
            depthAvg    = false;
        case 'T_air'
            longName    = '2 meter surface temperautre';
            shortName   = 'air';
            units       = 'K';
            monolevel   = true;
            depthAvg    = false;
        case 'a1'
            longName    = 'surface albedo';
            shortName   = 'alb';
            units       = 'unitless';
            monolevel   = true;
            depthAvg    = false;
        case 'compaction_dens'
            longName    = 'elevation lowering due to dry-snow densification';
            shortName   = 'compaction_dens';
            units       = 'm';
            monolevel   = true;
            depthAvg    = false;
        case 'compaction_melt'
            longName    = 'elevation lowering due to wet-snow densification and runoff';
            shortName   = 'compaction_melt';
            units       = 'm';
            mmonolevel   = true;
            depthAvg    = false;
        case 'd_10m'
            longName    = 'depth averaged 10 m density';
            shortName   = 'd_10m';
            units       = 'kg/m3';
            depth       = 10;
            monolevel   = true;
            depthAvg    = true;
        case 'd_50m'
            longName    = 'depth averaged 50 m density';
            shortName   = 'd_50m';
            units       = 'kg/m3';
            depth       = 50;
            monolevel   = true;
            depthAvg    = true;
        case 'd_100m'
            longName    = 'depth averaged 100 m density';
            shortName   = 'd_100m';
            units       = 'kg/m3';
            depth       = 100;
            monolevel   = true;
            depthAvg    = true;
        case 'T'
            longName = 'grid cell temperature';
            shortName = 'temp';
            units = 'K';
            monolevel   = true;
    end
    
    % initialize ouptut array
    if monolevel
        outSize = [size(mask.value), length(time)];
    else
        error('3-d dataset are too large to assemble on a PC')
        outSize = [size(mask.value), nLevel, length(time)];
    end
    OUT = nan(outSize);
    
    % check coordinate size and order
    if sum(size(mask.lat)>1) == 1
        constLat = true;
        if size(mask.lat,1) == 1
            latY = false;
        else
            latY = true;
        end
    else
        constLat= false;
    end
    
    % populate OUT array
    for i = 1:length(fName)
        
        if depthAvg
            % extract depth and density info
            O = load(fullfile('..','Output', fName{i}), 'd', 'S');
            O2 = load(fullfile('..','Output', fName{i}), 'dz', 'S');
            
            % calcualte depth averaged density
            O.d(isnan(O.d)) = 0;
            O2.dz(isnan(O2.dz)) = 0;
            Z = cumsum(O2.dz);
            O.(var) = zeros(size(time));
            for t = 1:length(time)
                idx = Z(:,t) ~= 0;
                O.(var)(t) = (interp1(Z(idx,t),...
                    cumsum(O.d(idx,t) .* O2.dz(idx,t),1), depth))/depth;
            end
            
            % cear uneanted variables
            O = rmfield(O, 'd');
            clear Z O2
   
        else
            O = load(fullfile('..','Output', fName{i}), var, 'S');
        end
        
        if constLat
            if latY
                rIdx = mask.lat == O.S.lat;
                cIdx = mask.lon == O.S.lon;
            else
                cIdx = mask.lat == O.S.lat;
                rIdx = mask.lon == O.S.lon;
            end
            
        else
            [rIdx, cIdx] = find(mask.lat == O.S.lat & mask.lon == O.S.lon);
        end
        
        if monolevel
            OUT(rIdx,cIdx,:) = O.(var);
        else
            OUT(rIdx,cIdx,:,:) = O.(var);
        end
    end
    
    switch var
        case {'elev', 'compaction_dens', 'compaction_melt'}
            OUT = OUT - repmat(OUT(:,:,1), [1,1,size(OUT,3)]);
    end
    
    %% place into netcdf file
    outFileName = fullfile('..','Output', [S.run_prefix '_' var '.nc']);
    [numrow, numcol]= size(mask.value);
    ncid = netcdf.create(outFileName,'NC_SHARE');
    
    if size(mask.lat,1) == size(OUT,1)
        latY = true;
    else
        latY = false;
    end
    
    % add netCDF dimensions
    if latY
        X = mask.lon;
        Y = mask.lat;
        xName = 'lon';
        yName = 'lat';
        xUnits = 'degrees_east';
        yUnits = 'degrees_north';
        yLongName = 'latitude';
        xLongName = 'longitude';
        
    else
        X = mask.lat;
        Y = mask.lon;
        xName = 'lat';
        yName = 'lon';
        yUnits = 'degrees_east';
        xUnits = 'degrees_north';
        xLongName = 'latitude';
        yLongName = 'longitude';
    end
    
    ncY     = netcdf.defDim(ncid, yName, length(Y));
    ncX     = netcdf.defDim(ncid, xName, length(X));
    
    % make time dimension unlimited
    ncT      = netcdf.defDim(ncid, 'time', ...
        netcdf.getConstant('NC_UNLIMITED'));
    
    if monolevel
        dims = [ncY ncX ncT];
    else
        ncL  = netcdf.defDim(ncid, xName, nLevel);
        varidL = netcdf.defVar(ncid, 'level', 'NC_INT', ncL);
        netcdf.putAtt(ncid, varidY, 'long_name', 'grid cell level [not constant depth]')
        netcdf.putAtt(ncid, varidY, 'units', 'NO UNITS')
        
        dims = [ncY ncX ncL ncT];
    end
    
    %  define variables in the new file
    varid  = netcdf.defVar(ncid, shortName, 'NC_DOUBLE', dims);
    varidY = netcdf.defVar(ncid,     yName, 'NC_DOUBLE',  ncY);
    varidX = netcdf.defVar(ncid,     xName, 'NC_DOUBLE',  ncX);
    varidT = netcdf.defVar(ncid,    'time',    'NC_INT',  ncT);
    
    % add attributes
    netcdf.putAtt(ncid, varid, 'long_name', longName)
    netcdf.putAtt(ncid, varid, 'units', units)
    netcdf.putAtt(ncid, varid, 'model_run', ['GEMB0.2_' S.run_prefix])
    netcdf.putAtt(ncid, varid, 'units', units)
    
    netcdf.putAtt(ncid, varidX, 'long_name', xLongName)
    netcdf.putAtt(ncid, varidX, 'units', xUnits)
    
    netcdf.putAtt(ncid, varidY, 'long_name', yLongName)
    netcdf.putAtt(ncid, varidY, 'units', yUnits)
    
    netcdf.putAtt(ncid, varidT, 'long_name', 'date')
    netcdf.putAtt(ncid, varidT, 'units', 'days since 0000-0-0 00:00'), 
    
    
    % leave define mode and enter data mode to write data.
    netcdf.endDef(ncid);
    
    if monolevel
        startIdx = [0 0 0];
        countIdx = [length(Y) length(X) length(time)];
    else
        startIdx = [0 0 0 0];
        countIdx = [length(Y) length(X) nLevel length(time)];
        netcdf.putVar(ncid, varidL, 0, nLevel, 1:nLevel);
    end
    
    netcdf.putVar(ncid, varidY, 0, length(Y), Y);
    netcdf.putVar(ncid, varidX, 0, length(X), X);
    
    netcdf.putVar(ncid, varidT, 0, length(time), time);
    netcdf.putVar(ncid, varid, startIdx, countIdx, OUT);
    netcdf.close(ncid);
end
