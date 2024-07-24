%Test Name: SquareShelfSMBGemb
ISSM_compare=false; %Set to true to run ISSM for comparison benchmark

addpath ../
mdt=load('GEMBrestart.mat');

% Use of Gemb method for SMB computation
% See the ISSM SMBgemb class for description of this class and inputs
% or after load call below, type 'smb'
smb = SMBgemb(mdt.mdmesh,mdt.mdgeometry);
smb.dsnowIdx = 1;
smb.swIdx = 1;

% Set to true if run should consider SZA forcing
issza=false;
% Set to true if run should be a restart of a previously saved GEMB output
isrestart=false;

% Other options for testing
%smb.eIdx = 0;
%smb.teValue = md.smb.teValue*0+.95;
%smb.adThresh = 700;
%smb.aIdx=0;

% load hourly surface forcing data from 1965 to 1967:
inputs=load('../../TEST_DATA/gemb_input.mat');

% setup the inputs (duplicated for all the elements of an ISSM shelf model in case these are needed):
smb.Ta=[repmat(inputs.Ta0',mdt.mdmesh.numberofelements,1);inputs.dateN'];
smb.V=[repmat(inputs.V0',mdt.mdmesh.numberofelements,1);inputs.dateN'];
smb.dswrf=[repmat(inputs.dsw0',mdt.mdmesh.numberofelements,1);inputs.dateN'];
smb.dlwrf=[repmat(inputs.dlw0',mdt.mdmesh.numberofelements,1);inputs.dateN'];
smb.P=[repmat(inputs.P0',mdt.mdmesh.numberofelements,1);inputs.dateN'];
smb.eAir=[repmat(inputs.eAir0',mdt.mdmesh.numberofelements,1);inputs.dateN'];
smb.pAir=[repmat(inputs.pAir0',mdt.mdmesh.numberofelements,1);inputs.dateN'];
smb.Vz=repmat(inputs.LP.Vz,mdt.mdmesh.numberofelements,1);
smb.Tz=repmat(inputs.LP.Tz,mdt.mdmesh.numberofelements,1);
smb.Tmean=repmat(inputs.LP.Tmean,mdt.mdmesh.numberofelements,1);
smb.C=repmat(inputs.LP.C,mdt.mdmesh.numberofelements,1);

smb.dswdiffrf=[repmat(0*inputs.Ta0',mdt.mdmesh.numberofelements,1);inputs.dateN'];
smb.cotValue=[repmat(0*inputs.Ta0',mdt.mdmesh.numberofelements,1);inputs.dateN'];
smb.ccsnowValue=[repmat(0*inputs.Ta0',mdt.mdmesh.numberofelements,1);inputs.dateN'];
smb.cciceValue=[repmat(0*inputs.Ta0',mdt.mdmesh.numberofelements,1);inputs.dateN'];
smb.szaValue=[repmat(0*inputs.Ta0',mdt.mdmesh.numberofelements,1);inputs.dateN'];

if issza
	addpath ../../TEST_DATA/SolarAzEl/

	xer=mean(mdt.mdmesh.x(mdt.mdmesh.elements),2);
	yer=mean(mdt.mdmesh.x(mdt.mdmesh.elements),2);
	[later,loner]=xy2ll(xer,yer,+1);
	surfe=mean(mdt.mdgeometry.surface(mdt.mdmesh.elements),2);

	DateVec=datenum([inputs.dateN(1),1,1,0,0,0])+180/1440*[0:(size(smb.Ta,2)-1)];
	for ne=1:(size(smb.szaValue,1)-1)
		[Az El]=SolarAzEl(DateVec,zeros(1,size(DateVec,2))+later(ne),zeros(1,size(DateVec,2))+loner(ne),zeros(1,size(DateVec,2))+surfe(ne)/1000);
		smb.szaValue(ne,:)=max(min(90-El,90),0);
	end
end

%In case of restart, take the last values from a previous run and set as *ini values in the smb class
if isrestart

	tendi=find(cell2mat({mdt.mdresults.time})==max(cell2mat({mdt.mdresults.time})));
	smb.Dzini=mdt.mdresults(tendi).SmbDz; %cell depth (m)
	smb.Dini=mdt.mdresults(tendi).SmbD; %snow density (kg m-3)
	smb.Reini=mdt.mdresults(tendi).SmbRe; %effective grain size (mm)
	smb.Gdnini=mdt.mdresults(tendi).SmbGdn;  %grain dendricity (0-1)
	smb.Gspini=mdt.mdresults(tendi).SmbGsp; %grain sphericity (0-1)
	smb.ECini=mdt.mdresults(tendi).SmbEC; %evaporation/condensation (kg m-2)
	smb.Wini=mdt.mdresults(tendi).SmbW; %Water content (kg m-2)
	smb.Aini=mdt.mdresults(tendi).SmbA; %albedo (0-1)
	smb.Tini=mdt.mdresults(tendi).SmbT; %snow temperature (K) 
	smb.Adiffini=mdt.mdresults(tendi).SmbAdiff; %albedo, diffusive radiation (0-1) 
	tend=mdt.mdresults(tendi).time
	smb.Vmean=mdt.mdsmb.Vmean; %mean annual wind velocity [m s-1]
	smb.C=mdt.mdsmb.C; %mean annual snow accumulation [kg m-2 yr-1]
	smb.Tmean=mdt.mdsmb.Tmean; %mean annual temperature [K]

	Size_ini=zeros(1,length(mdt.mdresults(tendi).SmbDz(:,1)));
	for i=1:length(mdt.mdresults(tendi).SmbDz(:,1))
		Size_ini(i)=0;
		for j=1:length(mdt.mdresults(tendi).SmbDz(1,:))
			if mdt.mdresults(tendi).SmbDz(i,j) > 0
				Size_ini(i)=Size_ini(i)+1;
			end
		end
	end

	Size_ini=transpose(Size_ini);
	smb.Sizeini=Size_ini;
end

addpath ./
S=struct(smb);
S.spinUp = 0;
S.runID = 'GEMBtest';
S.outputFreq='daily';
S.outDIR='./';

%For the GEMB matlab run, we only need the single point value as input
S.zTop=S.zTop(1);
S.dzTop=S.dzTop(1);
S.dzMin=S.dzMin(1);
S.zY=S.zY(1);
S.zMax=S.zMax(1);
S.zMin=S.zMin(1);
S.dswdiffrf=S.dswdiffrf(1,:);
S.szaValue=S.szaValue(1,:);
S.cotValue=S.cotValue(1,:);
S.ccsnowValue=S.ccsnowValue(1,:);
S.cciceValue=S.cciceValue(1,:);
S.Tmean=S.Tmean(1);
S.Vmean=S.Vmean(1);
S.C=S.C(1);
S.Tz=S.Tz(1);
S.Vz=S.Vz(1);
S.aValue=S.aValue(1);
S.teValue=S.teValue(1);
S.dulwrfValue=S.dulwrfValue(1);

if isrestart
	%Initialized values for restart
	S.Dzini=S.Dzini(1,:)';
	S.Dini=S.Dini(1,:)';
	S.Reini=S.Reini(1,:)';
	S.Gdnini=S.Gdnini(1,:)';
	S.Gspini=S.Gspini(1,:)';
	S.ECini=S.ECini(1,:)';
	S.Wini=S.Wini(1,:)';
	S.Aini=S.Aini(1,:)';
	S.Tini=S.Tini(1,:)';
	S.Adiffini=S.Adiffini(1,:)';
	S.Sizeini=S.Sizeini(1);
end

lenrun=8*365;
inputs.dateN=inputs.dateN(1:lenrun);
GEMB(inputs.P0(1:lenrun),inputs.Ta0(1:lenrun),inputs.V0(1:lenrun),[0:1/8:(length(inputs.dateN)/8-1/8)]',inputs.dlw0(1:lenrun),inputs.dsw0(1:lenrun),inputs.eAir0(1:lenrun),inputs.pAir0(1:lenrun),S,isrestart)

a=load('GEMBtest.mat');

if ISSM_compare

	md=triangle(model(),'../../TEST_DATA/Square.exp',350000.);
	md=setmask(md,'all','');
	md=parameterize(md,'../../TEST_DATA/SquareShelf.par');
	md=setflowequation(md,'SSA','all');
	md.materials.rho_ice=910;
	md.cluster=generic('name',oshostname(),'np',3);

	md.smb=smb;

	%smb settings
	md.smb.requested_outputs={'SmbDz','SmbT','SmbD','SmbRe','SmbGdn','SmbGsp','SmbEC','SmbW',...
		'SmbA','SmbAdiff','SmbMassBalance','SmbMAdd','SmbDzAdd','SmbFAC','SmbMeanSHF','SmbMeanLHF',...
		'SmbMeanULW','SmbNetLW','SmbNetSW','SmbWAdd','SmbRunoff','SmbRefreeze','SmbMelt',...
		'SmbEC','SmbPrecipitation','SmbRain','SmbAccumulatedMassBalance','SmbAccumulatedRunoff',...
		'SmbAccumulatedMelt','SmbAccumulatedEC','SmbAccumulatedPrecipitation','SmbAccumulatedRain',...
		'SmbAccumulatedPrecipitation','SmbAccumulatedRefreeze'};

	%only run smb core:
	md.transient.isstressbalance=0;
	md.transient.ismasstransport=0;
	md.transient.isthermal=0;

	%time stepping:
	md.timestepping.start_time=1965;
	md.timestepping.final_time=1966;
	md.timestepping.time_step=1/365.0;
	md.timestepping.interp_forcing=0;

	%Run transient
	md=solve(md,'Transient');
	M=[md.results.TransientSolution(1).SmbAccumulatedMelt diff([md.results.TransientSolution(:).SmbAccumulatedMelt],[],2)];

	%[(1:365)' M(1,:)'*1000*.910 a.M']
	max_meltbias_with_issm=max(abs([M(1,:)'*1000*.910-a.M']))

end

b=load('GEMBtest_output.mat');
max_meltbias_with_gemb_matlab=max(abs([b.M'-a.M']))
