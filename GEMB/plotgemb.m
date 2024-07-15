function plotgemb(results,fieldname,varargin);
	%Example call:
	% plotgemb(a,'T','numlevels',150,'zerolevel',0,'figure',1); colorbar

	options=pairoptions(varargin{:});

	zerolevel= getfieldvalue(options,'zerolevel',0);
	numlevels= getfieldvalue(options,'numlevels',-1);
	maxstep=getfieldvalue(options,'maxstep',length(results.time));
	fignum=getfieldvalue(options,'figure',1);

	%colormap
	c = getcolormap(options);

	%number of results: 
	nt=length(results.time); 
	Hadd=0;
	Madd=0;
	MeanSmb=0;
	if isfield(results,'P') & isfield(results,'EC') & isfield(results,'R');
		results.smb=[results.P + results.EC - results.R];
		MeanSmb=mean(results.smb);
	else
		display('WARNING: no components to calcuate SMB in results, plot surface height will not be adjusted accordingly.');
	end

	for i=1:nt, 

		z0=zerolevel;

		%retrieve time, and delta around time: 
		time=results.time(i)/365;
		if i<nt,
			deltat=results.time(i+1)/365-time;
		else
			deltat=time-results.time(i-1)/365;			
		end

		if i>1
			di=max(min(find(results.dz(:,i-1)==0))-1,1);
			if isfield(results,'mAdd') & isfield(results,'W') & isfield(results,'d');
				Hadd=Hadd+deltat*((results.mAdd(i)-results.W(di,i-1))/results.d(di,i-1));
			else
				display('WARNING: no mAdd, W, or d in results, plot surface height will not be adjusted accordingly.');
			end
			Madd=(Madd-MeanSmb)/1000;
		end

		%figure out number of levels: 
		%dz=results.dz(:,i);
		if Hadd~=0
			dz=[results.dz(:,i) Hadd];
		else
			dz=[results.dz(:,i)];
		end
		if numlevels==-1,
			nlevels=results.m(i); 
		else
			nlevels=numlevels;
		end
		dz=flipud(dz(end-results.m(i)+1:end-results.m(i)+nlevels));

		%retrieve values: 
		field=results.(fieldname)(:,i);
		%T=flipud(field(end-nlevels+1:end));
		if Hadd~=0
			T=[Hadd*nan; flipud(field(end-results.m(i)+1:end-results.m(i)+nlevels-1))];
		else
			T=[flipud(field(end-results.m(i)+1:end-results.m(i)+nlevels))];
		end

		%build vertical values: 
		nz=length(dz); 

		xA=(time-deltat/2)*ones(nz,1);
		xB=(time+deltat/2)*ones(nz,1);
		xC=(time+deltat/2)*ones(nz,1);
		xD=(time-deltat/2)*ones(nz,1);

		zA=zeros(nz,1);
		zB=zeros(nz,1);
		zC=zeros(nz,1);
		zD=zeros(nz,1);

		z0=z0+Madd;
		
		for j=1:nz, 
			zA(j)=z0;
			zB(j)=z0;
			zC(j)=z0+dz(j);
			zD(j)=z0+dz(j);
			z0=z0+dz(j);
      end       

		if i==1
			figure(fignum)
		end
		Xpatch{i}=[xA,xB,xC,xD]';
		Ypatch{i}=[zA,zB,zC,zD]';
		Tpatch{i}=[T,T,T,T]';
		patch(Xpatch{i},Ypatch{i},Tpatch{i},'EdgeColor','none');

		if i>=maxstep,
			break;
		end
    end

	 h = colormap(gca,c);
    
