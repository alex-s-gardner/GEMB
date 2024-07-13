function [d, T, dz, W, mAdd, dz_add, addE, a, adiff, m, EI, EW, re, gdn, gsp] = ...
    managelayers(T, d, dz, W, a, adiff, m, EI, EW, dzMin, zMax, zMin, re, gdn, gsp, zTop, zY, CI, LF, CtoK)
%% MANAGE LAYERS

% Description:
% Manages the depth and number of vertical layers in the model

%%%%%%%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%zY = 1.025;
%dzMin = zeros(size(dz))+dzMin;
%z = cumsum(dz);
%n1 = sum(z<10)+1;
%n2 = length(z);
%dzMin(n1:n2) = zY.^(1:(n2-n1+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dtol = 1e-11;
Wtol = 1e-13;

mAdd = 0.0;
addE = 0.0;
dz_add = 0.0;

n=length(T);
zY2=zY;
X1=1;
X2=1;
dzMin2=zeros(size(dz));

Delflag=-99999;

X=1;
Zcum = cumsum(dz); 
% check if depth is too small
dzMin2(1)=dzMin;
for i=2:n
    if (Zcum(i)<=zTop+Dtol)
        dzMin2(i)=dzMin;
        X=i;
    else
        dzMin2(i)=zY2*dzMin2(i-1);
    end
end

% Check to see if any cells are too small and need to be merged
for i=1:n

    if (i<=X && dz(i)<dzMin-Dtol) || (i>X && dz(i)<dzMin2(i)-Dtol)
        if i==n
            X2=i;
            %find closest cell to merge with
            for j=n-1:-1:1
                if m(j)~=Delflag
                    X1=j;
                    break;
                end
            end
        else
            X1=i+1;
            X2=i;
        end

        % adjust variables as a linearly weighted function of mass
        m_new = m(X2) + m(X1);
        T(X1) = (T(X2)*m(X2) + T(X1)*m(X1)) / m_new;
        a(X1) = (a(X2)*m(X2) + a(X1)*m(X1)) / m_new;
        adiff(X1) = (adiff(X2)*m(X2) + adiff(X1)*m(X1)) / m_new;

        %use grain properties from lower cell
        re(X1) = re(X2);
        gdn(X1) = gdn(X2);
        gsp(X1) = gsp(X2);

        %merge with underlying grid cell and delete old cell
        dz(X1) = dz(X2) + dz(X1);                 % combine cell depths
        d(X1) = m_new / dz(X1);                   % combine top densities
        W(X1) = W(X1) + W(X2);                    % combine liquid water
        m(X1) = m_new;                            % combine top masses

        % set cell to -99999 for deletion
        m(X2) = Delflag;
    end
end

% delete combined cells
D = (m <= Delflag+Wtol);
m(D)      = []; 
W(D)      = []; 
dz(D)     = []; 
d(D)      = []; 
T(D)      = []; 
a(D)      = [];
re(D)     = []; 
gdn(D)    = []; 
gsp(D)    = []; 
adiff(D)  = []; 
EI(D)     = []; 
EW(D)     = [];
dzMin2(D) = []; % <- EDIT This line added by Chad Greene, July 2024.

% check if any of the cell depths are too large
n = length(T);
dzMax2 = zeros(size(dz));
X = 1;
Zcum = cumsum(dz); 
dzMax2(1) = dzMin*2.0;
for i=2:n
    if Zcum(i)<=zTop+Dtol
        dzMax2(i)=dzMin*2.0;
        X=i;
    else
        dzMax2(i)=max(zY2*dzMin2(i-1),dzMin*2.0);
    end
end

%% Split cells
% 
% % This loop was in Alex Gardner's original code. Chad Greene vectorized
% % it in July 2024: 
%
% for j=n:-1:1
%     if (j<X && dz(j) > dzMax2(j)+Dtol) || (dz(j) > dzMax2(j)*zY2+Dtol)
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

% The rest of this code section is what Chad Greene wrote to replace the
% loop above: 

% Find the cells that exceed tolerances: 
f = find((1:n)'<X & ( dz > dzMax2+Dtol) | (dz > dzMax2*zY2+Dtol));

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
 
% calculate total model depth
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

end

