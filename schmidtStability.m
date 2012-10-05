function St = schmidtStability(wtr,depths,bthA,bthD,sal)
%----Author: Jordan S Read 2009 ----
% updated 7 March 2011, adding salinity
% updated 20 Nov 2011, adding timeseries bathymetric effects

%equation provided by:
%°° Idso, S.B., 1973. On the concept of lake stability. °°
%°° Limnol. Oceanogr. 18: 681–683.                      °°

%   St = (g/A0)* [Integral from 0 to hm: (hv-h)*A(h)*p(h)*dh]

% The Lake Analyzer program reverses the direction of this calculation so
% that z = 0 is the surface, and positive z is downwards

if ne(length(wtr),length(depths))
    error(['water temperature array must be the same '...
        'length as the depth array'])
elseif nargin < 4
    error('not enough input arguments')
elseif any(isnan([wtr depths bthA bthD]))
    error('input arguments must be numbers')
end
if eq(nargin,4)
    sal = wtr*0;
end

g = 9.81;   % (m s^-2)
dz = 0.1;   % (m)
% if bathymetry has negative values, intepolate to 0
if lt(min(bthD),0)
    useI = ge(bthD,0);
    if ~eq(bthD,0)
        depT = [0 bthD(useI)];
    else
        depT = bthD(useI);
    end
    bthA = interp1(bthD,bthA,depT);
    bthD = depT;
end

numD = length(wtr);
if max(bthD) > depths(numD)
    wtr(numD+1) = wtr(numD);
    sal(numD+1) = sal(numD);
    depths(numD+1) = max(bthD);
elseif max(bthD) < depths(numD)
    bthD = [bthD depths(numD)];
    bthA = [bthA 0];
end

    
if min(bthD)<depths(1)
    wtr = [wtr(1) wtr];
    sal = [sal(1) sal];
    depths = [min(bthD) depths];
end

[Zo,Io] = min(depths);
Zm = max(depths);
Ao = bthA(Io);
if eq(Ao,0)
    error('Surface area cannot be zero, check *.bth file')
end
rhoL = waterDensity(wtr,sal);
%interpolates the bathymetry data
layerD = Zo:dz:Zm;
layerP  = interp1(depths,rhoL,layerD);
layerA  = interp1(bthD,bthA,layerD);

%find depth to the center of volume
Zv = layerD.*layerA*dz;                    
Zcv = sum(Zv)/sum(layerA)/dz;

numInt = length(layerA);
st = NaN(numInt,1);
for i = 1:numInt
    z = layerD(i);
    A = layerA(i);
    st(i) = -(Zcv-z)*layerP(i)*A*dz;
end
St = g/Ao*sum(st);


end