function aveDensity = layerDensity(top,bottom,wtr,depths,bthA,bthD,sal)

%----Author: Jordan S Read, 2009----
% updated 20 Nov 2011, adding timeseries bathymetric effects
%
%Finds the average density thermal layer, bounded by "top" and "bottom"
% where distances are measured from the surface of the water column.
%
%Input:
%   -temps: temperature values (celsius)
%   -depths: corresponding depth values (m)
%   -metaDef: critical slope used to define the metalimnion
%   -top: surface to top of metalimnion
%   -bottom: surface to bottom of metalimnion
%
%Output:
%   -averageEpiDense: average density of the epilimnion (kg/m^3)
%   -thermoDense: density at the thermocline (kg/m^3)

if lt(nargin,7)
    sal = wtr*0;
end
if top > bottom
    error('bottom depth must be greater than top')
end
if ne(length(wtr),length(depths))
    error(['water temperature array must be the same '...
        'length as the depth array'])
elseif nargin < 4
    error('not enough input arguments')
elseif any(isnan([wtr depths bthA bthD]))
    error('input arguments must be numbers')
end

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

dz = 0.1;   % (m)

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
Ao = bthA(Io);
if eq(Ao,0)
    error('Surface area cannot be zero, check *.bth file')
end

%interpolates the bathymetry data
layerD = top:dz:bottom;
layerT  = interp1(depths,wtr,layerD);
layerS  = interp1(depths,sal,layerD);
layerA  = interp1(bthD,bthA,layerD);
layerP  = waterDensity(layerT,layerS);

mass = layerA.*layerP*dz;
aveDensity = sum(mass)/(sum(layerA))/dz;


end
