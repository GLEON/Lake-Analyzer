function aveTemperature = layerTemperature(top,bottom,wtr,depths,bthA,bthD)

%----Author: Jordan S Read, 2013----
% updated 20 Nov 2013, adding timeseries bathymetric effects
%
%Finds the average temperature of a thermal layer, bounded by "top" and "bottom"
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
%   -aveTemperature: average temperature the layer bounded by top & bottom

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
    depths(numD+1) = max(bthD);
elseif max(bthD) < depths(numD)
    bthD = [bthD depths(numD)];
    bthA = [bthA 0];
end
if min(bthD)<depths(1)
    wtr = [wtr(1) wtr];
    depths = [min(bthD) depths];
end

[Zo,Io] = min(depths);
Ao = bthA(Io);
if eq(Ao,0)
    error('Surface area cannot be zero, check *.bth file')
end

%interpolates the bathymetry data
layerD = top:dz:bottom;
layerA  = interp1(bthD,bthA,layerD);
layerT  = interp1(depths,wtr,layerD);


weightedT = layerA.*layerT*dz;
aveTemperature = sum(weightedT)/(sum(layerA))/dz;


end