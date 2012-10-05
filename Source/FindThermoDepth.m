function [ thermoD,thermoInd,drho_dz,SthermoD,SthermoInd ] = ...
    FindThermoDepth( rhoVar,depths,Smin )
%----Author: Jordan S Read 2009 ----
% updated 3 march 2011

% removed signal processing toolbox function 'findpeaks.m' replaced with
% 'locPeaks.m' which was written to provide the same functionality

seasonal = false;
if nargout > 3
    seasonal = 1;
end
if nargin < 3
    Smin = 0.1;
end
dRhoPerc = 0.15; %min percentage max for unique thermocline step
numDepths = length(depths);
drho_dz = NaN(1,numDepths-1);

for i = 1:numDepths-1
    drho_dz(i) = (rhoVar(i+1)-rhoVar(i))/...
        (depths(i+1)-depths(i));
end
if seasonal
    %look for two distinct maximum slopes, lower one assumed to be seasonal
    [mDrhoZ,thermoInd] = max(drho_dz);          %find max slope
        thermoD = mean([depths(thermoInd)...
        depths(thermoInd+1)]);                  %depth of max slope
    if thermoInd > 1 && thermoInd < numDepths-1 %if within range, 
        Sdn = -(depths(thermoInd+1)-depths(thermoInd))/...
            (drho_dz(thermoInd+1)-drho_dz(thermoInd));
        Sup = (depths(thermoInd)-depths(thermoInd-1))/...
            (drho_dz(thermoInd)-drho_dz(thermoInd-1));
        upD  = depths(thermoInd);
        dnD  = depths(thermoInd+1);
        if ~any([isinf(Sup) isinf(Sdn)])
            thermoD = dnD*(Sdn/(Sdn+Sup))+upD*(Sup/(Sdn+Sup));
        end
    end
    dRhoCut = max([dRhoPerc*mDrhoZ Smin]);
    [pks,locs] = locPeaks(drho_dz,dRhoCut);
    if isempty(pks)
        SthermoD = thermoD;
        SthermoInd = thermoInd;
    else
        mDrhoZ = pks(length(pks));
        SthermoInd = locs(length(pks));
        if SthermoInd > thermoInd+1
            SthermoD = mean([depths(SthermoInd)...
                depths(SthermoInd+1)]);
            if SthermoInd > 1 && SthermoInd < numDepths-1
                Sdn = -(depths(SthermoInd+1)-depths(SthermoInd))/...
                    (drho_dz(SthermoInd+1)-drho_dz(SthermoInd));
                Sup = (depths(SthermoInd)-depths(SthermoInd-1))/...
                    (drho_dz(SthermoInd)-drho_dz(SthermoInd-1));
                upD  = depths(SthermoInd);
                dnD  = depths(SthermoInd+1);
                if ~any([isinf(Sup) isinf(Sdn)])
                    SthermoD = dnD*(Sdn/(Sdn+Sup))+upD*(Sup/(Sdn+Sup));
                end
            end
        else
            SthermoD = thermoD;
            SthermoInd = thermoInd;
        end
    end
    if SthermoD < thermoD;
        SthermoD = thermoD;
        SthermoInd = thermoInd;
    end
    
else
    [mDrhoZ,thermoInd] = max(drho_dz);          %find max slope
        thermoD = mean([depths(thermoInd)...
        depths(thermoInd+1)]);                  %depth of max slope
    if thermoInd > 1 && thermoInd < numDepths-1 %if within range, 
        Sdn = -(depths(thermoInd+1)-depths(thermoInd))/...
            (drho_dz(thermoInd+1)-drho_dz(thermoInd));
        Sup = (depths(thermoInd)-depths(thermoInd-1))/...
            (drho_dz(thermoInd)-drho_dz(thermoInd-1));
        upD  = depths(thermoInd);
        dnD  = depths(thermoInd+1);
        if ~any([isinf(Sup) isinf(Sdn)])
            thermoD = dnD*(Sdn/(Sdn+Sup))+upD*(Sup/(Sdn+Sup));
        end
    end
end

end
