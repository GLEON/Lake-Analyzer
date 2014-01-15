function [ metaBot_depth ] = FindMetaBot(drho_dz,thermoD,depths,slope)
%----Author: Jordan S Read 2009 ----
%updated 12/15/2009 with thermoD pass

numDepths = length(depths);
metaBot_depth = depths(numDepths); %default as bottom
Tdepth = NaN(1,numDepths-1);
for i = 1:numDepths-1
    Tdepth(i) = mean([depths(i+1) depths(i)]);
end
[sortDepth,sortInd] = sort([Tdepth thermoD+1e-6]);
drho_dz = interp1(Tdepth,drho_dz,sortDepth);

thermo_index = 1; %check this...
thermoId = numDepths;
for i = 1:numDepths
    if eq(thermoId,sortInd(i))
        thermo_index = i;
        break
    end
end


for i = thermo_index:numDepths %moving down from thermocline index
    if drho_dz(i) < slope %top of metalimnion
        metaBot_depth = sortDepth(i);
        break
    end
end

if i-thermo_index > 1 && drho_dz(thermo_index) > slope
    metaBot_depth = interp1(drho_dz(thermo_index:i),....
        sortDepth(thermo_index:i),slope);
end
if isnan(metaBot_depth)
    metaBot_depth = max(depths);
end
    
end

