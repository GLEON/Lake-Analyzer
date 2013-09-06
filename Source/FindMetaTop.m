function [ metaTop_depth ] = FindMetaTop(drho_dz,thermoD,depths,slope)
%----Author: Jordan S Read 2009 ----
%updated 12/16/2009 with thermoD pass

numDepths = length(depths);
metaTop_depth = mean([depths(2) depths(1)]);
Tdepth = NaN(1,numDepths-1);
for i = 1:numDepths-1
    Tdepth(i) = mean([depths(i+1) depths(i)]);
end
[sortDepth,sortInd] = sort([Tdepth thermoD+1e-6]);
drho_dz = interp1(Tdepth,drho_dz,sortDepth);

thermo_index = 1;
thermoId = numDepths;
for i = 1:numDepths
    if eq(thermoId,sortInd(i))
        thermo_index = i;
        break
    end
end

for i = thermo_index:-1:1 %moving up from thermocline index
    if drho_dz(i) < slope %top of metalimnion
        metaTop_depth = sortDepth(i);
        break
    end
end
            % if interp can happen
if thermo_index-i > 1 && drho_dz(thermo_index) > slope  
    metaTop_depth = interp1(drho_dz(i:thermo_index),....
            sortDepth(i:thermo_index),slope);
end


if isnan(metaTop_depth)
    metaTop_depth = min(depths);
end
if metaTop_depth > thermoD
    metaTop_depth = thermoD;
end

end

