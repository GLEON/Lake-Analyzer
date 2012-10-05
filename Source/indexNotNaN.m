function [ indices ] = indexNotNaN( TS )
%----Author: Jordan S Read, 2009----
%returns indices of array values that are not NaNs

goodIndices = 1:length(TS);

goodIndices(isnan(TS)) = [];

indices = goodIndices;

end

