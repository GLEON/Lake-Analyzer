function [ indices ] = indexNaN( TS )
%----Author: Jordan S Read, 2009----
%returns indices of array values that are NaNs

indices = 1:length(TS);

indices(indexNotNaN(TS)) = [];


end

