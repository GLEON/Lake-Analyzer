function [matchIndices_1,matchIndices_2] = findMatches (TS_1, TS_2,res)
%----Author: Jordan S Read, 2009----
%res is how close the values need to be, in 1Eres
%default matchRes is 10e-6;

matchRes = 10e-6;
if nargin == 3
    matchRes = 1*10^res; %determine how close the values need to be
end


%sort
%assume sorted for now....
%doesn't matter which array has max length...

Length_1 = length(TS_1);
Length_2 = length(TS_2);

matchIndices_1 = NaN(Length_1,1);
matchIndices_2 = NaN(Length_2,1);
%should cut out all values greater than max, less than min

starter = 1;
for j = 1:Length_1
    for k = starter:Length_2
        if abs(TS_1(j)-TS_2(k)) < matchRes
            matchIndices_1(j) = j;
            matchIndices_2(k) = k;
            starter = k;
            %break
        end
    end
end
matchIndices_1 = matchIndices_1(indexNotNaN(matchIndices_1));
matchIndices_2 = matchIndices_2(indexNotNaN(matchIndices_2));

