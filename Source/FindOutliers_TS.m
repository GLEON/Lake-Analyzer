function [ outlierIndices ] = FindOutliers_TS( TS,varMx,varMn,multiplier )
%----Author: Jordan S Read 2009----
%must have at least 4 points
if nargin == 3
    multiplier = 2.6;
end
if nargin < 3
    varMn = -inf;
    varMx = inf;
    multiplier = 2.6;
end
res = 0.00001;
TS_Length = length(TS);
goodIndices = 1:TS_Length;
outlierIndices = NaN(TS_Length,1);
for i = 1:TS_Length
    if TS(i) < varMn || TS(i) > varMx
        outlierIndices(i) = i;
    end
    if isnan(TS(i))
        outlierIndices(i) = i;
    end  
end
tempOutInd = outlierIndices;
tempOutInd(indexNaN(outlierIndices)) = [];
goodIndices(tempOutInd) = [];

stdev = std(TS(goodIndices));
mu    = mean(TS(goodIndices));

lowBounds = mu-multiplier*stdev;
upBounds  = mu+multiplier*stdev;

for j = goodIndices
        
    if TS(j) - upBounds > res || lowBounds - TS(j) > res
        if isnan(outlierIndices(j))
            outlierIndices(j) = j;
        end
    end
end
outlierIndices(indexNaN(outlierIndices)) = [];