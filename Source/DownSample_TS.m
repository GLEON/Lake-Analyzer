function [ DS_dates,DS_varArray ] = DownSample_TS( Dates,outRs,VarArray )
%----Author: Jordan S Read, 2009----
%Dates in matDate format: (datenum('7/1/2009')) will be in matDate format
%outRs: seconds. If outRs < step between dates, original values are
%returned.
%if some dates are not uniuque, values will be averaged***
%

if nargin == 2
    VarArray = [];
    DS_varArray = [];
end

matStep = outRs/86400; %matdates based on day unit, conver seconds/day
u_matStep = 1/matStep;
sortedDates = Dates; %같�---같 assumes values passed are pre-sorted 같--같�

tempDates = floor(sortedDates*u_matStep)*matStep;

[DS_dates,uniStInd] = unique(tempDates,'last');
if nargin == 3
    DS_Length = length(DS_dates);
    DS_varArray = NaN(DS_Length,1);
    start = 1;
    for j = 1:DS_Length
        matchInd = start:uniStInd(j);
        if gt(uniStInd(j),length(VarArray)) || isempty(VarArray(matchInd))
            DS_varArray(j) = NaN;
        else
            DS_varArray(j) = mean(VarArray(matchInd));
        end
        start = uniStInd(j)+1;
    end
end

end