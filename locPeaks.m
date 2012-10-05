function [peaks,locs] = locPeaks(dataIn,dataMn)

%----Author: Jordan S Read 2011 ----

% this program attempts to mirror 'findpeaks.m' from the signal processing
% toolbox

% dataIn: vector of input data
% dataMn: threshold for peak height

% finds multiple peaks for dataIn
% peaks: peak values
% locs:  indices of peaks

% -- description --
% a peak is a peak if it represents a local maximum

varL = length(dataIn);
locs = false(1,varL);
peaks= NaN(1,varL);

for i = 2:varL-1
    [posPeak,pkI] = max(dataIn(i-1:i+1));
    if eq(pkI,2)
        peaks(i) = posPeak;
        locs(i)  = true;
    end
end

inds = 1:varL;
locs = inds(locs);
peaks= peaks(locs);

% remove all below threshold value

useI = ge(peaks,dataMn);
peaks= peaks(useI);
locs = locs(useI);