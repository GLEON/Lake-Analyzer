function [varArray,headings] = ...
    fileOpen(fileName,truDelim,headingNum)
%----Author: Jordan S Read, 2009----

    
fauDelim = '*';   %make sure this isn't in the first line of the file.
fID = fopen(fileName);

if eq(nargin,1)
    truDelim = ',';
    headingNum = 1;
end
    

if headingNum >= 0
    firstLine = textscan(fID,'%s',1,'delimiter',fauDelim);
    for i = 1:headingNum-1
        firstLine = textscan(fID,'%s',1,'delimiter',fauDelim);
    end
end
firstData = textscan(fID,'%s',1,'delimiter',fauDelim);

delimCount = 0;
hourDate = 0;
tossHeader = char(firstLine{1,1});
tossData = char(firstData{1,1});
indexer = 1;
for i = 1:length(tossHeader)
    if strcmp(tossHeader(i),truDelim)
        delimCount = delimCount + 1;
    end
end

numBins = delimCount;       %number of delimeters
numBins = numBins + 1;
fclose(fID);
fID = fopen(fileName);

reader = '%s';
for i = 1:numBins
    reader = [reader ' %s'];
end

headings = textscan(fID,reader,headingNum,'delimiter',truDelim);
if headingNum == 0
    headings = [];
end
reader = [];

for i = 1:numBins
    reader = [reader ' %f'];
end


varInfo = textscan(fID,reader,'delimiter',truDelim);

fclose all;

varLength = length(varInfo{1,1}); % length of file read out
varArray = NaN(varLength,numBins);

for j = 1:varLength
    for i = 1:numBins
        varArray(j,i) = varInfo{1,i}(j);
    end
end

end