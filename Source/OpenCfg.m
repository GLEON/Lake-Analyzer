function [ outPuts,outRs,maxZ,wndH,wndAv,lyrAv,outWn,wtrMx,wtrMn,...
    wndMx,wndMn,drhDz,Tdiff,plotYes,writeYes ] = OpenCfg( LakeName,Year )
%----Author: Jordan S Read 2009 ----

%trouble with commentStyle***
%need to get rid of whitespace in text inputs
fileName = [Year '/' LakeName '.lke'];
fID = fopen(fileName);
if fID ==-1
    fileName = [LakeName '.lke'];
end
fID = fopen(fileName);
if fID ==-1
    error([LakeName '.lke file not found'])
end

TPuts = textscan(fID,'%[^#]','delimiter',',','HeaderLines',2);
tOut = char(TPuts{1});
strL = length(tOut);
inD = [0 strfind((tOut),',') strL];
outPuts = cell(1,length(inD)-1);
chk = true(length(inD)-1,1);
cnt = 1;
for j = 1:length(inD)-1;
    if~isempty(strtrim(tOut(inD(j)+1:inD(j+1)-1)))
        outPuts{cnt} = strtrim(tOut(inD(j)+1:inD(j+1)-1));
        cnt = cnt+1;
    else
        chk(j) = false;
    end
end
fgets(fID); % advance one line
outPuts = outPuts(chk);

outRs = subString(fgets(fID));
maxZ  = subString(fgets(fID));
wndH  = subString(fgets(fID));
wndAv = subString(fgets(fID));
lyrAv = subString(fgets(fID));
outWn = subString(fgets(fID));
wtrMx = subString(fgets(fID));
wtrMn = subString(fgets(fID));
wndMx = subString(fgets(fID));
wndMn = subString(fgets(fID));
drhDz = subString(fgets(fID));
Tdiff = subString(fgets(fID));
plotRes = textscan(fID,'%s',2,'delimiter',',','commentStyle','#');
fclose all;
tempLine = char(plotRes{1,1}(1));
plotYes = false;
writeYes = false;
if strcmpi(tempLine(1),'y')
    plotYes = true;
end
tempLine = char(plotRes{1,1}(2));
if strcmpi(tempLine(1),'y')
    writeYes = true;
end
    % function for getting number for single line entry. Ignores spaces
    function [valOut] = subString(stringIn)
        val = regexp(stringIn,'#','split');
        val = regexprep(regexprep(val{1},' ',''),'\t','');
        if isempty(val)
            valOut = NaN;
        else
            valOut = str2double(val);
        end
    end
end
