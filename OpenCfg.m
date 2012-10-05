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
outPuts = outPuts(chk);
inPuts = textscan(fID,'%f',12,'delimiter',',','commentStyle','#');
outRs = inPuts{1,1}(1);
maxZ  = inPuts{1,1}(2);
wndH  = inPuts{1,1}(3);
wndAv = inPuts{1,1}(4);
lyrAv = inPuts{1,1}(5);
outWn = inPuts{1,1}(6);
wtrMx = inPuts{1,1}(7);
wtrMn = inPuts{1,1}(8);
wndMx = inPuts{1,1}(9);
wndMn = inPuts{1,1}(10);
drhDz = inPuts{1,1}(11);
Tdiff = inPuts{1,1}(12);
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


