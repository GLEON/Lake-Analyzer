function Run_LA( LakeName,Folder,skipLoad)
%----Author: Jordan S Read 2009 ----

if nargin < 3
    skipLoad = false;
end
clc; close all
done = false;
if ~skipLoad
    build_config(LakeName,Folder);
    while ~done
        pause(.4)
    end
    pause(0.1);
end

%% -- variables --
matSec = 86400;         % number of seconds in a day
smplTs = 50;            % samplerate test length
dateTl = 0.00001;       % tolerance for equal day (day fraction)
fprintf(['Reading ' LakeName '.lke file'])
[outPuts,outRs,maxZ,wndH,wndAv,lyrAv,outWn,wtrMx,wtrMn,...
    wndMx,wndMn,drhDz,Tdiff,plotYes,writeYes] = ...
    OpenCfg( LakeName,Folder );
Smin = 0.1;
if lt(drhDz,Smin)
    Smin = drhDz;
end
%% -- end variables --

if ~any([plotYes writeYes])
    
    error(['User must specify either to write results,'...
        ' plot results, or both ....... (see ' LakeName '.lke)'])

end

fprintf('...completed\n\n') ;

fprintf('****Building program structure****\n') ;

[TT,outputOptions,writeTable,plotTable,dateInput,dateOutput,delimI,delimO] = ...
    OutputConstructor(outPuts);
TT.useLvl = false;
fNms = fieldnames(TT);
for j = 1:length(fNms)
    if TT.(char(fNms{j}))
       fprintf([char(fNms{j}) '\n'])
    end
end
fprintf('****completed****\n\n') ;

if TT.openWtr
    TT.useSal = true;       % salinity use is default true
    fprintf(['Reading ' LakeName '.wtr file'])
    TT.errCkWtr = true;
    if isinf(wtrMx) && isinf(wtrMn)
        TT.errCkWtr = false;
    elseif isnan(wtrMx)||isnan(wtrMn)
        TT.errCkWtr = false;
    end
    wtrFileName  = [Folder '/' LakeName '.wtr'];
    oper = fopen(wtrFileName);
    if eq(oper,-1)
        error([Folder '/' LakeName '.wtr file not found']);
    end
    fclose all;
    [wtrD,wtr,wtrHead] = gFileOpen(wtrFileName);
    [wtrLength,numDepths] = size(wtr);
    depthAry = NaN(1,numDepths);
    wtrHead = textscan(wtrHead,'%s',numDepths+1,'delimiter',delimI);
    for i = 1:numDepths
        tempStr = char(wtrHead{1}(i+1));
        tempLen = length(tempStr);
        depthAry(i) = str2double(tempStr((5:tempLen))); 
    end
    if any(isnan(depthAry))
            error(['water temperature headers are not in the right format'...
                ' (e.g. temp0.5 for 0.5m thermistor)']);
    end
    fprintf('...completed\n\n') ;
end
if TT.openWtr
    % check for salinity file
    defHead = 'salinity';
    fprintf(['Checking for ' LakeName '.sal file'])
    salFileName  = [Folder '/' LakeName '.sal'];
    oper = fopen(salFileName);
    if eq(oper,-1)
        TT.useSal = false;      % salinity not used if file does exist
        fprintf('...not found\n\n')
    else
        fprintf('...file found')
        fprintf(['\n' LakeName ...
            ' salinity will be used for density calculations\n\n'])
        fclose all;
        [salD,sal,salHead] = gFileOpen(salFileName);
        [salLength,numSal] = size(sal);
        salDepths = NaN(1,numSal);
        salHead = textscan(salHead,'%s',numSal+1,'delimiter',delimI);
        for i = 1:numSal
            tempStr = char(salHead{1}(i+1));
            tempLen = length(tempStr);
            salDepths(i) = str2double(tempStr((length(defHead)+1:tempLen)));
        end
        if any(isnan(depthAry))
            error(['salinity headers are not in the right format'...
                ' (e.g. ' defHead '0.5 for 0.5m measurement)']);
        end
    end
    
    
end
if TT.openWnd
    fprintf(['Reading ' LakeName '.wnd file'])
    TT.errCkWnd = true;
    if isinf(wndMx) && isinf(wndMn)
        TT.errCkWnd = false;
    elseif isnan(wndMx)||isnan(wndMn)
        TT.errCkWnd = false;
    end
    wndFileName  = [Folder '/' LakeName '.wnd'];
    oper = fopen(wndFileName);
    if eq(oper,-1)
        error([Folder '/' LakeName '.wnd file not found']);
    end
    fclose all;
    
    [wndD,wnd] = gFileOpen(wndFileName);
    wndLength = length(wndD);
    fprintf('...completed\n\n') ;
end   
if TT.openBth
    fprintf(['Reading ' LakeName '.bth file'])
    oper = fopen([LakeName '.bth']);
    if eq(oper,-1)
        % check other directory
        oper = fopen([Folder '/' LakeName '.bth']);
        if eq(oper,-1)
            error([LakeName '.bth file not found']);
        end
        fclose all;
        bathArray = fileOpen([Folder '/' LakeName '.bth'],',',1);
    else
    fclose all;
    bathArray = fileOpen([LakeName '.bth'],',',1);
    end
    bthD(1,:)   = bathArray(:,1);
    bthA(1,:)   = bathArray(:,2);
    clear bathArray
    fprintf('...completed\n\n') ;
    % check for .lvl file here...
    fprintf(['Checking for ' LakeName '.lvl file'])
    lvlFileName  = [Folder '/' LakeName '.lvl'];
    oper = fopen(lvlFileName);
    if eq(oper,-1)
        TT.useSal = false;      % salinity not used if file does exist
        fprintf('...not found\n\n')
    else
        fprintf('...file found')
        fprintf(['\n' LakeName ...
            ' lake level file will be used for all hypsographically relevant functions\n'])
        fprintf(['\n' LakeName ...
            ' lake level is based on 0 depth of the *.bth file, '...
            '\nwith positive z in the downward direction.\n\n'])
        TT.useLvl = true;
        fclose all;
        [lvlD,lvl] = gFileOpen(lvlFileName);
        if any(lt(lvl,min(bthD)))   % test for values that are above the bth min
            dateOver = lvlD(lt(lvl,min(bthD)));
            error([LakeName '.lvl file contains values above the maximum '...
                'elevation of the *.bth file (i.e. values above ' ...
                num2str(min(bthD)) ...
                ') for reference, the first observation above the maximum *.bth '...
                'data occurred on ' datestr(dateOver(1),'yyyy-mm-dd')]);
        elseif any(ge(lvl,max(bthD)))   % test for values that are below the bth min
            dateOver = lvlD(ge(lvl,max(bthD)));
            error([LakeName '.lvl file contains values below the minimum '...
                'elevation of the *.bth file (i.e. values below ' ...
                num2str(max(bthD)) ...
                ') for reference, the first observation below the minimum *.bth '...
                'data occurred on ' datestr(dateOver(1),'yyyy-mm-dd')]);
        end
    end
end


%*** find samplerate of raw data *** (has to be at least wind or water)
if TT.openWtr
    if length(wtrD) < smplTs
        tLen = length(wtrD);
    else
        tLen = smplTs;
    end
    
    steps = NaN(1,tLen-1);
    for i = 1:tLen-1
        steps(i) = wtrD(i+1)-wtrD(i);
    end
else
    if length(wndD) < smplTs
        tLen = length(wndD);
    else
        tLen = smplTs;
    end
    
    steps = NaN(1,tLen-1);
    for i = 1:tLen-1
        steps(i) = wndD(i+1)-wndD(i);
    end
end
if eq(min(steps),0)
    matRs = mean(steps)*matSec;
else
    matRs = min(steps)*matSec; %current sample rate of raw data in seconds
end
clear vals ind numMx numCont steps tLen
if (outRs - matRs)>dateTl
    TT.dwnSmple = true;   %down sample if necessary
end

outWn = ceil(outWn/matRs); % *** find down sample window ***


% *** error checking *** // this is slow...NEED TO CHECK END OF FILE TOO...
if TT.errCkWtr
    for j = 1:ceil(outWn):wtrLength-outWn
        for i = 1:numDepths
            outInd = FindOutliers_TS(wtr(j:j+outWn,i),wtrMx,wtrMn);
            goodInd = j:j+outWn;
            goodInd(outInd) = [];
            outInd = outInd+j-1;
            if outInd > 0
                wtr(outInd,i) = mean(wtr(goodInd,i));
            end
        end
    end   
end
if TT.errCkWnd %make sure this is only true when useWnd is true
    indLow = find(wnd<wndMn);
    wnd(indLow) = ones(1,length(wnd(indLow)))*wndMn;
    for j = 1:ceil(outWn):wndLength-outWn
        outInd = FindOutliers_TS(wnd(j:j+outWn),wndMx,wndMn,2.5);
        goodInd = j:j+outWn;
        goodInd(outInd) = [];
        outInd = outInd+j-1;
        if outInd > 0
            wnd(outInd) = mean(wnd(goodInd));
        end
    end 
    clear goodInd outInd
end

% *** down sampling ***
if TT.dwnSmple
    fprintf('Down sampling data');
    if TT.openWtr
        for i = 1:numDepths
            [ DS_wtrD,DS_wtr(:,i) ] = DownSample_TS(wtrD,outRs,wtr(:,i));
        end
        wtrD = DS_wtrD;
        wtr = DS_wtr;
        clear DS_wtrD DS_wtr
    end
    if TT.openWnd
        if TT.openWnd
            [wndD,wnd] = DownSample_TS(wndD,outRs,wnd);
            if TT.openWtr
                [wtrMI,wndMI] = findMatches(wtrD,wndD);
                wndD = wndD(wndMI);
                wtrD = wtrD(wtrMI);
                wtr  = wtr(wtrMI,:);
                wnd  = wnd(wndMI);
            end
        end
        wndLength = length(wnd);
        clear wndMI wtrMI

        %at this point, use wind averaging
        %wind averaging will always be used (if not needed, specify a
        %window that is smaller than the downsampled resolution 
        %(i.e. outRs > wndAv)
        % *** this is backwards looking averging ***
        
        wndAv = ceil(wndAv/outRs);
        Awnd = NaN(wndLength,1);
        if isempty(wnd)
            error(['no matching wind and water temp dates - '...
                'check input files for formatting errors'])
        end
        for j = 1:wndAv
            Awnd(j) = mean(wnd(1:j));
        end
        for j = wndAv+1:wndLength
            Awnd(j) = mean(wnd(j-wndAv:j)); 
        end
        wnd = Awnd; 
        clear wndLength
    end
    if TT.openWtr
        varL = length(wtrD);
        dates = wtrD;
    else
        varL = length(wndD);
        dates = wndD;
    end
    fprintf('...completed\n\n');
else varL = length(wtrD);
end
% -- introduce salinity --
if TT.useSal
    fprintf('Introducing salinity')
    % extend minimum measurement to minimum thermistor
    % assumes depths are in order...
    salMn = min(salDepths);
    salMx = max(salDepths);
    thmMn = min(depthAry);
    thmMx = max(depthAry);
    if gt(thmMn,salMx) || lt(thmMx,salMn)
        error('no matching salinity dates')
    end
    if lt(thmMn,salMn)      % if the thermistors are closer to the surface
        fprintf(['Extending ' LakeName ' ' num2str(salDepths(1)) 'm'...
            ' measurement to ' num2str(thmMn) 'm\n']);
        sal = [sal(:,1) sal];
        salDepths = [thmMn salDepths];
        
    end
    if gt(thmMx,salMx)      % if the thermistors extend deeper than sal
        fprintf(['Extending ' LakeName ' ' num2str(salDepths(end)) 'm'...
            ' measurement to ' num2str(thmMx) 'm\n']);
        sal = [sal sal(:,end)];
        salDepths= [salDepths thmMx];
    end
    % now correct for time
    salSt = min(salD);
    salEn = max(salD);
    thmSt = min(wtrD);
    thmEn = max(wtrD);
    if lt(thmSt,salSt)      % if the thermal measuremets start earlier
        fprintf(['Extending ' LakeName ' ' datestr(salSt,'dd-mmm-yyyy HH:MM')...
            ' salinity measurement to ' datestr(thmSt,'dd-mmm-yyyy HH:MM') '\n']);
        sal = [sal(1,:); sal];
        salD = [thmSt; salD];
    end
    if gt(thmEn,salEn)      % if the thermal measuremets end later
        fprintf(['Extending ' LakeName ' ' datestr(salEn,'dd-mmm-yyyy HH:MM')...
            ' salinity measurement to ' datestr(thmEn,'dd-mmm-yyyy HH:MM') '\n']);
        sal = [sal; sal(end,:)];
        salD = [salD; thmEn];
    end
   % baseSal = sal;
   % sal = wtr*NaN;
    
    fprintf('\nInterpolating salinity to match thermistors\n\n')
    % first in Z direction:
    sal = interp1(salDepths',sal',depthAry)';   % interpolated in depth
    sal = interp1(salD,sal,wtrD);               % interpolated in time;
    
else
    sal = wtr*0;    % if salinity is not to be used, assumed zero
end



%% find layers
if TT.findLyr
    fprintf('Finding thermal layers');
    mixed = zeros(varL,1);
    thermoD = ones(varL,1)*depthAry(numDepths);metaT = thermoD;
    metaB = metaT;thermoInd = ones(varL,1);
    if TT.ssnLyr
        SthermoD = thermoD;SmetaT = thermoD;
        SmetaB = metaT;SthermoInd = thermoInd;
    end
    rho = zeros(varL,numDepths);
    for j = 1:varL
        rho(j,:) = waterDensity(wtr(j,:),sal(j,:));
        wtrT = wtr(j,:);
        % test shallowest depth with deepest depth (exclude NaNs)
        wtrT = wtrT(~isnan(wtrT));
        if abs(wtrT(1)-wtrT(end)) > Tdiff % not mixed... % GIVES mixed if NaN!!!!
            % remove NaNs, need at least 3 values
                rhoT = rho(j,:); 
                depT = depthAry; depT(isnan(rhoT)) = [];
                rhoT(isnan(rhoT)) = [];
                
            if TT.ssnLyr
                if length(depT)>2
                    
                    [thermoD(j),thermoInd(j),drho_dz,SthermoD(j),SthermoInd(j)] = ...
                        FindThermoDepth( rhoT,depT,Smin );
                	metaT(j)  = FindMetaTop(drho_dz,thermoD(j),depT,drhDz);
                    metaB(j)  = FindMetaBot(drho_dz,thermoD(j),depT,drhDz);
                    SmetaT(j) = FindMetaTop(drho_dz,SthermoD(j),depT,drhDz);
                    SmetaB(j) = FindMetaBot(drho_dz,SthermoD(j),depT,drhDz);
                end % or else, keep as NaN

            else
                if length(depT)>2
                    [thermoD(j),thermoInd(j),drho_dz] = ...
                        FindThermoDepth(rhoT,depT,Smin);
                    metaT(j) = FindMetaTop(drho_dz,thermoD(j),depT,drhDz);
                    metaB(j) = FindMetaBot(drho_dz,thermoD(j),depT,drhDz);
                end % or else, keep as NaN
            end
        else
            mixed(j) = 1;
        end
    end
    
    %layer averaging ***
    lyrAv = ceil(lyrAv/outRs);
    if TT.ssnLyr
        AsMetaB  = zeros(varL,1);
        AsMetaT  = zeros(varL,1);
        AsThermoD= zeros(varL,1);
    end
    AvMetaB  = zeros(varL,1);
    AvMetaT  = zeros(varL,1);
    AvThermoD= zeros(varL,1);
        %fill pre window size
    for j = 1:lyrAv
        AvMetaB(j) = mean(metaB(1:j));
        AvMetaT(j) = mean(metaT(1:j));
        AvThermoD(j) = mean(thermoD(1:j));
        if TT.ssnLyr
            AsMetaB(j) = mean(SmetaB(1:j));
            AsMetaT(j) = mean(SmetaT(1:j));
            AsThermoD(j) = mean(SthermoD(1:j));
        end
    end
    for j = lyrAv+1:varL
        AvMetaB(j) = mean(metaB(j-lyrAv+1:j));
        AvMetaT(j) = mean(metaT(j-lyrAv+1:j));
        AvThermoD(j) = mean(thermoD(j-lyrAv+1:j));
        if TT.ssnLyr
            AsMetaB(j) = mean(SmetaB(j-lyrAv+1:j));
            AsMetaT(j) = mean(SmetaT(j-lyrAv+1:j));
            AsThermoD(j) = mean(SthermoD(j-lyrAv+1:j));
        end        
    end
    if TT.ssnLyr
        SmetaB = AsMetaB;
        SmetaT = AsMetaT;
        SthermoD = AsThermoD;
        clear AsThermoD AsMetaT AsMetaB
    end
    
    metaB   = AvMetaB;
    metaT   = AvMetaT;
    thermoD = AvThermoD;
    clear AvMetaB AvMetaT AvThermoD
    
    fprintf('...completed\n\n');
end
%% make sure wind and temp dates match
if TT.openWnd && TT.openWtr
    [wtrMI,wndMI] = findMatches(wtrD,wndD);
    dates = wtrD(wtrMI);
    wtr  = wtr(wtrMI,:);
    wnd  = wnd(wndMI);
    [varL,nmD] = size(wtr);
    sal = sal(wtrMI,:);
elseif TT.openWnd
    dates = wndD;
    varL = length(dates);
    wtrMI = 1:varL;
elseif TT.openWtr
    dates = wtrD;
    varL = length(dates);
    wtrMI = 1:varL;
end
%% - interpolate level
if TT.useLvl
    lvlSt = min(lvlD);
    lvlEn = max(lvlD);
    wtrSt = min(dates);
    wtrEn = max(dates);
    if gt(wtrSt,lvlEn) || lt(wtrEn,lvlSt)
        error('no matching level dates')
    end
    if gt(lvlSt,wtrSt)      % if the thermal measuremets start earlier
        fprintf(['Extending ' LakeName ' ' datestr(lvlSt,'dd-mmm-yyyy HH:MM')...
            ' level measurement to ' datestr(wtrSt,'dd-mmm-yyyy HH:MM') '\n']);
        lvl = [lvl(1); lvl];
        lvlD = [wtrSt; lvlD];
    end
    if lt(lvlEn,wtrEn)      % if the thermal measuremets end later
        fprintf(['Extending ' LakeName ' ' datestr(lvlEn,'dd-mmm-yyyy HH:MM')...
            ' level measurement to ' datestr(wtrEn,'dd-mmm-yyyy HH:MM') '\n']);
        lvl = [lvl; lvl(end)];
        lvlD = [lvlD; wtrEn];
    end
    if any(isnan(lvl))
        rmvI = isnan(lvl);
        fprintf(['Removing ' LakeName ' ' datestr(lvlD(rmvI),'dd-mmm-yyyy HH:MM')...
            ' errant level measurements\n']);
    end
    fprintf('\nInterpolating level to match thermistors\n\n')
    lvl = interp1(lvlD,lvl,dates);
    lvl  = repmat(lvl,1,length(bthD));
    bthD = repmat(bthD,varL,1);
    bthD = bthD-lvl;    % now bthD is a time series, bthA stays relative to initial measurements
end


% - apply to layers -
if TT.findLyr
    metaB = metaB(wtrMI);
    metaT = metaT(wtrMI);
    thermoD = thermoD(wtrMI);
    thermoInd = thermoInd(wtrMI);
    if TT.useLvl
        mxZ   = bthD(:,end);
        if any(gt(metaB,mxZ))   % depth of indices is lower than max depth
            fprintf('One or more metaB depths appear below the bottom stage, replacing with bottom\n');
            rplcI = gt(metaB,mxZ);
            metaB(rplcI) = mxZ(rplcI);
        end
        if any(gt(metaT,mxZ))   % depth of indices is lower than max depth
            fprintf('One or more metaT depths appear below the bottom stage, replacing with bottom\n');
            rplcI = gt(metaT,mxZ);
            metaT(rplcI) = mxZ(rplcI);
        end
        if any(gt(thermoD,mxZ))   % depth of indices is lower than max depth
            fprintf('One or more thermocline depths appear below the bottom stage, replacing with bottom\n');
            rplcI = gt(thermoD,mxZ);
            thermoD(rplcI) = mxZ(rplcI);
        end
        % --- check thermoInd here --- incomplete
    end
        
end
if TT.ssnLyr
    SmetaB = SmetaB(wtrMI);
    SmetaT = SmetaT(wtrMI);
    SthermoD = SthermoD(wtrMI);
    SthermoInd = SthermoInd(wtrMI);
    if TT.useLvl
        mxZ   = bthD(:,end);
        if any(gt(SmetaB,mxZ))   % depth of indices is lower than max depth
            fprintf('One or more SmetaB depths appear below the bottom stage, replacing with bottom\n');
            rplcI = gt(SmetaB,mxZ);
            SmetaB(rplcI) = mxZ(rplcI);
        end
        if any(gt(SmetaT,mxZ))   % depth of indices is lower than max depth
            fprintf('One or more SmetaT depths appear below the bottom stage, replacing with bottom\n');
            rplcI = gt(SmetaT,mxZ);
            SmetaT(rplcI) = mxZ(rplcI);
        end
        if any(gt(SthermoD,mxZ))   % depth of indices is lower than max depth
            fprintf('One or more Sthermocline depths appear below the bottom stage, replacing with bottom\n');
            rplcI = gt(SthermoD,mxZ);
            SthermoD(rplcI) = mxZ(rplcI);
        end
        % --- check thermoInd here --- incomplete
    end
end

%****-----varL is the length of output files as of here-------*****

%
%% *** schmidt stability ***
if TT.StYes
    fprintf('Calculating Scmidt Stability');
    St = NaN(varL,1);
    for j = 1:varL
        wtrT = wtr(j,:);
        salT = sal(j,:);
        depT = depthAry; 
        depT(isnan(wtrT)) = [];
        salT(isnan(wtrT)) = [];
        wtrT(isnan(wtrT)) = [];
            
        if length(wtrT) > 2
            if TT.useLvl
                St(j) = schmidtStability(wtrT,depT,bthA,bthD(j,:),salT);
            else
                St(j) = schmidtStability(wtrT,depT,bthA,bthD,salT);
            end
        end % else keep as NaN
    end
    if writeTable.St
        writeTable.St = St;
    end
    fprintf('...completed\n\n');
end

%% *** uStar ***
if TT.uStYes
    fprintf('Calculating U-star');
    uSt = NaN(varL,1);
    minDepth = depthAry(1);
    for j = 1:varL
        wtrT = wtr(j,:);
        salT = sal(j,:);
        depT = depthAry;
        depT(isnan(wtrT)) = [];
        salT(isnan(wtrT)) = [];
        wtrT(isnan(wtrT)) = [];
        if ~any(isnan([minDepth metaT(j)])) && ~isempty(wtrT)
            if TT.useLvl
                AvEp_rho = layerDensity(minDepth,metaT(j),wtrT,depT,...
                    bthA,bthD(j,:),salT);
            else
                AvEp_rho = layerDensity(minDepth,metaT(j),wtrT,depT,...
                    bthA,bthD,salT); 
            end
            uSt(j) = uStar(wnd(j),wndH,AvEp_rho);
        end % else, keep as NaN
    end
    if writeTable.uSt
        writeTable.uSt = uSt;
    end
    fprintf('...completed\n\n');
end
    
%% *** lake number ***
if TT.wrt_Ln
    fprintf('Calculating Lake Number');
    Zm = max(bthD);
    LN = NaN(varL,1);
    for j = 1:varL
        wtrT = wtr(j,:);
        salT = sal(j,:);
        depT = depthAry;
        depT(isnan(wtrT)) = [];
        salT(isnan(wtrT)) = [];
        wtrT(isnan(wtrT)) = [];
        if ~any(isnan([Zm metaB(j) St(j) uSt(j)]))
            if TT.useLvl
                Zm = max(bthD(j,:));
                AvHyp_rho = layerDensity(metaB(j),Zm,wtrT,depT,...
                    bthA,bthD(j,:),salT);
                LN(j) = lakeNumber(bthA,bthD(j,:),uSt(j),...
                    St(j),metaT(j),metaB(j),AvHyp_rho);
            else
                AvHyp_rho = layerDensity(metaB(j),Zm,wtrT,depT,...
                    bthA,bthD,salT);
                LN(j) = lakeNumber(bthA,bthD,uSt(j),...
                    St(j),metaT(j),metaB(j),AvHyp_rho);
            end
            
        end
    end
    writeTable.Ln = LN;
    fprintf('...completed\n\n');
end
%% *** wedderburn number ***
if TT.wrt_W
    fprintf('Calculating Wedderburn Number');
    W = NaN(varL,1);
    [Zo,Io] = min(bthD);
    Zm = max(bthD);
    if ~TT.useLvl
        Ao = bthA(Io);
    end
    for j = 1:varL
        wtrT = wtr(j,:);
        salT = sal(j,:);
        depT = depthAry;
        depT(isnan(wtrT)) = [];
        salT(isnan(wtrT)) = [];
        wtrT(isnan(wtrT)) = [];
        if ~any(isnan([Zm minDepth metaT(j) metaB(j) uSt(j)]))
            if TT.useLvl
                [Zo,Io] = min(bthD(j,:));
                Zm = max(bthD(j,:));
                Ao = bthA(Io);
                AvEp_rho = layerDensity(minDepth,metaT(j),wtrT,depT,...
                    bthA,bthD(j,:),salT);
                AvHyp_rho = layerDensity(metaB(j),Zm,wtrT,depT,...
                    bthA,bthD(j,:),salT);
                del_rho = AvHyp_rho-AvEp_rho;
                W(j) = wedderburnNumber(del_rho,metaT(j),uSt(j),Ao,...
                    AvHyp_rho);
            else

                AvEp_rho = layerDensity(minDepth,metaT(j),wtrT,depT,...
                    bthA,bthD,salT);
                AvHyp_rho = layerDensity(metaB(j),Zm,wtrT,depT,...
                    bthA,bthD,salT);
                del_rho = AvHyp_rho-AvEp_rho;
                W(j) = wedderburnNumber(del_rho,metaT(j),uSt(j),Ao,...
                    AvHyp_rho);
            end
            
        end
    end
    writeTable.W = W;
    fprintf('...completed\n\n');
end

%% *** wind speed ***
if TT.wrt_wndSpd
    writeTable.wndSpd = wnd;
end

%% *** metalimnion top ***
if TT.wrt_metaT
    writeTable.metaT = metaT;
end

%% *** metalimnion bottom ***
if TT.wrt_metaB
    writeTable.metaB = metaB;
end

%% *** thermocline depth ***
if TT.wrt_thermD
    writeTable.thermD = thermoD;
end

%% *** seasonal thermocline depth ***
if TT.wrt_SthermD
    writeTable.SthermD = SthermoD;
end

%% *** seasonal metalimnion bottom ***
if TT.wrt_SmetaB
    writeTable.SmetaB = SmetaB;
end

%% *** seasonal metalimnion top ***
if TT.wrt_SmetaT
    writeTable.SmetaT = SmetaT;
end

%% *** seasonal uStar ***
if TT.SuStYes
    fprintf('Calculating Parent U-star');
    SuSt = NaN(varL,1);
    minDepth = depthAry(1);
    for j = 1:varL
        wtrT = wtr(j,:);
        salT = sal(j,:);
        depT = depthAry;
        depT(isnan(wtrT)) = [];
        salT(isnan(wtrT)) = [];
        wtrT(isnan(wtrT)) = [];
        if ~any(isnan([minDepth SmetaT(j)])) && ~isempty(wtrT)
            if TT.useLvl
                AvEp_rho = layerDensity(minDepth,SmetaT(j),wtrT,depT,...
                    bthA,bthD(j,:),salT);
            else
                AvEp_rho = layerDensity(minDepth,SmetaT(j),wtrT,depT,...
                    bthA,bthD,salT);
            end
            
            SuSt(j) = uStar(wnd(j),wndH,AvEp_rho);
        end % else, keep as NaN
    end
    if writeTable.SuSt
        writeTable.SuSt = SuSt;
    end
    fprintf('...completed\n\n');
end

%% *** seasonal lake number ***
if TT.wrt_SLn
    fprintf('Calculating Parent Lake Number');
    SLN = NaN(varL,1);
    Zm = max(bthD);
    for j = 1:varL
        wtrT = wtr(j,:);
        salT = sal(j,:);
        depT = depthAry;
        depT(isnan(wtrT)) = [];
        salT(isnan(wtrT)) = [];
        wtrT(isnan(wtrT)) = [];
        if ~any(isnan([Zm SmetaB(j) St(j) SuSt(j)])) && ~isempty(wtrT)
            if TT.useLvl
                Zm = max(bthD(j,:));
                AvHyp_rho = layerDensity(SmetaB(j),Zm,wtrT,depT,...
                    bthA,bthD(j,:),salT);
                SLN(j) = lakeNumber(bthA,bthD(j,:),SuSt(j),...
                    St(j),SmetaT(j),SmetaB(j),AvHyp_rho);
            else
                AvHyp_rho = layerDensity(SmetaB(j),Zm,wtrT,depT,...
                    bthA,bthD,salT);
                SLN(j) = lakeNumber(bthA,bthD,SuSt(j),...
                    St(j),SmetaT(j),SmetaB(j),AvHyp_rho);
            end
            
        end
    end
    writeTable.SLn = SLN;
    fprintf('...completed\n\n');
end

%% *** seasonal Wedderburn number ***
if TT.wrt_SW
    fprintf('Calculating Parent Wedderburn Number');
    SW = NaN(varL,1);
    [Zo,Io] = min(bthD);
    Zm = max(bthD);
    if ~TT.useLvl
        Ao = bthA(Io);
    end
    for j = 1:varL
        wtrT = wtr(j,:);
        salT = sal(j,:);
        depT = depthAry;
        depT(isnan(wtrT)) = [];
        salT(isnan(wtrT)) = [];
        wtrT(isnan(wtrT)) = [];
        if ~any(isnan([Zm SmetaT(j) SmetaB(j) minDepth SuSt(j)]))
            if TT.useLvl
                [Zo,Io] = min(bthD(j,:));
                Zm = max(bthD(j,:));
                Ao = bthA(Io);
                AvEp_rho = layerDensity(minDepth,SmetaT(j),wtrT,depT,...
                    bthA,bthD(j,:),salT);
                AvHyp_rho = layerDensity(SmetaB(j),Zm,wtrT,depT,...
                    bthA,bthD(j,:),salT);
            else
                AvEp_rho = layerDensity(minDepth,SmetaT(j),wtrT,depT,...
                    bthA,bthD,salT);
                AvHyp_rho = layerDensity(SmetaB(j),Zm,wtrT,depT,...
                    bthA,bthD,salT);
            end
            
            
            del_rho = AvHyp_rho-AvEp_rho;
            SW(j) = wedderburnNumber(del_rho,SmetaT(j),SuSt(j),Ao,...
                AvHyp_rho);
        end
    end
    writeTable.SW = SW;
    fprintf('...completed\n\n');
end

%% *** buoyancy frequency ***
if TT.wrt_N2
    fprintf('Calculating Buoyancy Frequency');
    g = 9.81; %m s-1
    N2 = NaN(varL,1);
    for j = 1:varL
        Pw = waterDensity(wtr(j,thermoInd(j)),sal(j,thermoInd(j)));
        P2 = waterDensity(wtr(j,thermoInd(j)+1),sal(j,thermoInd(j)+1));
        dw = depthAry(thermoInd(j));
        d2 = depthAry(thermoInd(j)+1);
        N2(j) = g/Pw*(P2-Pw)/(d2-dw);
    end
    writeTable.N2 = N2;
    fprintf('...completed\n\n');
end

%% *** seasonal buoyancy frequency ***
if TT.wrt_SN2
    fprintf('Calculating Parent Buoyancy Frequency');
    g = 9.81; %m s-1
    SN2 = NaN(varL,1);
    for j = 1:varL
        Pw = waterDensity(wtr(j,SthermoInd(j)),sal(j,SthermoInd(j)));
        P2 = waterDensity(wtr(j,SthermoInd(j)+1),sal(j,SthermoInd(j)+1));
        dw = depthAry(SthermoInd(j));
        d2 = depthAry(SthermoInd(j)+1);
        SN2(j) = g/Pw*(P2-Pw)/(d2-dw);
    end
    writeTable.SN2 = SN2;
    fprintf('...completed\n\n');
end

%% *** mode 1 vertical seiche 
if TT.wrt_T1
    fprintf('Calculating T1');
    g = 9.81; %m/s2
    T1 = NaN(varL,1);
    Zm = max(bthD);
    for j = 1:varL
        wtrT = wtr(j,:);
        salT = sal(j,:);
        depT = depthAry;
        depT(isnan(wtrT)) = [];
        salT(isnan(wtrT)) = [];
        wtrT(isnan(wtrT)) = [];
        if TT.useLvl
            Zm = max(bthD(j,:));
            A = interp1(bthD(j,:),bthA,thermoD(j));
            L = 2*sqrt(A/pi);
            AvEpi_rho = layerDensity(0,thermoD(j),wtrT,depT,...
                bthA,bthD(j,:),salT);
            AvHyp_rho = layerDensity(thermoD(j),Zm,wtrT,depT,...
                bthA,bthD(j,:),salT);
        else
            A = interp1(bthD,bthA,thermoD(j));
            L = 2*sqrt(A/pi);
            AvEpi_rho = layerDensity(0,thermoD(j),wtrT,depT,...
                bthA,bthD,salT);
            AvHyp_rho = layerDensity(thermoD(j),Zm,wtrT,depT,...
                bthA,bthD,salT);
        end
        
        delta_rho = AvHyp_rho-AvEpi_rho;
        go = g*delta_rho/AvHyp_rho;
        if lt(abs(wtr(j,1)-wtr(j,numDepths)),Tdiff)
            T1(j) = NaN;
        else
            T1(j) = 2*L/(go*thermoD(j)*(Zm-thermoD(j))/Zm);
        end
        if lt(T1(j),0)
            T1(j) = 0;
        end
    end
    writeTable.T1 = T1;
    fprintf('...completed\n\n');
end

%% *** Seasonal mode 1 vertical seiche 
if TT.wrt_ST1
    fprintf('Calculating Parent T1');
    g = 9.81; %m/s2
    ST1 = NaN(varL,1);
    Zm = max(bthD);
    for j = 1:varL
        wtrT = wtr(j,:);
        salT = sal(j,:);
        depT = depthAry;
        depT(isnan(wtrT)) = [];
        salT(isnan(wtrT)) = [];
        wtrT(isnan(wtrT)) = [];
        if TT.useLvl
            Zm = max(bthD(j,:));
            A = interp1(bthD(j,:),bthA,SthermoD(j));
            L = 2*sqrt(A/pi);
            AvEpi_rho = layerDensity(0,SthermoD(j),wtrT,depT,...
                bthA,bthD(j,:),salT);
            AvHyp_rho = layerDensity(SthermoD(j),Zm,wtrT,depT,...
                bthA,bthD(j,:),salT);
        else
            A = interp1(bthD(j,:),bthA,SthermoD(j));
            L = 2*sqrt(A/pi);
            AvEpi_rho = layerDensity(0,SthermoD(j),wtrT,depT,...
                bthA,bthD,salT);
            AvHyp_rho = layerDensity(SthermoD(j),Zm,wtrT,depT,...
                bthA,bthD,salT);
        end
        
        delta_rho = AvHyp_rho-AvEpi_rho;
        go = g*delta_rho/AvHyp_rho;
        if lt(abs(wtr(j,1)-wtr(j,numDepths)),Tdiff)
            ST1(j) = NaN;
        else
            ST1(j) = 2*L/(go*SthermoD(j)*(Zm-SthermoD(j))/Zm);
        end
        if lt(ST1(j),0)
            ST1(j) = 0;
        end
    end
    writeTable.ST1 = ST1;
    fprintf('...completed\n\n');
end

%% alter wtr values according to lvl file (if it is used)
if TT.useLvl
    for j = 1:varL
        nanI = gt(depthAry,bthD(j,end));
        wtr(j,nanI) = NaN;
    end
end
%% build plot array
if plotYes
    fprintf('Plotting results');
    if writeTable.wTemp
        plotLA_results(writeTable,plotTable,dates,LakeName,Folder,wtr,depthAry)
    else
        plotLA_results(writeTable,plotTable,dates,LakeName,Folder)
    end
    fprintf('...completed\n\n');
end

%% build file array
writeNames = {};
cnt = 1;
for k = 1:length(outputOptions)
    if ~islogical(writeTable.(char(outputOptions{k})))
        writeNames{cnt} = outputOptions{k};
        cnt = cnt+1;
    end
end

%% write to file
if writeYes
    fprintf('Writing results to file');
end
if writeYes && ~isempty(writeNames)
    
    
    outputFile = [Folder '/' LakeName '_results.txt'];
    outFile = fopen(outputFile,'w');
    if eq(outFile,-1)
        error([Folder '/' LakeName '_results.csv file in use, please close']);
    end
    wrt = @(writer)fprintf(outFile,writer); % build a subfunction that writes 
                                    % the contents of the input "writer" 
                                    % to the file everytime wrt is called
    wrt('DateTime');
    for i = 1:cnt-1
        wrt([delimO writeNames{i}]);
    end
    
    wrt('\r\n');
    for j = 1:varL
        wrt(datestr(dates(j),dateOutput)); %change 'dateOutput' 
                                    % in the 'OutputConstructor.m' file
        for i = 1:length(writeNames)
            wrt([delimO num2str(writeTable.(char(writeNames{i}))(j))]);
        end
        wrt('\r\n');
    end
    fclose all;
end
 %---- write water to separate file ----
if TT.wrt_wTemp && writeYes
    outputFile = [Folder '/' LakeName '_results_wtr.txt'];
    outFile = fopen(outputFile,'w');
    if eq(outFile,-1)
        error([Folder '/' LakeName '_results_wtr.csv file in use, please close']);
    end
    wrt = @(writer)fprintf(outFile,writer);
    wrt('DateTime');
    for i = 2:length(wtrHead{1,1})
        wrt([delimO char(wtrHead{1,1}(i))]);
    end
    wrt('\r\n');
    for j = 1:varL
        wrt(datestr(dates(j),dateOutput));
        for i = 1:numDepths
            wrt([delimO num2str(wtr(j,i))]);
        end
        wrt('\r\n');
    end
    fclose all;
end
if writeYes
    fprintf('...completed\n\n');
end
disp('Lake Analyzer is complete')
%profile off
%profile viewer
end
