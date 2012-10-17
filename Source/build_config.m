function [done] = build_config(LakeName,directory)

done = false;
spcFrac = 2; %number of spaces = char;
num2delim = 20;
defSlt = 'wTemp';
% see if file exists

fI = fopen([directory '/' LakeName '.lke']);
if ~lt(fI,0)
    loader = 'on';
else
    loader = 'off';
end

fclose all;
% define possible outputs
outputOptions = {'Schmidt Stability','St';'U-star','uSt';...
    'Lake Number','Ln';'Wedderburn Number','W';...
    'Water temperature','wTemp';'Wind speed','wndSpd';...
    'Metalimnion top','metaT';'Metalimnion bottom','metaB';...
    'Thermocline depth','thermD';'Thermocline depth (parent)','SthermD';...
    'Metalimnion bottom (parent)','SmetaB';...
    'Metalimnion top (parent)','SmetaT';'U-star (parent)','SuSt';...
    'Lake Number (parent)','SLn';'Wedderburn Number (parent)','SW';...
    'Buoyancy frequency','N2';'Buoyancy frequency (parent)','SN2';...
    'Mode 1 seiche period','T1';'Mode 1 seiche period (parent)','ST1'};
% defaults for parameters
defaults = {'86400','';'86400','';'86400','21600';'40','-12';'98','0';...
    '0.1','0.5';'N','Y'};
% long names for outputs
names    = {'output resolution (s)','total depth (m)';...
    'wind averaging (s)','wind height (m)';...
    'layer averaging (s)','outlier window (s)';...
    'max water temp (°C)','min water temp (°C)';...
    'max wind speed (m/s)','min wind speed (m/s)';...
    'metalimnion slope (kg/m4)','mixed temp differential (°C)';...
    'plot figure (Y/N)','write results (Y/N)'};
lines = [4,5;7,6;8,9;10,11;12,13;14,15;16,17];
[outs,b] = sortrows(outputOptions); vL = length(b);
outN = cell(vL,1);outA = outN;outUn = outN;
for k = 1:length(b)
    outUn{k}= outputOptions{b(k),1};
    outN{k} = outputOptions{b(k),1};
    outA{k} = outputOptions{b(k),2};
end

indx = strcmp(defSlt,outA);

if ~any(indx)
    slt{1} = '';
else
    slt{1} = outN{indx};
    
    txt = outA{indx};
    outN(indx) = [];
end
close all
bckColor = [.85 .85 .85];
figDims = [150 150 800 440];
lM = 15;
tM = 20;
rM = 15;
bM = 15;
pnS = 12;
bgPw = 180;
bgPh = 120;
btnH = 25;
txtH = 20;
txtW = 40;
txtS = 5;
spc = 7;
radioW = 150;
radioH = 20;

bigPanels(1,:) = [lM figDims(4)-bgPh-tM bgPw bgPh];
bigPanels(2,:) = bigPanels(1,:); bigPanels(2,1) = bigPanels(2,1)+pnS+bgPw;
filePanel = [bigPanels(2,1)+bgPw+pnS*2 bM+pnS+btnH ...
    figDims(3)-bgPw*2-lM-pnS*3-rM figDims(4)-tM-pnS-btnH-bM];
btns(1,:) = [lM bigPanels(1,2)-spc-btnH bgPw btnH];
btns(2,:) = [bigPanels(2,1) bigPanels(1,2)-spc-btnH bgPw btnH];
btns(3,:) = [filePanel(1)+filePanel(3)-bgPw bM bgPw btnH];
btns(4,:) = [filePanel(1) bM radioW radioH];

numTxtIn = 14;
rowsTxt = ceil(numTxtIn/2);
for p = 1:rowsTxt
    txtBoxL(p,:) = [lM bM+(p-1)*(pnS+txtH) txtW txtH];
    
    txtBoxR(p,:) = [lM+pnS+bgPw bM+(p-1)*(pnS+txtH) txtW txtH ];
    diaL(p,:) = [lM+txtW+txtS -txtS+bM+(p-1)*(pnS+txtH) txtBoxR(p,1)-...
        lM-txtW-pnS txtH ];
    diaR(p,:) = [txtBoxR(p,1)+txtW+txtS -txtS+bM+(p-1)*(pnS+txtH) ...
        txtBoxR(p,1)-lM-txtW-pnS txtH ];
end
S.fh = figure('units','pixels',...
              'position',figDims,...
              'menubar','none',...
              'resize','off',...
              'numbertitle','off',...
              'name','Configuration File','Color',bckColor);
          movegui(S.fh,'center')
outFrame=uipanel('Parent',S.fh,'BackgroundColor',bckColor,...
        'Title','Output options','FontSize',10);
        setpixelposition(outFrame,[lM-spc btns(1,2)-spc ...
        bgPw+spc*2 bgPh+btnH+spc*3+pnS]);

slcFrame=uipanel('Parent',S.fh,'BackgroundColor',bckColor,...
        'Title','Output selections','FontSize',10);
        setpixelposition(slcFrame,[bigPanels(2,1)-spc btns(1,2)-spc ...
        bgPw+spc*2 bgPh+btnH+spc*3+pnS]);

lkeFrame=uipanel('Parent',S.fh,'BackgroundColor',bckColor,...
        'Title',[LakeName '.lke preview'],'FontSize',10);
        setpixelposition(lkeFrame,[filePanel(1)-spc filePanel(2)-spc ...
        filePanel(3)+spc*2 filePanel(4)+spc*2+pnS]);
lkeInFrame=uipanel('Parent',S.fh,'BackgroundColor',bckColor,...
        'Title','','FontSize',10);
        setpixelposition(lkeInFrame,[filePanel(1)-1 filePanel(2)-2 ...
        filePanel(3)+3 filePanel(4)+3]);
ParamFrame=uipanel('Parent',S.fh,'BackgroundColor',bckColor,...
        'Title','User parameters','FontSize',10);
        setpixelposition(ParamFrame,[btns(1,1)-spc ...
            bM-spc bgPw*2+pnS+spc*2 btns(1,2)-2*spc-2]);

S.loader = uicontrol('style','checkbox',...
    'unit','pixel',...
    'position',btns(4,:),...
    'backgroundcolor',bckColor,...
    'HorizontalAlign','left',...
    'FontSize',10,...
    'string','load from existing?',...
    'enable',loader);
        
S.ed = uicontrol('style','list',...
                 'unit','pix',...
                 'position',bigPanels(1,:),...
                 'min',0,'max',2,...
                 'backgroundcolor','w',...
                 'HorizontalAlign','left',...
                 'string',outN);
             
S.ls = uicontrol('style','list',...
                 'units','pix',...
                 'min',0,'max',2,...
                 'position',bigPanels(2,:),...
                 'backgroundcolor','w',...
                 'HorizontalAlign','left',...
                 'string',slt);

S.file = uicontrol('style','text',...
                 'units','pix',...
                 'min',0,'max',2,...
                 'position',filePanel,...
                 'backgroundcolor','w',...
                 'HorizontalAlign','left','string',{...
                 ['Configuration file for ' LakeName];...
                 '';...
                 [txt ' '*ones(1,9) '#outputs'];...
                 [defaults{1,1} ' '*ones(1,11) '#output resolution (s)'];...
                 ['??' ' '*ones(1,17) '#total depth (m)'];...
                 ['??' ' '*ones(1,17) '#height from surface for wind measurement (m)'];...
                 [defaults{2,1} ' '*ones(1,11) '#wind averaging (s)'];...
                 [defaults{3,1} ' '*ones(1,11) '#thermal layer averaging (s)'];...
                 [defaults{3,2} ' '*ones(1,11) '#outlier window (s)'];...
                 [defaults{4,1} ' '*ones(1,17) '#max water temp (°C)    inf if none'];...
                 [defaults{4,2} ' '*ones(1,16) '#min water temp (°C)    -inf if none'];...
                 [defaults{5,1} ' '*ones(1,17) '#max wind speed (m/s)   inf if none'];...
                 [defaults{5,2} ' '*ones(1,19) '#min wind speed (m/s)   -inf if none'];...
                 [defaults{6,1} ' '*ones(1,14) '#meta min slope (drho/dz per m)'];...
                 [defaults{6,2} ' '*ones(1,16) '#mixed temp differential (°C)'];...
                 [defaults{7,1} ' '*ones(1,19) '#plot figure (Y/N)'];...
                 [defaults{7,2} ' '*ones(1,18) '#write results to file (Y/N)']});
             % ADD button
S.pb = uicontrol('style','push',...
                 'units','pix',...
                 'position',btns(1,:),...
                 'HorizontalAlign','center',...
                 'string','Add',...
                 'fontsize',12,'fontweight','bold',...
                 'callback',{@pb_call,S});
             
S.re = uicontrol('style','push',...
                 'units','pix',...
                 'position',btns(2,:),...
                 'HorizontalAlign','center',...
                 'string','Remove',...
                 'fontsize',12,'fontweight','bold',...
                 'callback',{@re_call,S});
             
S.publish = uicontrol('style','push',...
                 'units','pix',...
                 'position',btns(3,:),...
                 'HorizontalAlign','center',...
                 'string','Publish',...
                 'fontsize',12,'fontweight','bold',...
                 'callback',{@publish,S});

for lm = 1:rowsTxt
    indr = rowsTxt-lm+1;
S.outT(1,indr) = uicontrol('style','edit',...
    'units','pix',...
    'position',txtBoxL(indr,:),...
    'HorizontalAlign','center',...
    'string','Remove',...
    'fontsize',8,'String',defaults{lm,1},...
    'callback',{@txt_call,lines(lm,1)},...
    'BackgroundColor','w');
S.outT(2,indr) = uicontrol('style','edit',...
    'units','pix',...
    'position',txtBoxR(indr,:),...
    'HorizontalAlign','center',...
    'string','Remove',...
    'fontsize',8,'String',defaults{lm,2},...
    'callback',{@txt_call,lines(lm,2)}...
    ,'BackgroundColor','w');
S.outD(1,indr) = uicontrol('style','text',...
    'units','pix',...
    'position',diaL(indr,:),...
    'HorizontalAlign','left',...
    'string','Remove',...
    'fontsize',8,'String',names{lm,1},'BackgroundColor',bckColor);
S.outD(2,indr) = uicontrol('style','text',...
    'units','pix',...
    'position',diaR(indr,:),...
    'HorizontalAlign','left',...
    'string','Remove',...
    'fontsize',8,'String',names{lm,2},'BackgroundColor',bckColor);
end
set(S.loader,'Callback',{@loader_func,S});

function [] = rmv_elmnt(elmnt,fName)
    %remove element from selected field. Does nothing if element is not in
    %selected field
    exstSt = get(S.(char(fName)),'String');
    indI = strcmp(exstSt,elmnt);
    exstSt(indI) = [];
    set(S.(char(fName)),'String',exstSt);
end
function [] = add_elmnt(elmnt,fName)
    %add element from selected field. Does nothing if element is already
    %in selected field
    exstSt = get(S.(char(fName)),'String');
    indI = strcmp(exstSt,elmnt);
    if ~any(indI)
        exstSt{length(exstSt)+1} = elmnt;
        exstSt = sortrows(exstSt);
        set(S.(char(fName)),'String',exstSt);
    end
end
 
function [] = pb_call(varargin)
    % Callback for edit.
    S = varargin{3};
    Ced = get(S.ed,{'string','value'});
    newOut = cell(Ced{1}(Ced{2}));
    oldOut = cell(get(S.ls,'string'));
    L1 = length(newOut);
    L2 = length(oldOut);
    Out = cell(L1+L2,1);
    for j = 1:L1
        Out{j} = newOut{j};
    end
    for j = 1:L2
        Out{j+L1} = oldOut{j};
    end
    if isempty(Out)
        Out{1} = '';
    end
    
    set(S.ls,'string',sortrows(Out),'Value',1)  % Now set the listbox string to the value in E.
    Ced{1}(Ced{2}) = [];
    set(S.ed,'string',sortrows(Ced{1}),'Value',[])
    % add the selection
    lineT = 3;
    replaceFile(S,lineT,get(S.ls,'string'));
    %set(S.loader,'Value',0)
    fn = fieldnames(S);
    if any(strcmp(fn,'loader'))
        set(S.loader,'Value',0)
    end
end
function [] = re_call(varargin)
    % Callback for edit.
    S = varargin{3};
    Ced = get(S.ls,{'string','value'});
    newOut = cell(Ced{1}(Ced{2}));
    oldOut = cell(get(S.ed,'string'));
    L1 = length(newOut);
    L2 = length(oldOut);
    Out = cell(L1+L2,1);
    for j = 1:L1
        Out{j} = newOut{j};
    end
    for j = 1:L2
        Out{j+L1} = oldOut{j};
    end
    if isempty(Out)
        Out{1} = '';
    end

    Ced{1}(Ced{2}) = [];
    if isempty(Ced{1})
        Ced{1} =outUn(strcmp(outA,defSlt));
        indI = strcmp(outA,defSlt);
        nme  = outUn(indI);
        indI = strcmp(Out,nme);
        Out(indI) = [];
    end
    set(S.ed,'string',sortrows(Out),'Value',1)   
    set(S.ls,'string',sortrows(Ced{1}),'Value',1)
    lineT = 3;
    replaceFile(S,lineT,get(S.ls,'string'));
    % check for unique!!!
    fn = fieldnames(S);
    if any(strcmp(fn,'loader'))
        set(S.loader,'Value',0)
    end
end
function [] =  txt_call(varargin)
    if length(varargin)>3
        blip = false;
    else 
        blip = true;
    end
    txt = get(varargin{1},'String');
    rplcIn = varargin{3};
    %left to right from the bottom up
    tempStr = get(S.file,'String');
    rplcSt = char(tempStr(rplcIn));
    stI = strfind(rplcSt, ' ');
    if isempty(stI)
        error('wrong text format in default file')
    end
    stI = min(stI);
    rplcSt = [txt rplcSt(stI:length(rplcSt))];
    
    dlI = strfind(rplcSt, '#')-1;
    numS= length(strfind(rplcSt(1:dlI(1)),' '));
    numC= dlI-numS;
    len = numC*spcFrac+numS;
    
    ndSpc = num2delim-len;
    if gt(ndSpc,0)
        rplcSt = [rplcSt(1:dlI) ' '*ones(1,floor(ndSpc)) ...
            rplcSt(dlI:length(rplcSt))];
    elseif lt(ndSpc,0)
        rplcSt = [rplcSt(1:dlI+floor(ndSpc)) rplcSt(dlI:length(rplcSt))];
    end
    tempStr{rplcIn} = rplcSt;
    set(S.file,'String',tempStr)
    if blip
        set(S.loader,'Value',0)
    end
end

function[] = replaceFile(S,lineT,Out)
    indi = 1:length(Out);
    if eq(lineT,3)
        txt = '';
        for i = 1:length(Out)
            indi(i) = strmatch(Out{i},outUn, 'exact');
            txt = [txt outA{indi(i)} ', '];
        end
        txt = txt(1:length(txt)-2);
    else
        txt = Out;
    end
    temTx = get(S.file,'String');
    tt = temTx{lineT};
    strt = strfind(tt,'#');
    tt = tt(strt:length(tt));
    temTx{lineT} = [txt ' '*ones(1,9) tt];
    set(S.file,'String',temTx);
end

function[] = loader_func(varargin)
    S = varargin{3};
    if get(S.loader,'value')
        [ outputs,outRs,maxZ,wndH,wndAv,lyrAv,outWn,wtrMx,wtrMn,...
            wndMx,wndMn,drhDz,Tdiff,plotYes,writeYes] = ...
        OpenCfg(LakeName,directory);
        exstOut = get(S.ls,'String');
        for lk = 1:length(exstOut)
            add_elmnt(exstOut{lk},'ed');
            rmv_elmnt(exstOut{lk},'ls');
        end
        for lk = 1:length(outputs)
            indI = strcmp(outputs{lk},outA);
            add_elmnt(outUn{indI},'ls');
            rmv_elmnt(outUn{indI},'ed');
        end
        nOuts = get(S.ls,'String');
        replaceFile(S,3,nOuts)
    if plotYes
        plotYes = 'Y';
    else
        plotYes = 'N';
    end
    if writeYes
        writeYes = 'Y';
    else
        writeYes = 'N';
    end
    set(S.outT(1,1),'String',plotYes)
    set(S.outT(2,1),'String',writeYes)
    set(S.outT(1,2),'String',num2str(drhDz))
    set(S.outT(2,2),'String',num2str(Tdiff))
    set(S.outT(1,3),'String',num2str(wndMx))
    set(S.outT(2,3),'String',num2str(wndMn))
    set(S.outT(1,4),'String',num2str(wtrMx))
    set(S.outT(2,4),'String',num2str(wtrMn))
    set(S.outT(1,5),'String',num2str(lyrAv))
    set(S.outT(2,5),'String',num2str(outWn))
    set(S.outT(1,6),'String',num2str(wndAv))
    set(S.outT(2,6),'String',num2str(wndH))
    set(S.outT(1,7),'String',num2str(outRs))
    set(S.outT(2,7),'String',num2str(maxZ))
    for m = 1:2
        for n = 1:7
            txt_call(S.outT(m,n),[],lines(7-n+1,m),true)
        end
    end
    end
end

function[] =  publish(varargin)
    S = varargin{3};
    tempStr = get(S.file,'String');
    fileName = [directory '/' LakeName '.lke'];
    outFile = fopen(fileName,'w');
    if lt(outFile,0)
        error('file directory invalid')
    end
    wrt = @(writer)fprintf(outFile,writer); % build a subfunction that writes 
                                    % the contents of the input "writer" 
                                    % to the file everytime wrt is called
   
    for ldn = 1:length(tempStr)
        lin = char(tempStr{ldn});
        ind = strfind(lin,'#');
        if isempty(ind)
            wrt([lin '\r\n']);
        else
            %replace spaces with tabs
            spcI = strfind(lin(1:ind-1),'  ');
            lin = [lin(1:min(spcI)-1) '\t\t' lin(ind:length(lin))];
            wrt([lin '\r\n']);
        end
    end
    fclose all;
    close all;
    assignin('caller','done',true);
    
end

        
end
