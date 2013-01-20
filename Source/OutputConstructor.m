function [truthTable,outputOptions,writeTable,plotTable,...
    dateInput,dateOutput,delimI,delimO] ...
    = OutputConstructor(outputNames,pltMods)

outputOptions = {'St','uSt','Ln','W','wTemp','wndSpd','metaT',...
     'metaB','thermD','SthermD','SmetaB','SmetaT','SuSt','SLn','SW',...
     'N2','SN2','T1','ST1'};
                        % writable outputs
programOptions= {'openWtr','openWnd','openBth','errCkWtr','errCkWnd',...
    'findLyr','ssnLyr','uStYes','SuStYes','StYes','dwnSmple'};
                        % program flow specifiers

%% Error Check
if ne(length(outputNames),sum(ismember(outputNames,outputOptions)))
    error(['output "' ...
        outputNames{~ismember(outputNames,outputOptions)} ...
        '" not recognized'])
end
%% Other defaults
dateInput  = 'yyyy-mm-dd HH:MM';
delimI = '\t';

dateOutput  = 'yyyy-mm-dd HH:MM';
delimO = '\t';

%% Figure defaults
isStringMod = {'figUnits','figType','fontName','figRes'};
figUnits    = 'inches';
figWidth    = 6;   % relative to fig_units
figHeight   = 3;   % relative to fig_units
leftMargin  = .75; % relative to fig_units
rightMargin = .1;  % relative to fig_units
topMargin   = .4;  % relative to fig_units
botMargin   = .4;  % relative to fig_units
figType     = 'png';
figRes      = '150'; % dots per inch (not relative to units?)
fontName    = 'Arial';
fontSize    = 12;
heatMapMin  = 0;
heatMapMax  = 30;


plt = struct('figUnits',figUnits,'figWidth',figWidth,'figHeight',figHeight,...
    'leftMargin',leftMargin,'rightMargin',rightMargin,'topMargin',topMargin,...
    'botMargin',botMargin,'figType',figType,'figRes',figRes,...
    'fontName',fontName,'fontSize',fontSize,'heatMapMin',heatMapMin,...
    'heatMapMax',heatMapMax);

if ~isempty(pltMods)
    % use plot mods to modify plotting defaults
    fN = fieldnames(pltMods);
    for n = 1:length(fN)
        if ~any(strcmp(fN{n},isStringMod))
            if ne(str2double(pltMods.(fN{n})),plt.(fN{n}))
                fprintf(['>>User plot modification; replacing ' fN{n} '=' ...
                    num2str(plt.(fN{n})) ' with ' pltMods.(fN{n})  '<<<\n']);
                plt.(fN{n}) = str2double(pltMods.(fN{n}));
            end
        else
            if ~strcmp(pltMods.(fN{n}),plt.(fN{n}))
                fprintf(['>>User plot modification; replacing ' fN{n} '=' ...
                    plt.(fN{n}) ' with ' pltMods.(fN{n})  '<<<\n']);
                plt.(fN{n}) = pltMods.(fN{n});
            end
        end
    end
end


fig_Defaults = struct('Units',plt.figUnits,'Color','w',...
    'PaperUnits',plt.figUnits,...
    'PaperPosition',[0 0 plt.figWidth plt.figHeight],...
    'Position',[1 1 plt.figWidth plt.figHeight],...
    'PaperPositionMode','manual',...
    'PaperSize',[plt.figWidth plt.figHeight]);

print_Defaults = struct('format',['-d' plt.figType],'res',['-r' plt.figRes],...
    'toClose',true);

position = [plt.leftMargin/plt.figWidth plt.botMargin/plt.figHeight ...
    (plt.figWidth-plt.leftMargin-plt.rightMargin)/plt.figWidth ...
    (plt.figHeight-plt.topMargin-plt.botMargin)/plt.figHeight];

axes_Defaults = struct('FontName',plt.fontName,'FontSize',plt.fontSize,...
    'Layer','top','Position',position,...
    'Box','on',...
    'YLabel','','Title','','YDir','normal','YScale','linear','CLim',[0 1]);
%% build Structures
plotTable = struct('FigD',fig_Defaults,'PrintD',print_Defaults);
for j = 1:length(programOptions)
    truthTable.(char(programOptions{j})) = false;
end
for j = 1:length(outputOptions)
    truthTable.(['wrt_' char(outputOptions{j})]) = false;
    writeTable.(char(outputOptions{j})) = false;
    plotTable.(char(outputOptions{j})) = axes_Defaults;
    
end
for j = 1:length(outputOptions)
    outputConstruct.(char(outputOptions{j})) = truthTable;
end

%% Schmidt Stability
name = 'St';
RunNeed.(char(name)) = {'openWtr','openBth','StYes','wrt_St'};
RunAxes.(char(name)) = struct('YLabel','Schmidt Stability (J m^{-2})',...
    'Title','Schmidt Stability');

%% U-Star
name = 'uSt';
RunNeed.(char(name)) = {'openWtr','openWnd','openBth','findLyr','uStYes',...
    'wrt_uSt'};
RunAxes.(char(name)) = struct('Ylabel','U* (m s^{-1})','Title','Ustar');

%% Lake Number
name = 'Ln';
RunNeed.(char(name)) = {'openWtr','openWnd','openBth','findLyr','uStYes',...
    'StYes','wrt_Ln'};
RunAxes.(char(name)) = struct('YLabel','Lake Number (dim)',...
    'YScale','log','Title','Lake Number');

%% Wedderburn Number
name = 'W';
RunNeed.(char(name)) = {'openWtr','openWnd','openBth','findLyr','uStYes',...
    'wrt_W'};
RunAxes.(char(name)) = struct('YLabel','Wedderburn Number (dim)',...
    'YScale','log','Title','Wedderburn Number');

%% Water Temperature
name = 'wTemp';
RunNeed.(char(name)) = {'openWtr','wrt_wTemp'};
RunAxes.(char(name)) = struct('YLabel','Depth (m)','YDir','reverse',...
    'Title','Water Temperature','Clim',[plt.heatMapMin plt.heatMapMax]);

%% Wind Speed
name = 'wndSpd';
RunNeed.(char(name)) = {'openWnd','wrt_wndSpd'};
RunAxes.(char(name)) = struct('YLabel','Wind Speed (m s^{-1})',...
    'Title','Wind Speed');

%% Metalimnion Top
name = 'metaT';
RunNeed.(char(name)) = {'openWtr','findLyr','wrt_metaT'};
RunAxes.(char(name)) = struct('YLabel','Metalimnion Top (m)',...
    'YDir','reverse','Title','Metalimnion Top');

%% Metalimnion Bottom
name = 'metaB';
RunNeed.(char(name)) = {'openWtr','findLyr','wrt_metaB'};
RunAxes.(char(name)) = struct('YLabel','Metalimnion Bottom (m)',...
    'YDir','reverse','Title','Metalimnion Bottom');

%% Thermocline Depth
name = 'thermD';
RunNeed.(char(name)) = {'openWtr','findLyr','wrt_thermD'};
RunAxes.(char(name)) = struct('YLabel','Thermocline (m)',...
    'YDir','reverse','Title','Thermocline');

%% Parent Thermocline
name = 'SthermD';
RunNeed.(char(name)) = {'openWtr','findLyr','ssnLyr','wrt_SthermD'};
RunAxes.(char(name)) = struct('YLabel','Parent Thermocline (m)',...
    'YDir','reverse','Title','Parent Thermocline');

%% Parent Metalimnion Bottom
name = 'SmetaB';
RunNeed.(char(name)) = {'openWtr','findLyr','ssnLyr','wrt_SmetaB'};
RunAxes.(char(name)) = struct('YLabel','Parent Metalimnion Bottom (m)',...
    'YDir','reverse','Title','Parent Metalimnion Bottom');

%% Parent Metalimnion Top
name = 'SmetaT';
RunNeed.(char(name)) = {'openWtr','findLyr','ssnLyr','wrt_SmetaT'};
RunAxes.(char(name)) = struct('YLabel','Parent Metalimnion Top (m)',...
    'YDir','reverse','Title','Parent Metalimnion Top');

%% Parent U-Star
name = 'SuSt';
RunNeed.(char(name)) = {'openWtr','openWnd','openBth','findLyr','ssnLyr',...
    'SuStYes','wrt_SuSt'};
RunAxes.(char(name)) = struct('YLabel','Seasonal U* (m s^{-1})',...
    'Title','Seasonal Ustar');

%% Parent Lake Number
name = 'SLn';
RunNeed.(char(name)) = {'openWtr','openWnd','openBth','findLyr','ssnLyr',...
    'StYes','SuStYes','wrt_SLn'};
RunAxes.(char(name)) = struct('YLabel','Seasonal Lake Number (dim)',...
    'Title','Seasonal Lake Number','YScale','log');

%% Parent Wedderburn Number
name = 'SW';
RunNeed.(char(name)) = {'openWtr','openWnd','openBth','findLyr','ssnLyr',...
    'SuStYes','wrt_SW'};
RunAxes.(char(name)) = struct('YLabel','Seasonal Wedderburn Number (dim)',...
    'Title','Seasonal Wedderburn Number','YScale','log');

%% Buoyancy Frequency
name = 'N2';
RunNeed.(char(name)) = {'openWtr','findLyr','wrt_N2'};
RunAxes.(char(name)) = struct('YLabel','Buoyancy Frequency (s^{-2})',...
    'Title','Buoyancy Frequency');

%% Parent Buoyancy Frequency
name = 'SN2';
RunNeed.(char(name)) = {'openWtr','findLyr','ssnLyr','wrt_SN2'};
RunAxes.(char(name)) = struct('YLabel','Parent Buoyancy Frequency (s^{-2})',...
    'Title','Parent Buoyancy Frequency');

%% Mode 1 vertical seiche
name = 'T1';
RunNeed.(char(name)) = {'openWtr','openBth','findLyr','wrt_T1'};
RunAxes.(char(name)) = struct('YLabel','Mode 1 Period (s)',...
    'Title','Mode 1 Vertical Seiche');

%% Parent Mode 1 vertical seiche
name = 'ST1';
RunNeed.(char(name)) = {'openWtr','openBth','findLyr','ssnLyr','wrt_ST1'};
RunAxes.(char(name)) = struct('YLabel','Mode 1 Period (s)',...
    'Title','Parent Mode 1 Vertical Seiche');

%% Fill Construct
fNms = fieldnames(outputConstruct);
for k = 1:length(fNms)
    strN = RunNeed.(char(fNms{k}));
    for j = 1:length(strN)
        outputConstruct.(char(fNms{k})).(char(strN{j})) = true;
    end
    axN = fieldnames(RunAxes.(char(fNms{k})));
    for j = 1:length(axN)
        plotTable.(char(fNms{k})).(char(axN{j})) = ...
            RunAxes.(char(fNms{k})).(char(axN{j}));
    end
    
end

fNms = fieldnames(truthTable);
for k = 1:length(fNms)
    for j = 1:length(outputNames)
        if outputConstruct.(char(outputNames{j})).(char(fNms{k}))
            truthTable.(char(fNms{k})) = true;
        end
    end
end

fNms = fieldnames(writeTable);
for k = 1:length(fNms)
    for j = 1:length(outputNames)
        if outputConstruct.(char(outputNames{j})).(['wrt_' char(fNms{k})])
            writeTable.(char(fNms{k})) = true;
        end
    end
end




end

