function [dates,dat,headers] = gFileOpen(fileName)
% author: Luke A Winslow March 2011
    fid = fopen(fileName);

    headers = fgetl(fid);
    format = '%s';
    for i=2:length(regexp(headers,'\t','split'));
        format = strcat(format,'%f');
    end
    d = textscan(fid,format,'delimiter','\t','treatAsEmpty',{'na','NA',...
        '#VALUE!','#NAME?'});
    fclose(fid);
    
    
    try dates = datenum(d{1},'yyyy-mm-dd HH:MM');
        
    catch mssg
        
        dates = NaN(length(d{1}),1);
        for i = 1:length(d{1})
            dates(i) = datenum(char(d{1}(i)),'yyyy-mm-dd HH:MM');
        end
    end
    dat = NaN(length(dates),length(d)-1);
    for j = 1:length(dates)
        dat(l,1:length(d{j+1})) = d{j+1};
    end
    
end