function [ pltMods ] = pltFileOpen( fileName ) 

% author: Jordan S Read Jan 2013
format = '%s %s';
fid = fopen(fileName);
d = textscan(fid,format,'delimiter','\t');
fclose(fid);

if ne(length(d{1}),length(d{2}))
    error('plot parameters must be given as tab delimited pairs')
end

pltMods = struct('name',[]);

for fN = 1:length(d{1})
    pltMods.(char(d{1}(fN))) = char(d{2}(fN));
end


pltMods = rmfield(pltMods,'name');
end

