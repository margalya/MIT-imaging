function [ ] = PicMean( appData )
%calculate average picture outside of imaging software. Does not load all
%files to memory.

[nums, files] = readCurrentFolder(appData.analyze.readDir, appData.consts.cameras{appData.options.cameraType}.fileFormat);

files = files( nums~=1000 ); %remove data-1000.mat from the list, in case it exist
[atoms, back, dark] = appData.data.camera.fileReadFunction( [appData.analyze.readDir '\' files(1).name] ); %load first file

[~,name,~] = fileparts(files(1).name);
load([appData.analyze.readDir '\' name '.mat'], 'savedData'); %load first mat file
isAbs = 1; % Absorption
if isfield(savedData.consts.cameras{appData.options.cameraType}, 'photonPerADU')
   isAbs = 0; % Fluorescence
end
newatoms = zeros(size( atoms ));
newback = zeros(size( back ));
newdark = zeros(size( dark ));
progressbar(0);
for j = 1 : length(files)
    %     load( [appData.analyze.readDir '\' files(j).name] );
    [atomsTemp, backTemp, darkTemp] = appData.data.camera.fileReadFunction( [appData.analyze.readDir '\' files(j).name]);
    %     atomsTemp = atomsTemp .* ( ~(atomsTemp<0)); % set all pixelvalues<0 to 0
    %     backTemp = backTemp .* ( ~(backTemp<0)); % set all pixelvalues<0 to 0
    
    %add new file into the sum
    %     pic = pic + absorption;
    newatoms = newatoms + double(atomsTemp);
    newback = newback + double(backTemp);
    newdark = newdark + double(darkTemp);
    
    %create absorption and add to the sum
    %     absorption = log( (newback + 1)./ (newatoms + 1)  );
    %     newabsorption = newabsorption + absorption;
    
    progressbar(j/length(files));
end

% Plot the resulting absorption image:
% newatoms = newatoms - newdark;                           % subtract the dark background from the atom pic
% newatoms = newatoms .* ( ~(newatoms<0));                                                   % set all pixelvalues<0 to 0
% newback =  newback - newdark;                              % subtract the dark background from the background pic
% newback = newback .* ( ~(newback<0));                                                         % set all pixelvalues<0 to 0
% 
% absorption = log( (newback + 1)./ (newatoms + 1)  );
% figure; imagesc(absorption(400:700,560:900))

clear atoms back dark absoprtion
%divide by number of files (to average), and convert back to unit16
atoms = uint16(newatoms / length(files));
back = uint16(newback / length(files));
dark = uint16(newdark / length(files));
% absopriton = newabsorption / length(matFiles); %#ok<NASGU>

savedData.data.date = datestr(now);
savedData.save.picNo = 1000;

save( [appData.analyze.readDir '\data-' num2str(1000) '.mat'] ,'savedData', 'atoms', 'back', 'dark') %, 'absopriton' can also be saved, if needed for averaging over absorption images and not on 'atoms' and 'back' images
end

function [nums, files] = readCurrentFolder(readDir, fileFormat)

% readDir = appData.analyze.readDir;
% fileFormat = appData.consts.cameras{appData.options.cameraType}.fileFormat;
files = dir([readDir '\*.' fileFormat]);
matFiles = dir([readDir '\*.mat']);
if isempty(files)
    nums = [];
    files = [];
    return
end
if any(strcmp({'data-', 'Data-'}, files(1).name(1:5) )) %new file name format
    nums = zeros(1, length(files));
    for j = 1 : length(files)
        nums(j) = str2double(files(j).name( 6 : end-4 ));
    end
    [nums, sortIndx] = sort(nums); %sort according to file number
    files = files(sortIndx);
else %old file name format
    files = dir([readDir '\*.' fileFormat]);
    [~, sortIndx] = sort([files.datenum], 'descend');
    files = files(sortIndx); %sort files according to modification date
    nums = zeros(1, length(files));
    for j = 1 : length(files)
        nums(j) = str2double(files(j).name( 12 : end-4 )); % file number used is only the time (HHMMSS), date omitted
    end
end

if any(arrayfun(@(x) strcmp(x.name,'data-1000.mat'),matFiles)) %if 1000.mat exist, for images saved inside mat files
    matFiles = matFiles( arrayfun(@(x) strcmp(x.name,'data-1000.mat'),matFiles) ); %keep only data-1000.mat in the list
    nums = [nums 1000];
    files = [files; matFiles];
end

end