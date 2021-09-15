function [ ] = setSaveParam( folder, values, varargin)
% goes over all .mat files in the given folder, and changes the saveParamVal according to a given vector of values

% matFiles = dir([folder '\*.tif']); %sort by tif! and not by mat files
[nums, matFiles] = readCurrentFolder(folder);
% [~, sortIndx] = sort([matFiles.datenum], 'descend');
% matFiles = matFiles(sortIndx);
    
% [matFiles, nums] = matFilesListInternal(folder);
if isempty(matFiles)
    disp('No mat files found');
    return
end
if ~isempty(varargin)
    startIndex = varargin{1};
    startIndex = find(nums == startIndex); %change to picture number
    if isempty(startIndex)
        disp('Cannot find start index file, aborting')
        return
    end
    if numel(varargin)>1 %user has specified an end index
        endIndex = varargin{2};
        endIndex = find(nums == endIndex, 1); %change to picture number
        if isempty(endIndex)
            disp('Cannot find end index file, aborting')
            return
        end
    else
        endIndex = length(nums);
    end
else
    startIndex = 1;
    endIndex = length(nums);
end

% if length(startIndex : length(nums)) ~= length(values)
if length(startIndex : endIndex) ~= length(values)
    disp('Different number of files and values, aborting')
    return
end

fprintf(1,'[PicNum, New saveParamVal] = \n');
progressbar(0);
for j = startIndex : length(nums)
    [~, fileNameOnly] = fileparts(matFiles(j).name);
    load( [ folder '\'  fileNameOnly '.mat'], 'savedData' );
    savedData.save.saveParamVal = values(j-startIndex+1);
%     savedData.options.detuning = 10;
    fprintf(1,'[%u, ', nums(j));
    fprintf(1,'%u, ]\n', savedData.save.saveParamVal);
    save([ folder '\'  fileNameOnly '.mat'], 'savedData');
    clear savedData
    progressbar(j/length(nums));
end

% fprintf(1,'\n');

end

function [nums, files] = readCurrentFolder(readDir)

files = dir([readDir '\*.tif']);
if isempty(files)
    files = dir([readDir '\*.FITS']);
end
nums = zeros(1, length(files));
for j = 1 : length(files)
    nums(j) = str2double(files(j).name( 6 : end-4 ));
end
[nums, sortIndx] = sort(nums); %sort according to file number
files = files(sortIndx);

end

% function [ matFiles, nums ] = matFilesListInternal( folder )
% % return a sorted list of mat files in 'folder', while excluding pic 1000,1001.
% % list is sorted according to picNo
% matFiles = dir( [folder '\*.mat'] );
% % matFiles = matFiles( ~arrayfun(@(x) strcmp(x.name,'data-1000.mat'),matFiles) ); %remove data-1000.mat from the list
% 
% nums = zeros(1, length(matFiles));
% for j = 1 : length(matFiles)
%     dotIndex = find(matFiles(j).name == '.');
%     dashIndex = find(matFiles(j).name == '-');
%     if ( length(dashIndex) == 1 )
%         nums(j) = str2double(matFiles(j).name(dashIndex(1)+1 : dotIndex(end)-1));
%     else
%         nums(j) = str2double(matFiles(j).name(dashIndex(1)+1 : dashIndex(2)-1));
%     end
% end
% 
% [nums, sortIndex] = sort(nums);
% matFiles = matFiles(sortIndex);
% matFiles = matFiles(nums<1000);
% nums = nums(nums<1000);
% end

