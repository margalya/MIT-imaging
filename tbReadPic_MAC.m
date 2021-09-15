function savedData = tbReadPic_MAC(appData)

comment = [];
fileName = dir([ appData.analyze.readDir '/data-' num2str(appData.analyze.showPicNo) '-*.mat']);
if ( size(fileName, 1) > 1 )
    fileName = strtrim(fileName(1, :));
end
if ( isempty(fileName) )
    fileName = dir([appData.analyze.readDir '/data-' num2str(appData.analyze.showPicNo) '.mat']);
    if ( isempty(fileName) )
        warndlg({'File doesnt exist:', [appData.analyze.readDir '/data-' num2str(appData.analyze.showPicNo) '.mat']}, 'Warning', 'modal');
        return
    end
else
    dotIndex = find(fileName.name == '.');
    dashIndex = find(fileName.name == '-');
    %         comment = fileName(dashIndex(end):dotIndex-1);
    comment = fileName.name(dashIndex(2):dotIndex(end)-1);
end

load([appData.analyze.readDir '/' fileName.name], 'savedData');
if ( ~isempty(comment) )
    savedData.save.commentStr = comment;
end
if (length(appData.data.fits) > length(savedData.data.fits) )
    savedData.consts.fitTypes = appData.consts.fitTypes;
    for j=length(savedData.data.fits)+1 : length(appData.data.fits)
        savedData.data.fits{j} = appData.consts.fitTypes.fits{j};
    end
end
if (length(appData.consts.availableAnalyzing.str) > length(savedData.consts.availableAnalyzing.str) )
    savedData.consts.availableAnalyzing = appData.consts.availableAnalyzing;
    %         for i=length(savedData.data.fits)+1 : length(appData.data.fits)
    %             savedData.data.fits{i} = appData.data.fits{i};
    %         end
end
% appData = createAppData(savedData, appData);
% appData = analyzeAndPlot(appData);
    
end