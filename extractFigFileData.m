function [ x, y, U, L] = extractFigFileData( fileName )
% Extract figure y data from .fig file, without opening the file


tempFig = load( fileName , '-mat'); %load figure as mat file
close(gcf);
try
    x = tempFig.hgS_070000.children.children(1,1).properties.XData; %exctract data
    y = tempFig.hgS_070000.children.children(1,1).properties.YData; %exctract data
    if sum(strcmp(fieldnames(tempFig.hgS_070000.children.children(1,1).properties), 'UData')) == 1 %check for field exsitence
        U = tempFig.hgS_070000.children.children(1,1).properties.UData; %upper errorbar, usually the same as lower one
    end
    if sum(strcmp(fieldnames(tempFig.hgS_070000.children.children(1,1).properties), 'LData')) == 1
        L = tempFig.hgS_070000.children.children(1,1).properties.LData; %lower errorbar
    end
catch
    try
        x = tempFig.hgS_070000.children(1,1).children(1,1).properties.XData;
        y = tempFig.hgS_070000.children(1,1).children(1,1).properties.YData;
        if sum(strcmp(fieldnames(tempFig.hgS_070000.children(1,1).children(1,1).properties), 'UData')) == 1
            U = tempFig.hgS_070000.children(1,1).children(1,1).properties.UData;
        end
        if sum(strcmp(fieldnames(tempFig.hgS_070000.children(1,1).children(1,1).properties), 'LData')) == 1
            L = tempFig.hgS_070000.children(1,1).children(1,1).properties.LData;
        end
    catch
        x = tempFig.hgS_070000.children(2).children(1).properties.XData;
        y = tempFig.hgS_070000.children(2).children(1).properties.YData;
        if isfield(tempFig.hgS_070000.children(2).children(1).properties, 'UData')
            U = tempFig.hgS_070000.children(2).children(1).properties.UData;
        end
        if isfield(tempFig.hgS_070000.children(2).children(1).properties, 'LData')
            L = tempFig.hgS_070000.children(2).children(1).properties.LData;
        end
    end
end

end

