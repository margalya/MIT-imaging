
function plotImage(appData)

lineColor = 0.85;
colors = ['r', 'k', 'g', 'c'];

% setWinName(appData);

%
% plot the main image
%
% [pic, xFit, yFit] = appData.data.fits{appData.data.fitType}.getPicData(appData);
% [xData, yData] = appData.data.plots{appData.data.plotType}.getXYData(appData);
% [h w] = size(pic);
% x = [1 : w];
% y = [1 : h];
[pic x0 y0] = appData.data.plots{appData.data.plotType}.getPic();
if isempty(pic)
    cla(appData.ui.plot);
    cla(appData.ui.xPlot);
    cla(appData.ui.yPlot);
    return;
elseif numel(pic)==1 %single pixel case (e.g. super-pixel of fluoresence image)
    set(0, 'CurrentFigure', appData.ui.win)
    set(appData.ui.win,'CurrentAxes',appData.ui.plot);
    image( pic);
    appData.data.fits{appData.data.fitType}.plotFitResults(appData);
else
    pic = appData.data.plots{appData.data.plotType}.normalizePic(appData, pic);
%     pic = LowPassFilter(pic, 3, 1);
    % assignin('base','pic',pic);
    % pic = appData.data.fits{appData.data.fitType}.normalizePic(pic);
    [h w] = size(pic);
    x = [1 : w];
    y = [1 : h];
    chipStart = appData.data.camera.chipStart;
    % [xCenter yCenter] = appData.data.ROITypes{appData.data.ROIUnits}.getCenter(appData, appData.data.fits{appData.data.fitType});
    %            x0 = fitObj.xCenter; %TODO: round, min max
    %            y0 = fitObj.yCenter;
    xCenter = appData.data.fits{appData.data.fitType}.xCenter;
    yCenter = appData.data.fits{appData.data.fitType}.yCenter;
    xCenter = max(xCenter, x0);
    yCenter = max(yCenter, y0);
    xCenter = min(xCenter, x0+w-1-1);
    yCenter = min(yCenter, y0+h-1-1);
    xSz = round(appData.data.fits{appData.data.fitType}.xUnitSize)*2;
    ySz = round(appData.data.fits{appData.data.fitType}.yUnitSize)*2;
    if appData.data.fitType == appData.consts.fitTypes.onlyMaximum
        xSz = min(xCenter-x0-10, w-(xCenter-x0)-10);%round(w/2-10);
        ySz = min(yCenter-y0-10, h-(yCenter-y0)-10);%round(h/2-10);
    else
        xSz = min(xSz, round(w/2-10));
        ySz = min(ySz, round(h/2-10));
    end
    % xSz = max(xCenter-x0+1-xSz, 1);
    % ySz = max(yCenter-y0+1-ySz, 1);
    % xSz = min(xCenter-x0+1+xSz, w);
    % ySz = min(yCenter-y0+1+ySz, h);
    [xFit yFit] = appData.data.fits{appData.data.fitType}.getXYFitVectors(x+x0-1, y+y0-1);
    [xData yData] = appData.data.plots{appData.data.plotType}.getXYDataVectors(xCenter, yCenter, appData.options.avgWidth);
    % assignin('base','yData',yData);
    % xlswrite('E:\Dropbox\MATLAB\Fringes_raw_data.xlsx', yData', 1, [char(appData.save.picNo-985+'A'-1) '1']) %export yData into excel for Yoni
    % pic( [yCenter : yCenter+1] -y0+1, :) = ones(2, w) * lineColor;  % a line at the center (x axis)
    % pic( [yCenter : yCenter+1] -y0+1, [1:xCenter-x0+1-xSz xCenter-x0+1+xSz:w]) = ones(2, w-2*xSz+1) * lineColor;
    % pic(:, [xCenter : xCenter+1] -x0+1 ) = ones(h, 2) * lineColor; % a line at the center (y axis)
    % pic([1:yCenter-y0+1-ySz yCenter-y0+1+ySz:h], [xCenter : xCenter+1] -x0+1 ) = ones(h-2*ySz+1, 2) * lineColor;
    if all(size(pic)>200) %don't plot ROI square for small images
        pic = appData.data.plots{appData.data.plotType}.createSquare(appData, pic);
    end
    
    appData.data.imageWidth = round(appData.consts.maxplotSize*(w / h));
    appData.data.imageHeight = appData.consts.maxplotSize;
    if ( appData.data.imageWidth > appData.consts.maxplotSize )
        appData.data.imageWidth = appData.consts.maxplotSize;
        appData.data.imageHeight = round(appData.consts.maxplotSize*( h / w));
    end
    set(appData.ui.plot, 'Position', [5 5 appData.data.imageWidth appData.data.imageHeight]);
    set(appData.ui.xPlot, 'Position', [5 5+appData.data.imageHeight+appData.consts.strHeight appData.data.imageWidth appData.consts.xyPlotsHeight]);
    set(appData.ui.yPlot, 'Position', [5+appData.data.imageWidth+appData.consts.strHeight*1.5 5 appData.consts.xyPlotsHeight appData.data.imageHeight]);
    
%     figure(appData.ui.win);
    set(0, 'CurrentFigure', appData.ui.win)
    set(appData.ui.win,'CurrentAxes',appData.ui.plot);
    colormap(jet(256));
    % image( ([x(1) x(end)]+x0-1)*appData.data.camera.xPixSz * 1000, ...
    %     ([y(1) y(end)]+y0-1-chipStart-1)*appData.data.camera.yPixSz * 1000, pic*256);
    
    plotUnits = get(appData.ui.pmPlotUnits, 'Value');
    switch plotUnits
        case 1 % plot using mm units
            image( ([x(1) x(end)]+x0-1)*appData.data.camera.xPixSz * 1000, ...
                ([y(1) y(end)]+y0-1-chipStart-1)*appData.data.camera.yPixSz * 1000, pic*256);
%             hold on;
%             N = appData.data.fits{6}.atomsNo;
%             EF = 1.05e-34*2*pi*(600*25e3*25e3*6*N)^(1/3);
%             TOF = 0.5e-3;
%             RF = sqrt(2*EF/9.9883414e-27)*TOF*1e3; %Fermi Radius in mm
%             viscircles([ mean([x(1) x(end)]+x0-1)*appData.data.camera.xPixSz*1000 mean([y(1) y(end)]+y0-1-chipStart-1)*appData.data.camera.yPixSz*1000], RF,...
%                 'Color', 'White', 'LineWidth', 0.5);
            
%             plot(4.5859,3.3475,'+white', 'markersize', 34, 'linewidth', 1)
%             plot(3.9926,2.909,'+white', 'markersize', 34, 'linewidth', 1)
        case 2 % plot using pixel units
            image(pic*256);
%             hold on;
%             viscircles([yCenter xCenter],0.4359*50)
            
%             hold on;
%             plot(635,450,'+white', 'markersize', 34, 'linewidth', 1.5)
%             text(890,750,' X marks the spot', 'Color', 'white')
    end
    set( appData.ui.plot, 'XAxisLocation', 'top');
    set( appData.ui.plot, 'YAxisLocation', 'right');
    
    %
    % plot x plot
    %
    set(0, 'CurrentFigure', appData.ui.win)
    % figure(appData.ui.win);
    set(appData.ui.win,'CurrentAxes',appData.ui.xPlot);
    plot((x+x0-1)*appData.data.camera.xPixSz * 1000, xData, 'b');
    hold on
    % create and plot the xFit vector
    if ( length(xFit) == length(x) )
        for i = 1 : size(xFit, 1)
            plot((x+x0-1)*appData.data.camera.xPixSz * 1000, xFit(i, :), colors(i));
        end
    end
    hold off
    set( appData.ui.xPlot, 'XAxisLocation', 'top');
    set( appData.ui.xPlot, 'YAxisLocation', 'right');
    switch plotUnits
        case 1 % plot using mm units
            xlabel(appData.ui.xPlot, 'Distance [mm]', 'FontSize', appData.consts.fontSize);
        case 2 % plot using pixel units
            xlabel(appData.ui.xPlot, 'Distance [pixels]', 'FontSize', appData.consts.fontSize);
    end
    ylabel(appData.ui.xPlot, 'Optical Density', 'FontSize', appData.consts.fontSize);
    xtick = get(appData.ui.plot, 'XTick');
    set(appData.ui.xPlot, 'XTick', xtick);%/(appData.data.camera.xPixSz * 1000)-xFit(1));%(ytick/(appData.data.camera.yPixSz* 1000)-yFit(1)+appData.data.camera.chipStart)
    set(appData.ui.xPlot, 'XTickLabel', []);
    if length(x)>1 % avoid single pixel plotting problem
        set(appData.ui.xPlot, 'XLim', ([x(1) x(end)]+x0-1)*appData.data.camera.xPixSz * 1000);%[1 width]);
    end
    % ylim = get(appData.ui.xPlot, 'YLim');
    if min([xData yData]) ~= max([xData yData]) % avoid single pixel plotting problem
        set(appData.ui.xPlot, 'YLim', [min([xData yData]) max([xData yData])]);%[min([yData yFit 0]) ylim(2)]);
    end
    
    %
    % plot y plot
    %
    set(0, 'CurrentFigure', appData.ui.win)
    % figure(appData.ui.win);
    set(appData.ui.win,'CurrentAxes',appData.ui.yPlot);
    % plot(xData, fliplr(y), 'b');
    plot(yData, fliplr(y-y0+1+chipStart+1)*appData.data.camera.yPixSz * 1000, 'b');
    hold on
    % create and plot they Fit vector
    if ( length(yFit) == length(y) )
        %     plot(xFit, fliplr(y), 'r');
        for i = 1 : size(yFit, 1)
            %         plot((x+x0-1)*appData.data.camera.xPixSz * 1000, xFit(i, :), colors(i));
            plot(yFit(i, :), fliplr(y-y0+1+chipStart+1)*appData.data.camera.yPixSz * 1000,  colors(i));
        end
    end
    hold off
    set( appData.ui.yPlot, 'XAxisLocation', 'top');
    set( appData.ui.yPlot, 'YAxisLocation', 'right');
    switch plotUnits
        case 1 % plot using mm units
            ylabel(appData.ui.yPlot, 'Distance [mm]', 'FontSize', appData.consts.fontSize);
        case 2 % plot using pixel units
            ylabel(appData.ui.yPlot, 'Distance [pixels]', 'FontSize', appData.consts.fontSize);
    end
    
    xlabel(appData.ui.yPlot, 'Optical Density', 'FontSize', appData.consts.fontSize);
    ytick = get(appData.ui.plot, 'YTick');
    % set(appData.ui.yPlot, 'YTick', fliplr(ytick));%fliplr(length(yFit)-(ytick/(appData.data.camera.yPixSz* 1000)-yFit(1)+appData.data.camera.chipStart)));
    set(appData.ui.yPlot, 'YTickLabel', []);
    if length(y)>1 % avoid single pixel plotting problem
        set(appData.ui.yPlot, 'YLim', ([y(1) y(end)]-y0+1+chipStart+1)*appData.data.camera.yPixSz * 1000);%[1 length(yFit)]);%yFit(end)-yFit(1)]);
    end 
    % xlim = get(appData.ui.yPlot, 'XLim');
    if min([xData yData]) ~= max([xData yData]) % avoid single pixel plotting problem
        set(appData.ui.yPlot, 'XLim', [min([xData yData]) max([xData yData])]);
    end
    set(appData.ui.yPlot, 'YTick', y(end)*appData.data.camera.yPixSz * 1000 - fliplr(ytick));
    
    
    %
    % plot results
    %
    set(0, 'CurrentFigure', appData.ui.win)
    % figure(appData.ui.win);
    set(appData.ui.win,'CurrentAxes',appData.ui.tmp);
    cla(appData.ui.tmp, 'reset');
    set(appData.ui.tmp, 'XLim', [0 1]);
    set(appData.ui.tmp, 'YLim', [0 1]);
    appData.data.fits{appData.data.fitType}.plotFitResults(appData);
    % linkaxes([appData.ui.plot, appData.ui.yPlot], 'y')
    % linkaxes([appData.ui.plot, appData.ui.xPlot], 'x')
    dcm_obj = datacursormode(appData.ui.win);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn,pic})
end
end

function txt = myupdatefcn(~,event_obj,t)
% Customizes text of data tips
pos = get(event_obj,'Position');
I = get(event_obj, 'DataIndex');
txt = {['X: ',num2str(pos(1))],...
    ['Y: ',num2str(pos(2))],...
    ['I: ',num2str(I)],...
    ['T: ',num2str(t(I))]};
end

% print(appData.ui.win, 'D:\My Documents\Documents\PhD and Fellowships\PhD Research Proposal\pics\MOT_analysis.eps');
% print('-depsc2', '-opengl', ['-f' num2str(appData.ui.win)], '-r864', 'D:\My Documents\Documents\PhD and Fellowships\PhD Research Proposal\pics\MOT_analysis.eps');

% assignin('base','xPlot',appData.ui.xPlot)

% assignin('base','mainPlot',appData.ui.plot)
% assignin('base','yPlot',appData.ui.yPlot)

% run('E:\Dropbox\Presentations\FRISNO 2017\T1_4us_average_eps.m')
% run('E:\Dropbox\MATLAB\imaging working copy\T1_4us_average.m')
