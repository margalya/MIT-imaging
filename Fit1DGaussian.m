classdef Fit1DGaussian < FitTypes
%FIT1DGAUSSIAN Summary of this class goes here
%   Detailed explanation goes here

   properties ( Constant = true )
       ID = 'Fit1DGaussian';
   end
   properties (SetAccess = private )
       xRes = [];
       xGof = [];
       xOutput = [];
       yRes = [];
       yGof = [];
       yOutput = [];
%        ODx = -1;
%        ODy = -1;
%        x0 = -1;
%        y0 = -1;
%        sigmaX = -1;
%        sigmaY = -1;
%        Cx = -1;
%        Cy = -1;
   end

   methods 
       function appData = analyze(obj, appData) % do the analysis
           [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
           [h w] = size(pic);           %binning and max
           %            binnedData = binning ( pic, appData.options.avgWidth*2+1);
           binnedData = LowPassFilter(pic, 30, appData.options.avgWidth);
           [maxes, indexes] = max(binnedData);                     % find maximum
           [maxValue, xPosMax] = max(maxes);
           yPosMax = indexes(xPosMax);   
%            obj.maxVal = maxValue;% / (appData.options.avgWidth*2)^2; %no need, already done in binning.m
           % center
           xCenter = xPosMax + x0; %(appData.options.avgWidth*2+1) * (xPosMax ) + x0-appData.options.avgWidth;
           yCenter = yPosMax + y0; %(appData.options.avgWidth*2+1) * (yPosMax ) + y0-appData.options.avgWidth; 
%            obj.xCenter = appData.options.avgWidth*2 * (xPosMax - 0.5) + x0;
%            obj.yCenter = appData.options.avgWidth*2 * (yPosMax - 0.5) + y0; 
%            % unit size
%            obj.xUnitSize = w/2;
%            obj.yUnitSize = h/2;   

            % fitting
%             [xData yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, xCenter, yCenter, appData.options.avgWidth);
% 
%           
%             [obj.xRes, obj.xGof, obj.xOutput] = fit(x', xData', 'gauss1'); 
%             [obj.yRes, obj.yGof, obj.yOutput] = fit(y', yData', 'gauss1');  
%             xCenter = round(obj.xRes.b1);
%             yCenter = round(obj.yRes.b1);
%             if ( xCenter<0 || yCenter<0)
%                 a=1;
%             end
            
%             [xData yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, xCenter, yCenter, appData.options.avgWidth);
% 
%             s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint, 'Lower', lower, 'Upper', upper);
%             f = fittype('Ax*exp(-(x-x0)^2/(2*sigmaX^2))', 'coefficients', {'Ax', 'x0', 'sigmaX'}, 'independent', 'x', 'dependent', 'y', 'options', s);
%             [obj.xRes, obj.xGof, obj.xOutput] = fit(x', xData', f); 
%             s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint, 'Lower', lower, 'Upper', upper);
%             f = fittype('Ay*exp(-(y-y0)^2/(2*sigmaY^2))', 'coefficients', {'Ay', 'y0', 'sigmaY'}, 'independent', 'y', 'dependent', 'x', 'options', s);
%             [obj.yRes, obj.yGof, obj.yOutput] = fit(y', yData', f); 
%             xCenter = round(obj.xRes.x0);
%             yCenter = round(obj.yRes.y0);            
%             if ( xCenter<0 || yCenter<0)
%                 a=1;
%             end
            
            [xData, yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, xCenter, yCenter, appData.options.avgWidth);
            x = [1 : w] + x0-1;
            y = [1 : h] + y0-1;
            xData = LowPassFilter(xData, 15, appData.options.avgWidth);
            yData = LowPassFilter(yData, 15, appData.options.avgWidth);
            [ startPoint, lower, upper ] = calcInitialGuess( x, xData, appData.options.avgWidth );
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint, 'Lower', lower, 'Upper', upper);
            f = fittype('Ax*exp(-(x-x0)^2/(2*sigmaX^2))+Cx', 'coefficients', {'Ax', 'x0', 'sigmaX', 'Cx'}, 'independent', 'x', 'dependent', 'y', 'options', s);
            [obj.xRes, obj.xGof, obj.xOutput] = fit(x', xData', f);
            [ startPoint, lower, upper ] = calcInitialGuess( y, yData, appData.options.avgWidth );
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint, 'Lower', lower, 'Upper', upper);
            f = fittype('Ay*exp(-(y-y0)^2/(2*sigmaY^2))+Cy', 'coefficients', {'Ay', 'y0', 'sigmaY', 'Cy'}, 'independent', 'y', 'dependent', 'x', 'options', s);
            [obj.yRes, obj.yGof, obj.yOutput] = fit(y', yData', f); 

            [xData, yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, xCenter, yCenter, appData.options.avgWidth);
            x = [1 : w] + x0-1;
            y = [1 : h] + y0-1;
            lower = [0 0 0 -10];
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [ obj.xRes.Ax obj.xRes.x0 abs(obj.xRes.sigmaX) obj.xRes.Cx ], 'Lower', lower);
            f = fittype('Ax*exp(-(x-x0)^2/(2*sigmaX^2))+Cx', 'coefficients', {'Ax', 'x0', 'sigmaX', 'Cx'}, 'independent', 'x', 'dependent', 'y', 'options', s);
            [obj.xRes, obj.xGof, obj.xOutput] = fit(x', xData', f);
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [ obj.yRes.Ay obj.yRes.y0 abs(obj.yRes.sigmaY) obj.yRes.Cy ], 'Lower', lower);
            f = fittype('Ay*exp(-(y-y0)^2/(2*sigmaY^2))+Cy', 'coefficients', {'Ay', 'y0', 'sigmaY', 'Cy'}, 'independent', 'y', 'dependent', 'x', 'options', s);
            [obj.yRes, obj.yGof, obj.yOutput] = fit(y', yData', f); 
            
            i = 0;
            while ( ( (abs(round(obj.xRes.x0)-xCenter) > 1) || (abs(round(obj.yRes.y0)-yCenter) > 1) ) && i < 10 )
                i = i + 1;
                xCenter = round(obj.xRes.x0);
                yCenter = round(obj.yRes.y0);
%                 appData.data.xPosMax = round(xRes.x0);
%                 appData.data.yPosMax = round(yRes.y0);

                [xData yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, xCenter, yCenter, appData.options.avgWidth);
                s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [obj.xRes.Ax obj.xRes.x0 abs(obj.xRes.sigmaX) obj.xRes.Cx]);
                f = fittype('Ax*exp(-(x-x0)^2/(2*sigmaX^2))+Cx', 'coefficients', {'Ax', 'x0', 'sigmaX', 'Cx'}, 'independent', 'x', 'dependent', 'y', 'options', s);
                [obj.xRes, obj.xGof, obj.xOutput] = fit(x', xData', f); 
                s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [obj.yRes.Ay obj.yRes.y0 abs(obj.yRes.sigmaY) obj.yRes.Cy]);
                f = fittype('Ay*exp(-(y-y0)^2/(2*sigmaY^2))+Cy', 'coefficients', {'Ay', 'y0', 'sigmaY', 'Cy'}, 'independent', 'y', 'dependent', 'x', 'options', s);
                [obj.yRes, obj.yGof, obj.yOutput] = fit(y', yData', f); 
            end
%             xCenter = round(xRes.x0);
%             yCenter = round(yRes.y0);
           
            % set fit 1D params
%             obj.ODx = xRes.Ax;
%             obj.ODy = yRes.Ay;
%             obj.x0 = xRes.x0;
%             obj.y0 = yRes.y0;
%             obj.sigmaX = xRes.sigmaX;
%             obj.sigmaY = yRes.sigmaY;
%             obj.Cx = xRes.Cx;
%             obj.Cy = yRes.Cy;
            
            %set fit params
            obj.xCenter = round(obj.xRes.x0); % should be indexes (integers)
            obj.yCenter= round(obj.yRes.y0);
            obj.xUnitSize = obj.xRes.sigmaX;
            obj.yUnitSize = obj.yRes.sigmaY;
            obj.maxVal = mean([obj.xRes.Ax obj.yRes.Ay]);
           
           % calc ROI size (use ROIUnits.m) - MUST be after fit
           [obj.ROILeft obj.ROITop obj.ROIRight obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, obj);
           
           obj.atomsNo = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
               [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1); 
%            obj.atomsNo = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
%                [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1); 
           
            
           obj.xData = xData;
           obj.xStart = x0;
           obj.yData = yData;
           obj.yStart = y0;
           
           % last 
           appData.data.fits{appData.consts.fitTypes.oneDGaussian} = obj;
           appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);
           
           function  [ startPoint, lower, upper ] = calcInitialGuess( x, xData, filterSigma )
               if filterSigma==0
                   xData = xData( round(1.5): (end-round(1.5)) ); % omit points which are at the edges, to avoid filter edge effects
               else
                   xData = xData( round(1.5*filterSigma): (end-round(1.5*filterSigma)) ); % omit points which are at the edges, to avoid filter edge effects
               end
               startPoint =[]; % {'Ax', 'x0', 'sigmaX', 'Cx'}
%                xDataS = smooth(xData, round(length(xData)/5) )'; %
%                smoothing dropped, since data is passed through a low pass
%                filter
               [temp, Ix] = max(xData);
               CStartPoint = mean( xData(xData<(mean(xData)+std(xData))) ); %background start point is the mean value of all points which are no further than 1 sigma from mean value of xData
               AStartPoint = max(xData) - CStartPoint; %Gaussian amplitude
               xData(xData<(mean(xData)+std(xData))) = 0; % set points outside of peak to zero - to evaluate sigma
               sigmaStartPoint = find(xData, 1, 'last' )-find(xData, 1 );

               startPoint = [ AStartPoint x(Ix)+round(1.5*filterSigma) sigmaStartPoint CStartPoint]; % [Amplitude, x0, sigma, background]
               lower(1) = 0.8 * startPoint(1); 
               upper(1) = 1.5 * startPoint(1);
               lower(2) = startPoint(2) - 10;
               upper(2) = startPoint(2) + 10;
               lower(3) = 0.25 * startPoint(3);
               upper(3) = 2*startPoint(3);
               lower(4) = startPoint(4) - abs(startPoint(4));
               upper(4) = startPoint(4) + abs(startPoint(4));

           end
           
       end
       
       function normalizedROI = getNormalizedROI(obj, pic, x, y) % return normalized ROI (to the fitting constant)
           normalizedROI = pic(y, x) - mean([obj.xRes.Cx obj.yRes.Cy]);
       end
       
       function normalizedROI = getTheoreticalROI(obj, pic, x, y)
           [X, Y] = meshgrid(x, y);
           normalizedROI = mean([obj.xRes.Ax obj.yRes.Ay]) * exp( -0.5 * ((X-obj.xRes.x0).^2 / obj.xRes.sigmaX^2 + (Y-obj.yRes.y0).^2 / obj.yRes.sigmaY^2) );
       end
       
       function normalizedPic = normalizePic(obj, pic)
           normalizedPic = (pic-mean([obj.xRes.Cx obj.yRes.Cy]))/obj.maxVal;
       end
       
       function [xFit yFit] = getXYFitVectors(obj, x, y)
           xFit = obj.xRes.Ax * exp( -0.5 * (x-obj.xRes.x0).^2 / obj.xRes.sigmaX^2 ) + obj.xRes.Cx;
           yFit = obj.yRes.Ay * exp( -0.5 * (y-obj.yRes.y0).^2 / obj.yRes.sigmaY^2 ) + obj.yRes.Cy;
       end
       
       function  plotFitResults(obj, appData)  % plots the text           
%            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getPic();
           chipStart = appData.data.camera.chipStart;
           cla(appData.ui.pFitResults)
%            text( 10, 190, ['Atoms Num: ' num2str(obj.atomsNo/1e6) '*10^6'], 'fontSize', 20);
%            text( 50, 140, {['x_0 = ' num2str((obj.xCenter-x0+1) * appData.data.camera.xPixSz * 1000) ' mm'], ...
%                ['y_0 = ' num2str((obj.yCenter-y0+1) * appData.data.camera.yPixSz * 1000) ' mm']}, ...
%                'fontsize', 12);
           text( 0, 1, ['Atoms Num: ' addCommas(obj.atomsNo)], 'fontSize', 20, 'Parent', appData.ui.pFitResults);
           text( 0, 0.8, 'fit function = A_x * {\ite}^{-(x-x_0)^2 / 2\sigma_x^2}+C_x + A_y * {\ite}^{-(y-y_0)^2 / 2\sigma_y^2}+C_y',...
               'fontsize', 12, 'Parent', appData.ui.pFitResults);
           text( 0, 0.29, {['A_x = ' num2str(obj.xRes.Ax)], ...
               ['x_0 = ' num2str((obj.xRes.x0) * appData.data.camera.xPixSz * 1000) ' mm;'],...
               ['\sigma_x = ' num2str(obj.xRes.sigmaX * appData.data.camera.xPixSz * 1e6) ' um;'], ... %['\sigma_x = ' num2str(obj.xRes.sigmaX * appData.data.camera.xPixSz * 1000) ' mm'], ...
               ['C_x = ' num2str(obj.xRes.Cx) ';'], ...
               ['R^2 (x) = ' num2str(obj.xGof.rsquare) ';'],...
               ['X/Y size aspect ratio = ' num2str(obj.xRes.sigmaX./obj.yRes.sigmaY)]}, ...
               'fontsize', 12, 'Parent', appData.ui.pFitResults);
           text( 0.3, 0.35, {['A_y = ' num2str(obj.yRes.Ay)], ...
               ['y_0 = ' num2str((obj.yRes.y0-chipStart) * appData.data.camera.yPixSz * 1000) ' mm (from the chip)'], ...
               ['\sigma_y = ' num2str(obj.yRes.sigmaY * appData.data.camera.yPixSz * 1e6) ' um'], ... %['\sigma_y = ' num2str(obj.yRes.sigmaY * appData.data.camera.yPixSz * 1000) ' mm'], ...
               ['C_y = ' num2str(obj.yRes.Cy)], ...
               ['R^2 (y)= ' num2str(obj.yGof.rsquare)]}, ...
               'fontsize', 12, 'Parent', appData.ui.pFitResults);
       end
       
   end
end 
