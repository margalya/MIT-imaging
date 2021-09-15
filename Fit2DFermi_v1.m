classdef Fit2DFermi < FitTypes
    %FIT1DGAUSSIAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties ( Constant = true )
        ID = 'Fit2DFermi';
    end
    properties (SetAccess = private )
        xRes1DG = []; %1DGaussian x fit results
        yRes1DG = []; %1DGaussian y fit results
        OD = -1;
        x0 = -1;
        y0 = -1;
        Rx = -1;
        Ry = -1;
        q = -1;
        C = -1;
        ToverTF = [];
%         ToverTFError = [];
        fval = [];
        exitflag = [];
        output = [];
    end
    
    methods
        function appData = analyze(obj, appData) % do the analysis
            
            % 1D Gaussian fit
            fitObj = appData.data.fits{appData.consts.fitTypes.oneDGaussian};
            if ( isempty( fitObj.xRes ) )
                tmpFitType = appData.data.fitType;
                appData.data.fitType = appData.consts.fitTypes.oneDGaussian;
                appData = appData.data.fits{appData.consts.fitTypes.oneDGaussian}.analyze(appData);
                fitObj = appData.data.fits{appData.consts.fitTypes.oneDGaussian};
                appData.data.fitType = tmpFitType;
            end
            obj.xRes1DG = fitObj.xRes; % assign 1DGaussian x fit results, for use in later functions (getNormalizedROI, normalizePic, getXYFitVectors).
            obj.yRes1DG = fitObj.yRes; % 1DGaussian y fit results
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);  % the whole pic
            [pic, x0, y0] = appData.data.plots{appData.consts.plotTypes.ROI}.getAnalysisPic(appData); %#ok<PROPLC> % only ROI
            [h, w] = size(pic);
            
            % 2D fit
            firstGuess = [0.5*(fitObj.xRes.Ax+fitObj.yRes.Ay) 0 fitObj.xRes.x0 1.8*fitObj.xRes.sigmaX fitObj.yRes.y0 ...
                1.8*fitObj.yRes.sigmaY 0.5*(fitObj.xRes.Cx+fitObj.yRes.Cy)];
            binW = appData.options.avgWidth;
            if ( numel(pic) > 250^2 )
                binnedPic = binning(pic, binW);
                [h, w] = size(binnedPic);
                [X, Y] = meshgrid( (1 : binW : binW*w) +x0-1, (1 : binW : binW*h) +y0-1);% - appData.data.camera.chipStart);
                [fitRes, obj.fval, obj.exitflag, obj.output] = ...
                    fminsearch(@(p) fit2DFermi_scalar( p, X, Y, binnedPic), firstGuess, optimset('TolX', 1e-4, 'MaxFunEvals', 2000) );
            else
                [X, Y] = meshgrid((1  : w) +x0-1, (1  : h) +y0-1);
                [fitRes, obj.fval, obj.exitflag, obj.output] = ...
                    fminsearch(@(p) fit2DFermi_scalar( p, X, Y, pic), firstGuess, optimset('TolX', 1e-4, 'MaxFunEvals', 2000) );
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % set 2D params
            obj.OD = fitRes(1);
            obj.q = fitRes(2);
            obj.x0 = fitRes(3);
            obj.Rx = fitRes(4);
            obj.y0 = fitRes(5);
            obj.Ry = fitRes(6);
            obj.C = fitRes(7);
            obj.ToverTF = (-6*polylog(3,-exp(obj.q)))^(-1/3);
%             obj.ToverTFError = (-6*polylog(3,-exp(obj.q - obj.conf(2))))^(-1/3) - (-6*polylog(3,-exp(obj.q + obj.conf(2))))^(-1/3);           
            
            %%% test
            picFit = (Fermi2D( obj.OD, obj.q, X - obj.x0, obj.Rx, Y - obj.y0, obj.Ry ) + obj.C);
            picFit = reshape(picFit,size(binnedPic));
            figure;
            imagesc(picFit)
            %%%
            
            %set fit params
            obj.xCenter = round(obj.x0); % should be indexes (integers)
            obj.yCenter= round(obj.y0);
            obj.xUnitSize = obj.Rx;
            obj.yUnitSize = obj.Ry;
            obj.maxVal = obj.OD;
           
           % calc ROI size (use ROIUnits.m) - MUST be after fit
           [obj.ROILeft, obj.ROITop, obj.ROIRight, obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, obj);
           appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);
           [pic, x0, y0] = appData.data.plots{appData.consts.plotTypes.ROI}.getAnalysisPic(appData);
%            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
           
            obj.atomsNo = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
               [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1); 
%            obj.atomsNo = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
%                [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1); 
           
            [xData, yData] = appData.data.plots{appData.data.plotType }.getXYDataVectors(obj.xCenter, obj.yCenter, binW);

            obj.xData = xData;
            obj.xStart = x0; %#ok<PROPLC>
            obj.yData = yData;
            obj.yStart = y0; %#ok<PROPLC>

            % last
            appData.data.fits{appData.consts.fitTypes.Fermi2D} = obj;
%            appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);

            
            
        end
        
       function normalizedROI = getNormalizedROI(obj, pic, x, y) % return normalized ROI (to the fitting constant)
           normalizedROI = pic(y, x) - mean([obj.xRes1DG.Cx obj.yRes1DG.Cy]);
       end
       
       function normalizedROI = getTheoreticalROI(obj, pic, x, y)
           [X, Y] = meshgrid(x, y);
           normalizedROI = mean([obj.xRes1DG.Ax obj.yRes1DG.Ay]) * exp( -0.5 * ((X-obj.xRes1DG.x0).^2 / obj.xRes1DG.sigmaX^2 + (Y-obj.yRes1DG.y0).^2 / obj.yRes1DG.sigmaY^2) );
       end
       
       function normalizedPic = normalizePic(obj, pic)
           normalizedPic = (pic-mean([obj.xRes1DG.Cx obj.yRes1DG.Cy]))/obj.maxVal;
       end
        
        function [xFit, yFit] = getXYFitVectors(obj, x, y)
%             Fermi2D( n0, q, x, Rx, y, Ry ) + c
            xFit = (Fermi2D( obj.OD, obj.q, x - obj.x0, obj.Rx, 0*obj.y0, obj.Ry ) + obj.C)';
            yFit = (Fermi2D( obj.OD, obj.q, 0*obj.x0 , obj.Rx, y - obj.y0, obj.Ry ) + obj.C)';
%             [xFit2, yFit2] = appData.data.fits{appData.consts.fitTypes.twoDGaussian}.getXYFitVectors(x, y);
%             xFit2 = obj.OD * exp( -0.5 * (x-obj.x0).^2 / obj.sigmaX^2 ) + obj.C;
%             yFit2 = obj.OD * exp( -0.5 * (y-obj.y0).^2 / obj.sigmaY^2 ) + obj.C;
            
%            xFit = [xFit2; xFit1];
%            yFit = [yFit2; yFit1];
        end
        
        function  plotFitResults(obj, appData)  % plots the text
            chipStart = appData.data.camera.chipStart;
            text( 10, 220, ['Atoms Num: ' addCommas(obj.atomsNo)], 'fontSize', 20);
            text( 10, 190, 'fit function = Fermi2D( n0, q, x, Rx, y, Ry ) + c', 'fontsize', 12);
            text( 50, 80, {['n_0 = ' num2str(obj.OD)], ... %                ['A_G = ' num2str(obj.rRes1DG.A)], ...
                ['q = ' num2str(obj.q)], ... %' +/- ' num2str(obj.conf(2))
                ['T/T_F = ' num2str(obj.ToverTF)], ... %' +/- ' num2str(obj.ToverTFError)
                ['(x_0, y_0) = (' num2str( obj.x0 * appData.data.camera.xPixSz * 1000) '; '...
                num2str((obj.y0-chipStart) * appData.data.camera.yPixSz * 1000) ') mm'], ...
                ['(R_x, R_y) = ' num2str(obj.Rx * appData.data.camera.xPixSz * 1e6) '; '...
                num2str(obj.Ry * appData.data.camera.yPixSz * 1e6) ' um'],...
                ['C = ' num2str(obj.C)]},...
                'fontsize', 12);
        end
        
    end
end

function ret = fit2DFermi_scalar( p, X, Y, pic )
% p(1) = n0;
% p(2) = q ;
% p(3) = x0;
% p(4) = Rx;
% p(5) = y0;
% p(6) = Ry;
% p(7) = C;
% 1D function: 'Fermi1D( n0, q, x, R ) + C'
% 2D function:  Fermi2D( n0, q, x, Rx, y, Ry )
ret = Fermi2D( p(1), p(2), (X - p(3)), p(4), (Y - p(5)), p(6) ) + p(7) - pic(:); %returns a vector for a matrix input, but is seems that it dosen't matter as we just need the sum
ret = sum(ret.^2);
end



