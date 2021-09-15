classdef FitSternGerlach < FitTypes
    %FITSTERNGERLACH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties ( Constant = true )
        ID = 'FitSternGerlach';
    end
    properties (SetAccess = private )
        
        xCenter1 = [];
        yCenter1 = [];
        xRes1 = [];
        xGof1 = [];
        xOutput1 = [];
        yRes1 = [];
        yGof1 = [];
        yOutput1 = [];
        
        OD1 = -1;
        x01 = -1;
        y01 = -1;
        sigmaX1 = -1;
        sigmaY1 = -1;
        C1 = [];
        
        fval1 = [];
        exitflag1 = [];
        output1 = [];
        
        xCenter2 = [];
        yCenter2 = [];
        xRes2 = [];
        xGof2 = [];
        xOutput2 = [];
        yRes2 = [];
        yGof2 = [];
        yOutput2 = [];
        
        OD2 = -1;
        x02 = -1;
        y02 = -1;
        sigmaX2 = -1;
        sigmaY2 = -1;
        C2 = [];
        
        fval2 = [];
        exitflag2 = [];
        output2 = [];
        
        atomsNo1 = [];
        atomsNo2 = [];
        fAtomNo = 0;
        mF1 = [];
    end
    
    methods
        function appData = analyze2(obj, appData)
            %            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
            %            [h w] = size(pic);
            
            tmpPlotType = appData.data.plotType;
            
            
            % FIRST FIT
            appData.data.fits{appData.consts.fitTypes.onlyMaximum} = FitOnlyMax;
            appData.data.fits{appData.consts.fitTypes.oneDGaussian} = Fit1DGaussian;
            appData.data.fits{appData.consts.fitTypes.twoDGaussian} = Fit2DGaussian;
            % only maximum and ROI
            appData.data.fitType = appData.consts.fitTypes.onlyMaximum;
            appData = appData.data.fits{appData.data.fitType}.analyze(appData);
            appData.data.plotType = appData.consts.plotTypes.ROI;
            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
            [h w] = size(pic);
            %2D fit in ROI
            appData.data.fitType = appData.consts.fitTypes.twoDGaussian;
            appData = appData.data.fits{appData.data.fitType}.analyze(appData);
            fitObj_first = appData.data.fits{appData.consts.fitTypes.twoDGaussian};
            
            %SECOND FIT
            appData.data.fits{appData.consts.fitTypes.onlyMaximum} = FitOnlyMax;
            appData.data.fits{appData.consts.fitTypes.oneDGaussian} = Fit1DGaussian;
            appData.data.fits{appData.consts.fitTypes.twoDGaussian} = Fit2DGaussian;
            appData.data.plotType = appData.consts.plotTypes.absorption;
            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
            [h w] = size(pic);
            tmpAbsorptionPic = pic;
            pic([fitObj_first.ROITop : fitObj_first.ROIBottom] -y0, [fitObj_first.ROILeft : fitObj_first.ROIRight] ) = nan;
            pic(1 : 2*appData.data.camera.chipStart, :) = nan; %Erase the upper area of the picture, of height 2*chipStart, to avoid false maxima idenification
%             pic(1 : 1000, :) = 0; %Erase a larger area - to avoid identifying mF=0 atoms (for certain pictures).
            ODNoise = nanstd(pic(:)); %OD noise is calcualted here, before pic is changed
            OD_DC_Level = nanmean(pic(:));
            newPic = zeros(h+y0-1, w+x0-1) ;
            newPic(y0:end, x0:end) = pic;
            appData = appData.data.plots{appData.consts.plotTypes.absorption}.setPic(appData, newPic);
            % only maximum and ROI
            appData.data.fitType = appData.consts.fitTypes.onlyMaximum;
            appData = appData.data.fits{appData.data.fitType}.analyze(appData);
            appData.data.plotType = appData.consts.plotTypes.ROI;
            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
            [h w] = size(pic);
            %2D fit in ROI
            try 
            appData.data.fits{appData.consts.fitTypes.oneDGaussian} = Fit1DGaussian;
            appData.data.fitType = appData.consts.fitTypes.twoDGaussian;
            appData = appData.data.fits{appData.data.fitType}.analyze(appData);
            fitObj_second = appData.data.fits{appData.consts.fitTypes.twoDGaussian};
            twoClouds = 1; %fit obtained two clouds
            catch
                fitObj_second.OD = ODNoise;
                fitObj_second.C = OD_DC_Level;
                cloudArea = 2.0510e-10; %estimation of cloud area, taken from data - this needs to be fixed every time running the code
                fitObj_second.atomsNo = round(fitObj_second.OD * cloudArea / appData.consts.scatcross0 *2*pi / 2); % half of the evaluated minimum detecable number of atoms, for the given sigma
                % we take half since we claim - the undetected number of atoms is 0 < N < N_minimum
                fitObj_second.x0 = nan;
                fitObj_second.y0 = nan;
                fitObj_second.sigmaX = nan;
                fitObj_second.sigmaY = nan;
                twoClouds = 0; %fit did not obtain two clouds
            end
            
            %reset absorption
            appData.data.plotType = appData.consts.plotTypes.absorption;
            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
            [h w] = size(pic);
            newPic = zeros(h+y0-1, w+x0-1) ;
            newPic(y0:end, x0:end) = tmpAbsorptionPic;
            appData = appData.data.plots{appData.consts.plotTypes.absorption}.setPic(appData, newPic);
            pic = tmpAbsorptionPic;
            
            %set params
            if twoClouds %two clouds were fitted - work as usual
                if ( fitObj_first.y0 > fitObj_second.y0 ) %first is mF2
                    obj.OD2 = fitObj_first.OD;
                    obj.x02 = fitObj_first.x0;
                    obj.y02 = fitObj_first.y0;
                    obj.sigmaX2 = fitObj_first.sigmaX;
                    obj.sigmaY2 = fitObj_first.sigmaY;
                    obj.C2 = fitObj_first.C;
                    
                    obj.OD1 = fitObj_second.OD;
                    obj.x01 = fitObj_second.x0;
                    obj.y01 = fitObj_second.y0;
                    obj.sigmaX1 = fitObj_second.sigmaX;
                    obj.sigmaY1 = fitObj_second.sigmaY;
                    obj.C1 = fitObj_second.C;
                    
                    obj.xCenter = round(obj.x02); % should be indexes (integers)
                    obj.yCenter= round(obj.y02);
                    obj.xUnitSize = obj.sigmaX2;
                    obj.yUnitSize = obj.sigmaY2;
                    obj.maxVal = obj.OD2;
                    
                    obj.atomsNo2 = fitObj_first.atomsNo;
                    obj.atomsNo1 = fitObj_second.atomsNo;
                    %                 obj.atomsNo2 = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, fitObj_first, pic, ...
                    %                     [fitObj_first.ROILeft : fitObj_first.ROIRight] - x0+1, ...
                    %                     [fitObj_first.ROITop : fitObj_first.ROIBottom] - y0+1);
                    %                 obj.atomsNo1 = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, fitObj_second, pic, ...
                    %                     [fitObj_second.ROILeft : fitObj_second.ROIRight] - x0+1, ...
                    %                     [fitObj_second.ROITop : fitObj_second.ROIBottom] - y0+1);
                    
                else
                    obj.OD2 = fitObj_second.OD;
                    obj.x02 = fitObj_second.x0;
                    obj.y02 = fitObj_second.y0;
                    obj.sigmaX2 = fitObj_second.sigmaX;
                    obj.sigmaY2 = fitObj_second.sigmaY;
                    obj.C2 = fitObj_second.C;
                    
                    obj.OD1 = fitObj_first.OD;
                    obj.x01 = fitObj_first.x0;
                    obj.y01 = fitObj_first.y0;
                    obj.sigmaX1 = fitObj_first.sigmaX;
                    obj.sigmaY1 = fitObj_first.sigmaY;
                    obj.C1 = fitObj_first.C;
                    
                    obj.xCenter = round(obj.x01); % should be indexes (integers)
                    obj.yCenter= round(obj.y01);
                    obj.xUnitSize = obj.sigmaX1;
                    obj.yUnitSize = obj.sigmaY1;
                    obj.maxVal = obj.OD1;
                    
                    obj.atomsNo1 = fitObj_first.atomsNo;
                    obj.atomsNo2 = fitObj_second.atomsNo;
                    %                 obj.atomsNo1 = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, fitObj_first, pic, ...
                    %                     [fitObj_first.ROILeft : fitObj_first.ROIRight] - x0+1, ...
                    %                     [fitObj_first.ROITop : fitObj_first.ROIBottom] - y0+1);
                    %                 obj.atomsNo2 = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, fitObj_second, pic, ...
                    %                     [fitObj_second.ROILeft : fitObj_second.ROIRight] - x0+1, ...
                    %                     [fitObj_second.ROITop : fitObj_second.ROIBottom] - y0+1);
                    
                end
            else %only one cloud was fitted - check height to determine which mF is it
                
                if ( fitObj_first.y0 >= 1400 ) %first is mF2
                    obj.OD2 = fitObj_first.OD;
                    obj.x02 = fitObj_first.x0;
                    obj.y02 = fitObj_first.y0;
                    obj.sigmaX2 = fitObj_first.sigmaX;
                    obj.sigmaY2 = fitObj_first.sigmaY;
                    obj.C2 = fitObj_first.C;
                    
                    obj.OD1 = fitObj_second.OD;
                    obj.x01 = fitObj_second.x0;
                    obj.y01 = fitObj_second.y0;
                    obj.sigmaX1 = fitObj_second.sigmaX;
                    obj.sigmaY1 = fitObj_second.sigmaY;
                    obj.C1 = fitObj_second.C;
                    
                    obj.xCenter = round(obj.x02); % should be indexes (integers)
                    obj.yCenter= round(obj.y02);
                    obj.xUnitSize = obj.sigmaX2;
                    obj.yUnitSize = obj.sigmaY2;
                    obj.maxVal = obj.OD2;
                    
                    obj.atomsNo2 = fitObj_first.atomsNo;
                    obj.atomsNo1 = fitObj_second.atomsNo;
                    %                 obj.atomsNo2 = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, fitObj_first, pic, ...
                    %                     [fitObj_first.ROILeft : fitObj_first.ROIRight] - x0+1, ...
                    %                     [fitObj_first.ROITop : fitObj_first.ROIBottom] - y0+1);
                    %                 obj.atomsNo1 = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, fitObj_second, pic, ...
                    %                     [fitObj_second.ROILeft : fitObj_second.ROIRight] - x0+1, ...
                    %                     [fitObj_second.ROITop : fitObj_second.ROIBottom] - y0+1);
                    
                elseif ( fitObj_first.y0 < 1400 )
                    obj.OD2 = fitObj_second.OD;
                    obj.x02 = fitObj_second.x0;
                    obj.y02 = fitObj_second.y0;
                    obj.sigmaX2 = fitObj_second.sigmaX;
                    obj.sigmaY2 = fitObj_second.sigmaY;
                    obj.C2 = fitObj_second.C;
                    
                    obj.OD1 = fitObj_first.OD;
                    obj.x01 = fitObj_first.x0;
                    obj.y01 = fitObj_first.y0;
                    obj.sigmaX1 = fitObj_first.sigmaX;
                    obj.sigmaY1 = fitObj_first.sigmaY;
                    obj.C1 = fitObj_first.C;
                    
                    obj.xCenter = round(obj.x01); % should be indexes (integers)
                    obj.yCenter= round(obj.y01);
                    obj.xUnitSize = obj.sigmaX1;
                    obj.yUnitSize = obj.sigmaY1;
                    obj.maxVal = obj.OD1;
                    
                    obj.atomsNo1 = fitObj_first.atomsNo;
                    obj.atomsNo2 = fitObj_second.atomsNo;
                    %                 obj.atomsNo1 = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, fitObj_first, pic, ...
                    %                     [fitObj_first.ROILeft : fitObj_first.ROIRight] - x0+1, ...
                    %                     [fitObj_first.ROITop : fitObj_first.ROIBottom] - y0+1);
                    %                 obj.atomsNo2 = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, fitObj_second, pic, ...
                    %                     [fitObj_second.ROILeft : fitObj_second.ROIRight] - x0+1, ...
                    %                     [fitObj_second.ROITop : fitObj_second.ROIBottom] - y0+1);
                    
                end
                
            end
            
            obj.atomsNo=obj.atomsNo2+obj.atomsNo1;
            obj.mF1 = obj.atomsNo1/obj.atomsNo;
            
            appData.data.plotType = tmpPlotType;
            appData.data.fitType = appData.consts.fitTypes.SG;
            appData.data.fits{appData.consts.fitTypes.onlyMaximum} = FitOnlyMax;
            appData.data.fits{appData.consts.fitTypes.oneDGaussian} = Fit1DGaussian;
            appData.data.fits{appData.consts.fitTypes.twoDGaussian} = Fit2DGaussian;
            
            [obj.ROILeft obj.ROITop obj.ROIRight obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, obj);
            
            [xData yData] = appData.data.plots{appData.data.plotType }.getXYDataVectors(obj.xCenter, obj.yCenter, appData.options.avgWidth);
            
            obj.xData = xData;
            obj.xStart = x0;
            obj.yData = yData;
            obj.yStart = y0;
            
            % last
            appData.data.fits{appData.consts.fitTypes.SG} = obj;
            appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);
        end
        
        function appData = analyze(obj, appData) % do the analysis
            appData = obj.analyze2(appData);
            return
            
            %             % 2 1D fits - only Y
            %             fitObj = appData.data.fits{appData.consts.fitTypes.twoOneDGaussianOnlyY};
            %             if ( isempty( fitObj.yRes ) )
            %                 tmpFitType = appData.data.fitType;
            %                 appData.data.fitType = appData.consts.fitTypes.twoOneDGaussianOnlyY;
            %                 appData = appData.data.fits{appData.consts.fitTypes.twoOneDGaussianOnlyY}.analyze(appData);
            %                 fitObj = appData.data.fits{appData.consts.fitTypes.twoOneDGaussianOnlyY};
            %                 appData.data.fitType = tmpFitType;
            %            end
            
            
            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
            [h w] = size(pic);
            
            % 1D fit - FIRST - mF2
            % binning and max
            binnedData = binning ( pic, appData.options.avgWidth*2+1);
            [maxes, indexes] = max(binnedData);                     % find maximum
            [maxValue, xPosMax] = max(maxes);
            yPosMax = indexes(xPosMax);
            % center
            %             xCenter = (appData.options.avgWidth*2+1) * (xPosMax ) + x0-appData.options.avgWidth;
            %            yCenter = (appData.options.avgWidth*2+1) * (yPosMax ) + y0-appData.options.avgWidth;
            obj.xCenter1 =(appData.options.avgWidth*2+1) * (xPosMax ) + x0-appData.options.avgWidth;
            obj.yCenter1 = (appData.options.avgWidth*2+1) * (yPosMax ) + y0-appData.options.avgWidth;
            
            % fitting
            [xData yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, obj.xCenter1, obj.yCenter1, appData.options.avgWidth);
            x = [1 : w] + x0-1;
            y = [1 : h] + y0-1;
            [obj.xRes1, obj.xGof1, obj.xOutput1] = fit(x', xData', 'gauss1');
            [obj.yRes1, obj.yGof1, obj.yOutput1] = fit(y', yData', 'gauss1');
            obj.xCenter1 = round(obj.xRes1.b1);
            obj.yCenter1 = round(obj.yRes1.b1);
            
            [xData yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, obj.xCenter1, obj.yCenter1, appData.options.avgWidth);
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [obj.xRes1.a1 obj.xRes1.b1 obj.xRes1.c1/sqrt(2)]);%[1 width/2 width/10]);
            f = fittype('Ax*exp(-(x-x0)^2/(2*sigmaX^2))', 'coefficients', {'Ax', 'x0', 'sigmaX'}, 'independent', 'x', 'dependent', 'y', 'options', s);
            [obj.xRes1, obj.xGof1, obj.xOutput1] = fit(x', xData', f);
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [obj.yRes1.a1 obj.yRes1.b1 obj.yRes1.c1/sqrt(2)]);%[1 width/2 width/10]);
            f = fittype('Ay*exp(-(y-y0)^2/(2*sigmaY^2))', 'coefficients', {'Ay', 'y0', 'sigmaY'}, 'independent', 'y', 'dependent', 'x', 'options', s);
            [obj.yRes1, obj.yGof1, obj.yOutput1] = fit(y', yData', f);
            obj.xCenter1 = round(obj.xRes1.x0);
            obj.yCenter1 = round(obj.yRes1.y0);
            
            [xData yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, obj.xCenter1, obj.yCenter1, appData.options.avgWidth);
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [obj.xRes1.Ax obj.xRes1.x0 obj.xRes1.sigmaX 0], 'Lower', [0 -Inf 0 -Inf]);
            f = fittype('Ax*exp(-(x-x0)^2/(2*sigmaX^2))+Cx', 'coefficients', {'Ax', 'x0', 'sigmaX', 'Cx'}, 'independent', 'x', 'dependent', 'y', 'options', s);
            [obj.xRes1, obj.xGof1, obj.xOutput1] = fit(x', xData', f);
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [obj.yRes1.Ay obj.yRes1.y0 obj.yRes1.sigmaY 0], 'Lower', [0 -Inf 0 -Inf]);
            f = fittype('Ay*exp(-(y-y0)^2/(2*sigmaY^2))+Cy', 'coefficients', {'Ay', 'y0', 'sigmaY', 'Cy'}, 'independent', 'y', 'dependent', 'x', 'options', s);
            [obj.yRes1, obj.yGof1, obj.yOutput1] = fit(y', yData', f);
            obj.xCenter1 = round(obj.xRes1.x0);
            obj.yCenter1 = round(obj.yRes1.y0);
            
            % 2D fit - FIRST - mF2
            appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);
            %              [obj.ROILeft obj.ROITop obj.ROIRight obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, obj);
            [pic x0 y0] = appData.data.plots{appData.consts.plotTypes.ROI}.getAnalysisPic(appData); % only ROI
            [h w] = size(pic);
            firstGuess = [mean([obj.xRes1.Ax obj.yRes1.Ay]) obj.xRes1.x0 obj.yRes1.y0 ...
                obj.xRes1.sigmaX obj.yRes1.sigmaY mean([obj.xRes1.Cx obj.yRes1.Cy])];
            if ( numel(pic) > 350^2 )
                binnedPic = binning(pic, 2);
                [h w] = size(binnedPic);
                [X, Y] = meshgrid([1 : 2 : 2*w] +x0-1, [1 : 2 : 2*h] +y0-1);% - appData.data.camera.chipStart);
                [fitRes1, obj.fval1, obj.exitflag1, obj.output1] = ...
                    fminsearch(@(p) fitGauss2D_scalar( p, X, Y, binnedPic), firstGuess, optimset('TolX',1e-6, 'MaxFunEvals', 1000, 'MaxIter', 500) );
            else
                [X, Y] = meshgrid([1  : w] +x0-1, [1  : h] +y0-1);
                [fitRes1, obj.fval1, obj.exitflag1, obj.output1] = ...
                    fminsearch(@(p) fitGauss2D_scalar( p, X, Y, pic), firstGuess, optimset('TolX',1e-6, 'MaxFunEvals', 1000, 'MaxIter', 500) );
            end
            % set 2D params
            obj.OD1 = fitRes1(1);
            obj.x01 = fitRes1(2);
            obj.y01 = fitRes1(3);
            obj.sigmaX1 = fitRes1(4);
            obj.sigmaY1 = fitRes1(5);
            obj.C1 = fitRes1(6);
            
            obj.xCenter = round(obj.x01); % should be indexes (integers)
            obj.yCenter= round(obj.y01);
            obj.xUnitSize = obj.sigmaX1;
            obj.yUnitSize = obj.sigmaY1;
            obj.maxVal = obj.OD1;
            % calc ROI size (use ROIUnits.m) - MUST be after fit
            [obj.ROILeft obj.ROITop obj.ROIRight obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, obj);
            
            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
            [h w] = size(pic);
            obj.fAtomNo = 1;
            obj.atomsNo1 = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
                [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1);
            
            % 1D fit - SECOND - mF1
            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
            [h w] = size(pic);
            
            %             pic = pic(1:obj.yCenter1-y0-round(obj.sigmaY1)*3, :);
            pic([obj.ROITop : obj.ROIBottom] -y0, [obj.ROILeft : obj.ROIRight] ) = 0;
            [h w] = size(pic);
            %binning and max
            binnedData = binning ( pic, appData.options.avgWidth*2+1);
            [maxes, indexes] = max(binnedData);                     % find maximum
            [maxValue, xPosMax] = max(maxes);
            yPosMax = indexes(xPosMax);
            % center
            obj.xCenter2 = (appData.options.avgWidth*2+1) * (xPosMax ) + x0-appData.options.avgWidth;
            obj.yCenter2 = (appData.options.avgWidth*2+1) * (yPosMax ) + y0-appData.options.avgWidth;
            
            % fitting
            [xData yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, obj.xCenter2, obj.yCenter2, appData.options.avgWidth);
            yData([obj.ROITop : obj.ROIBottom]-y0) = 0;
            x = [1 : w] + x0-1;
            y = [1 : h] + y0-1;
            [obj.xRes2, obj.xGof2, obj.xOutput2] = fit(x', xData', 'gauss1');
            [obj.yRes2, obj.yGof2, obj.yOutput2] = fit(y', yData', 'gauss1');
            obj.xCenter2 = round(obj.xRes2.b1);
            obj.yCenter2 = round(obj.yRes2.b1);
            
            [xData yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, obj.xCenter2, obj.yCenter2, appData.options.avgWidth);
            yData([obj.ROITop : obj.ROIBottom]-y0) = 0;
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [obj.xRes2.a1 obj.xRes2.b1 obj.xRes2.c1/sqrt(2)]);%[1 width/2 width/10]);
            f = fittype('Ax*exp(-(x-x0)^2/(2*sigmaX^2))', 'coefficients', {'Ax', 'x0', 'sigmaX'}, 'independent', 'x', 'dependent', 'y', 'options', s);
            [obj.xRes2, obj.xGof2, obj.xOutput2] = fit(x', xData', f);
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [obj.yRes2.a1 obj.yRes2.b1 obj.yRes2.c1/sqrt(2)]);%[1 width/2 width/10]);
            f = fittype('Ay*exp(-(y-y0)^2/(2*sigmaY^2))', 'coefficients', {'Ay', 'y0', 'sigmaY'}, 'independent', 'y', 'dependent', 'x', 'options', s);
            [obj.yRes2, obj.yGof2, obj.yOutput2] = fit(y', yData', f);
            obj.xCenter2 = round(obj.xRes2.x0);
            obj.yCenter2 = round(obj.yRes2.y0);
            
            [xData yData] = appData.data.plots{appData.data.plotType}.getAnalysisXYDataVectors(appData, obj.xCenter2, obj.yCenter2, appData.options.avgWidth);
            yData([obj.ROITop : obj.ROIBottom]-y0) = 0;
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [obj.xRes2.Ax obj.xRes2.x0 obj.xRes2.sigmaX 0], 'Lower', [0 -Inf 0 -Inf]);
            f = fittype('Ax*exp(-(x-x0)^2/(2*sigmaX^2))+Cx', 'coefficients', {'Ax', 'x0', 'sigmaX', 'Cx'}, 'independent', 'x', 'dependent', 'y', 'options', s);
            [obj.xRes2, obj.xGof2, obj.xOutput2] = fit(x', xData', f);
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [obj.yRes2.Ay obj.yRes2.y0 obj.yRes2.sigmaY 0], 'Lower', [0 -Inf 0 -Inf]);
            f = fittype('Ay*exp(-(y-y0)^2/(2*sigmaY^2))+Cy', 'coefficients', {'Ay', 'y0', 'sigmaY', 'Cy'}, 'independent', 'y', 'dependent', 'x', 'options', s);
            [obj.yRes2, obj.yGof2, obj.yOutput2] = fit(y', yData', f);
            obj.xCenter2 = round(obj.xRes2.x0);
            obj.yCenter2 = round(obj.yRes2.y0);
            
            % 2D fit - SECOND - mF1
            obj.xCenter = round(obj.xCenter2); % should be indexes (integers)
            obj.yCenter= round(obj.yCenter2);
            obj.xUnitSize = obj.xRes2.sigmaX;
            obj.yUnitSize = obj.yRes2.sigmaY;
            obj.maxVal = obj.xRes2.Ax;
            %             [obj.ROILeft obj.ROITop obj.ROIRight obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, obj);
            [obj.ROILeft obj.ROITop obj.ROIRight obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, obj);
            %             [pic x0 y0] =  appData.data.plots{appData.consts.plotTypes.ROI}.getAnalysisPic(appData); % only ROI
            pic = pic([obj.ROITop : obj.ROIBottom] -y0+1, [obj.ROILeft : obj.ROIRight] );
            x0 = obj.ROILeft;
            y0 = obj.ROITop;
            [h w] = size(pic);
            
            firstGuess = [mean([obj.xRes2.Ax obj.yRes2.Ay]) obj.xRes2.x0 obj.yRes2.y0 ...
                obj.xRes2.sigmaX obj.yRes2.sigmaY mean([obj.xRes2.Cx obj.yRes2.Cy])];
            if ( numel(pic) > 350^2 )
                binnedPic = binning(pic, 2);
                [h w] = size(binnedPic);
                [X, Y] = meshgrid([1 : 2 : 2*w]+x0-1 , [1 : 2 : 2*h]+y0-1);% - appData.data.camera.chipStart);
                [fitRes2, obj.fval2, obj.exitflag2, obj.output2] = ...
                    fminsearch(@(p) fitGauss2D_scalar( p, X, Y, binnedPic), firstGuess, optimset('TolX',1e-4) );
            else
                [X, Y] = meshgrid([1  : w] +x0-1, [1  : h]+y0-1 );
                [fitRes2, obj.fval2, obj.exitflag2, obj.output2] = ...
                    fminsearch(@(p) fitGauss2D_scalar( p, X, Y, pic), firstGuess, optimset('TolX',1e-4) );
            end
            % set 2D params
            obj.OD2 = fitRes2(1);
            obj.x02 = fitRes2(2);
            obj.y02 = fitRes2(3);
            obj.sigmaX2 = fitRes2(4);
            obj.sigmaY2 = fitRes2(5);
            obj.C2 = fitRes2(6);
            
            obj.xCenter = round(obj.x02); % should be indexes (integers)
            obj.yCenter= round(obj.y02);
            obj.xUnitSize = obj.sigmaX2;
            obj.yUnitSize = obj.sigmaY2;
            obj.maxVal = obj.OD2;
            
            [obj.ROILeft obj.ROITop obj.ROIRight obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, obj);
            
            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
            [h w] = size(pic);
            obj.fAtomNo = 2;
            obj.atomsNo2 = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
                [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1);
            
            obj.fAtomNo = 0;
            
            obj.atomsNo=obj.atomsNo2+obj.atomsNo1;
            if (obj.y01 > obj.y02)
                %                t = obj.atomsNo1;
                %                obj.atomsNo1 = obj.atomsNo2;
                %                obj.atomsNo2 = t;
                obj.mF1 = obj.atomsNo2/obj.atomsNo;
            else
                obj.mF1 = obj.atomsNo1/obj.atomsNo;
            end;
            %            obj.mF1 = obj.atomsNo1/obj.atomsNo; %atomNo2 = mF1
            
            %set fit params
            if ( obj.OD1 > obj.OD2 )
                obj.xCenter = round(obj.x01); % should be indexes (integers)
                obj.yCenter= round(obj.y01);
                obj.xUnitSize = obj.sigmaX1;
                obj.yUnitSize = obj.sigmaY1;
                obj.maxVal = obj.OD1;
            else
                obj.xCenter = round(obj.x02); % should be indexes (integers)
                obj.yCenter= round(obj.y02);
                obj.xUnitSize = obj.sigmaX2;
                obj.yUnitSize = obj.sigmaY2;
                obj.maxVal = obj.OD2;
            end
            
            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getAnalysisPic(appData);
            [h w] = size(pic);
            
            % calc ROI size (use ROIUnits.m) - MUST be after fit
            [obj.ROILeft obj.ROITop obj.ROIRight obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, obj);
            
            %            obj.atomsNo = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
            %                [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1);
            %            obj.atomsNo=obj.atomsNo2+obj.atomsNo1;
            
            % last
            appData.data.fits{appData.consts.fitTypes.SG} = obj;
            appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);
        end
        
        function normalizedROI = getNormalizedROI(obj, pic, x, y) % return normalized ROI (to the fitting constant)
            switch obj.fAtomNo
                case 0
                    normalizedROI = pic(y, x) - mean([obj.C1 obj.C2]);
                case 1
                    normalizedROI = pic(y, x) - obj.C1;
                case 2
                    normalizedROI = pic(y, x) - obj.C2;
            end
            %            normalizedROI = pic(y, x);
        end
        
        function normalizedROI = getTheoreticalROI(obj, pic, x, y) %#ok<INUSL>
            [X, Y] = meshgrid(x, y);
            %            normalizedROI = mean([obj.xRes.Ax obj.yRes.Ay]) * exp( -0.5 * ((X-obj.xRes.x0).^2 / obj.xRes.sigmaX^2 + (Y-obj.yRes.y0).^2 / obj.yRes.sigmaY^2) );
            switch obj.fAtomNo
                case 0
                    %                    normalizedROI = pic(y, x) - mean([obj.C1 obj.C2]);
                    normalizedROI = obj.OD1* exp( -0.5 * ((X-obj.x01 ).^2 / obj.sigmaX1^2 + (Y-obj.y01).^2 / obj.sigmaY1^2) ) + ...
                        obj.OD2* exp( -0.5 * ((X-obj.x02 ).^2 / obj.sigmaX2^2 + (Y-obj.y02).^2 / obj.sigmaY2^2) );
                case 1
                    %                    normalizedROI = pic(y, x) - obj.C1;
                    normalizedROI = obj.OD1* exp( -0.5 * ((X-obj.x01 ).^2 / obj.sigmaX1^2 + (Y-obj.y01).^2 / obj.sigmaY1^2) );
                case 2
                    %                    normalizedROI = pic(y, x) - obj.C2;
                    normalizedROI = obj.OD2* exp( -0.5 * ((X-obj.x02 ).^2 / obj.sigmaX2^2 + (Y-obj.y02).^2 / obj.sigmaY2^2) );
            end
            %             normalizedROI = pic(y, x);
        end
        
        function normalizedPic = normalizePic(obj, pic)

            if ( (size(pic,2)==length(obj.xData)) && (size(pic,1)==length(obj.yData)) ) %if it is a full absorption picture (not ROI), set different scaling for each area (ROI and rest of picture)
                normalizedPic = (pic - mean([obj.C1 obj.C2])) / min([ obj.OD1 obj.OD2]);
                normalizedPic((obj.ROITop:obj.ROIBottom)-obj.yStart+1, obj.ROILeft:obj.ROIRight) = (pic((obj.ROITop:obj.ROIBottom)-obj.yStart+1, obj.ROILeft:obj.ROIRight)- mean([obj.C1 obj.C2])) / obj.maxVal;
            else
                normalizedPic = (pic - mean([obj.C1 obj.C2])) / obj.maxVal;
            end
%             if obj.OD1 < obj.OD2
        end
        
        function [xFit yFit] = getXYFitVectors(obj, x, y)
            xFit1 = obj.OD1 * exp( -0.5 * (x-obj.x01).^2 / obj.sigmaX1^2 ) + obj.C2;            
            xFit2 = obj.OD2 * exp( -0.5 * (x-obj.x02).^2 / obj.sigmaX2^2 ) + obj.C2;
            yFit1 = obj.OD1 * exp( -0.5 * (y-obj.y01).^2 / obj.sigmaY1^2 ) + obj.C1;% + obj.OD2 * exp( -0.5 * (y-obj.y02).^2 / obj.sigmaY2^2 ) ;% + obj.yRes.Cy;
            yFit2 = obj.OD2 * exp( -0.5 * (y-obj.y02).^2 / obj.sigmaY2^2 ) + obj.C2;%+ obj.OD2 * exp( -0.5 * (y-obj.y02).^2 / obj.sigmaY2^2 )% + obj.yRes.Cy;
            
           xFit = [xFit2; xFit1];
           yFit = [yFit2; yFit1];
        end
        
        function  plotFitResults(obj, appData)  % plots the text
            %            [pic x0 y0] = appData.data.plots{appData.data.plotType}.getPic();
            chipStart = appData.data.camera.chipStart;
            %            text( 10, 190, ['Atoms Num: ' num2str(obj.atomsNo/1e6) '*10^6'], 'fontSize', 20);
            %            text( 50, 140, {['x_0 = ' num2str((obj.xCenter-x0+1) * appData.data.camera.xPixSz * 1000) ' mm'], ...
            %                ['y_0 = ' num2str((obj.yCenter-y0+1) * appData.data.camera.yPixSz * 1000) ' mm']}, ...
            %                'fontsize', 12);
            text( 10, 190, ['Atoms Num: (' num2str(obj.atomsNo1/1e3) ' + ' num2str(obj.atomsNo2/1e3) ')*10^3'], 'fontSize', 16);
            text( 10, 145, {'fit function = A1 * {\ite}^{-(x-x_{01})^2 / 2\sigma_{x1}^2 - (y-y_{01})^2 / 2\sigma_{y1}^2} + C1' ...
                '                  + A2 * {\ite}^{-(x-x_{02})^2 / 2\sigma_{x2}^2 - (y-y_{02})^2 / 2\sigma_{y2}^2} + C2'}, 'fontsize', 10);
            %            text( 10, 155, 'fit function = A_{y1} * {\ite}^{-(y-y_{01})^2
            %            / 2\sigma_{y1}^2}+A_{y2} * {\ite}^{-(y-y_{02})^2 / 2\sigma_{y2}^2}', 'fontsize', 12);
            %            text( 50, 80, {['A_x = ' num2str(obj.xRes.Ax)], ...
            %                ['x_0 = ' num2str((obj.xRes.x0) * appData.data.camera.xPixSz * 1000) ' mm'], ...
            %                ['\sigma_x = ' num2str(obj.xRes.sigmaX * appData.data.camera.xPixSz * 1000) ' mm'], ...
            %                ['C_x = ' num2str(obj.xRes.Cx)], ...
            %                ['R^2 = ' num2str(obj.xGof.rsquare)]}, ...
            %                'fontsize', 12);
            text( 50, 60, {['A_{y1} = ' num2str(obj.OD1)], ...
                ['x_{01} = ' num2str(obj.x01 * appData.data.camera.xPixSz * 1000) ' mm'], ...
                ['\sigma_{x1} = ' num2str(obj.sigmaX1 * appData.data.camera.xPixSz * 1000) ' mm'] ...
                ['y_{01} = ' num2str((obj.y01-chipStart) * appData.data.camera.yPixSz * 1000) ' mm'], ...
                ['\sigma_{y1} = ' num2str(obj.sigmaY1 * appData.data.camera.yPixSz * 1000) ' mm']}, ...
                'fontsize', 10);
            text( 200, 60, {['A_{y2} = ' num2str(obj.OD2)], ...
                ['x_{02} = ' num2str(obj.x02 * appData.data.camera.xPixSz * 1000) ' mm'], ...
                ['\sigma_{x2} = ' num2str(obj.sigmaX2 * appData.data.camera.xPixSz * 1000) ' mm'] ...
                ['y_{02} = ' num2str((obj.y02-chipStart) * appData.data.camera.yPixSz * 1000) ' mm'], ...
                ['\sigma_{y2} = ' num2str(obj.sigmaY2 * appData.data.camera.yPixSz * 1000) ' mm']}, ...
                'fontsize', 10);
            text( 350, 60, {['C_{1} = ' num2str(obj.C1)] ...
                ['C_{2} = ' num2str(obj.C2)] ...
                ' ' ...
                ['m_F1 = ' num2str(100*obj.mF1) '%'] ...
                '           (of total)'}, ...
                'fontsize', 10);
        end
        
    end
end


function ret = fitGauss2D_matrix( p, X, Y, g ) %#ok<DEFNU>
% amp = p(1);
% cx = p(2);
% cy = p(3);
% wx = p(4);
% wy = p(5);
% C = p(6)
ret = p(1)*( exp( -0.5 * (X - p(2)).^2./p(4).^2 - 0.5 * (Y - p(3)).^2./p(5).^2 ) ) + p(6) - g;
end

function ret = fitGauss2D_scalar( p, X, Y, g )
% amp = p(1);
% cx = p(2);
% cy = p(3);
% wx = p(4);
% wy = p(5);
% C = p(6)
ret = p(1)*( exp( -0.5 * (X - p(2)).^2./p(4).^2 - 0.5 * (Y - p(3)).^2./p(5).^2 ) ) + p(6) - g;
ret = sum(sum(ret.^2));
end
