classdef Fit2DFermi < FitTypes
    %FIT1DGAUSSIAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties ( Constant = true )
        ID = 'Fit2DFermi';
    end
    properties (SetAccess = private )
        xRes1DG = []; %1DGaussian x fit results
        yRes1DG = []; %1DGaussian y fit results
        res = [];
        gof = [];
%         output = []; %eliminated, for some reason it takes 4MBytes size
        conf = [];
        OD = -1;
        x0 = -1;
        y0 = -1;
        Rx = -1;
        Ry = -1;
        q = -1;
        C = -1;
        ToverTF = [];
        ToverTFError = [];
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
            firstGuess = [0.5*(fitObj.xRes.Ax+fitObj.yRes.Ay) 0 obj.xRes1DG.x0, obj.yRes1DG.y0 2*fitObj.xRes.sigmaX ...
                2*fitObj.yRes.sigmaY 0.5*(fitObj.xRes.Cx+fitObj.yRes.Cy)];
            lower = [0 -10 x0 y0 0.5*fitObj.xRes.sigmaX 0.5*fitObj.yRes.sigmaY -10]; %#ok<PROPLC>
            upper = [10 10 x0+w y0+h 4*fitObj.xRes.sigmaX 4*fitObj.yRes.sigmaY 10 ];%#ok<PROPLC>
            binW = appData.options.avgWidth;
            if binW==0
                binW = 1;
            end
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', firstGuess, 'Lower', lower, 'Upper', upper, 'TolFun', 1e-010);
            f = fittype('Fermi2D( n0, q, x - x0, Rx, y - y0, Ry ) + C', 'coefficients', {'n0', 'q', 'x0', 'y0', 'Rx', 'Ry', 'C'},...
                'independent', {'x','y'}, 'dependent', 'z', 'options', s);
            if ( numel(pic) > 450^2 )
                binnedPic = binning(pic, binW);
                [h, w] = size(binnedPic);
                X = (1 : binW : binW*w) +x0-1; %#ok<PROPLC>
                Y = (1 : binW : binW*h) +y0-1; %#ok<PROPLC>
                [X, Y, binnedPic] = prepareSurfaceData(X, Y, binnedPic);
                [obj.res, obj.gof, output] = fit([X, Y], binnedPic, f);
            else
                X = (1 : w) +x0-1; %#ok<PROPLC>
                Y = (1 : h) +y0-1; %#ok<PROPLC>
                [X, Y, binnedPic] = prepareSurfaceData(X, Y, pic);
                [obj.res, obj.gof, output] = fit([X, Y], binnedPic, f);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % set 2D params
            obj.conf = confint(obj.res);
            obj.conf = (obj.conf(2,:)-obj.conf(1,:))/2;
           
%             obj.x0 = obj.xRes1DG.x0; % center taken from from 1DGaussian
%             obj.y0 = obj.yRes1DG.y0; % center taken from from 1DGaussian
            
            obj.x0 = obj.res.x0; 
            obj.y0 = obj.res.y0;
            
            obj.OD = obj.res.n0;
            obj.q = obj.res.q;
            obj.Rx = obj.res.Rx;
            obj.Ry = obj.res.Ry;
            obj.C = obj.res.C;
            obj.ToverTF = real((-6*polylog(3,-exp(obj.q)))^(-1/3)); %polylog sometimes gives a very small imaginary part to the result
            obj.ToverTFError = real((-6*polylog(3,-exp(obj.q - obj.conf(2))))^(-1/3) - (-6*polylog(3,-exp(obj.q + obj.conf(2))))^(-1/3))./2;
            
            %%% test
%             picFit = (Fermi2D( obj.OD, obj.q, X - obj.x0, obj.Rx, Y - obj.y0, obj.Ry ) + obj.C);
%             picFit = reshape(picFit,size(binnedPic));
%             figure;
%             imagesc(picFit)
            
            %%%%%%%%%%% plot the residuals:
% %             disp(['Mean 2D Fermi-Dirac polylog fit resiuals RMS = ' num2str(sqrt(mean(output.residuals(:).^2))) ])
%             figure;
% %             imagesc(LowPassFilter(reshape(output.residuals,size(pic)),30,9) )
%             imagesc( reshape(output.residuals,size(pic)) )
%             colorbar
            %%%%%%%%%%%%%
            
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
           normalizedROI = pic(y, x) - obj.C;
       end
       
       function normalizedROI = getTheoreticalROI(obj, pic, x, y)
           [X, Y] = meshgrid(x, y);
%            normalizedROI = mean([obj.xRes1DG.Ax obj.yRes1DG.Ay]) * exp( -0.5 * ((X-obj.xRes1DG.x0).^2 / obj.xRes1DG.sigmaX^2 + (Y-obj.yRes1DG.y0).^2 / obj.yRes1DG.sigmaY^2) );
       end
       
       function normalizedPic = normalizePic(obj, pic)
           normalizedPic = (pic-obj.C)/obj.maxVal;
       end
        
       function [xFit, yFit] = getXYFitVectors(obj, x, y)
           %             Fermi2D( n0, q, x, Rx, y, Ry ) + c
           xFit = (Fermi2D( obj.OD, obj.q, x - obj.x0, obj.Rx, 0*obj.y0, obj.Ry ) + obj.C)';
           yFit = (Fermi2D( obj.OD, obj.q, 0*obj.x0 , obj.Rx, y - obj.y0, obj.Ry ) + obj.C)';
           
           %%%%%%%%%%%%%%%%%%%
% %            calculate the expected profile, from atom number constraint and wings temperature fit
%            kB = 1.3806504e-23; %J/K
% %            hbar = 1.054571628e-34;
%            m = 6/6.022e26;
%            pixSz = 6.45e-6;
%            sigmaAbs = 3*670.977338e-9^2/2/pi;
%            TOF = 0.8e-3;
%            out = calcXODTDensity(1*obj.atomsNo, 0.5, TOF*1e3, obj.xRes1DG.sigmaX * pixSz * 1000, obj.yRes1DG.sigmaY * pixSz * 1000);
%            % Temperature from wings of the FD fit:
%            Tx = (obj.Rx * pixSz./TOF).^2 *m ./ (2*kB) .* exp(obj.q)./( (1+exp(obj.q)).*log(1+exp(obj.q)) ) ;
%            Ty = (obj.Ry * pixSz./TOF).^2 *m ./ (2*kB) .* exp(obj.q)./( (1+exp(obj.q)).*log(1+exp(obj.q)) ) ;
%            % Temperature from wings of the FD fit, also include finite TOF effect:
% %            Tx = m ./ (2*kB) * out.omega_r^2 * (obj.Rx * pixSz)^2 ./ ( 1 + (out.omega_r*TOF)^2 ) .* exp(obj.q)./( (1+exp(obj.q)).*log(1+exp(obj.q)) ) ;
% %            Ty = m ./ (2*kB) * out.omega_z^2 * (obj.Ry * pixSz)^2 ./ ( 1 + (out.omega_z*TOF)^2 ) .* exp(obj.q)./( (1+exp(obj.q)).*log(1+exp(obj.q)) ) ;
%            T = mean([Tx Tx Ty]);
%            %             lambdaDB = sqrt( 2*pi*hbar^2 /(m*kB*out.T) );
%            f = @(x) (1+x)./x.*log(1+x);
%            %             fq = (1+exp(obj.q))./exp(obj.q).*log(1+exp(obj.q));
%            n2D0 = obj.atomsNo .* polylog(2,-exp(obj.q)) .* f(exp(obj.q))./( pi* obj.Rx * pixSz * obj.Ry * pixSz * polylog(3,-exp(obj.q)) ) ; %2D peak column density [m^2]
%            OD0 = sigmaAbs * n2D0; % sanity check - this returns exactly the experimental OD value, as it should
%            
%            % expected value for peak OD from Temperature, Fermi energy, nmber of atoms, and trap frequencies:
%            mu = kB*out.Tf*( 1 - pi^2/3*(T./out.Tf)^2 ); %Chemical potential, Sommerfeld expansion
%            qTh = mu./(kB*T); %theoretical value of the logarighmic fugacity. Using T here from the fit seems like bootstrapping, but the result for the temperature agrees with the Gaussian wings fit (should double check);
%            sigmax = sqrt(2*kB*Tx./(m*out.omega_r.^2)) * sqrt( 1 + out.omega_r^2*TOF^2);
%            sigmay = sqrt(2*kB*Ty./(m*out.omega_z.^2)) * sqrt( 1 + out.omega_z^2*TOF^2);
%            Rx = sigmax * sqrt(f(exp(qTh))); %#ok<PROPLC> % Cloud size in x [m]
%            Ry = sigmay * sqrt(f(exp(qTh))); %#ok<PROPLC> % Cloud size in y [m]
%            OD0Theretical = sigmaAbs * obj.atomsNo * polylog(2,-exp(qTh))* f(exp(qTh))./( pi* Rx  * Ry  * polylog(3,-exp(qTh)) ) ; %#ok<PROPLC> %2D peak column density [m^2], calcualted using the calibrated atoms number, and the expected Fermi temprature from the number of atoms
%            OD0Theretical = real(OD0Theretical);
%            
%            %             n2D0 = -m * (kB*T).^2./(2*pi*out.omega_r * hbar^3) ./ sqrt( 1 + out.omega_r^2*TOF^2) ; %2D peak column density [m^2]
%            %             argX = -exp(mu./(kB*T) - ((x- obj.x0)*pixSz).^2./2./(out.sigma_r*sqrt(1 + out.omega_r^2.*TOF.^2))^2);
%            %             xFit2 = OD0 .* real(polylog(2, argX ))./polylog(2,-exp(obj.q)) + obj.C; %polylog returns imaginary values for some values, e.g. polylog(2,-1.0604)
%            
%            % calcualte vectors using measured Rx,Ry
%            %             xFit2 = (Fermi2D( OD0Theretical, obj.q, x - obj.x0, obj.Rx, 0*obj.y0, obj.Ry ))' + obj.C;
%            %             yFit2 = (Fermi2D( OD0Theretical, obj.q, 0*obj.x0 , obj.Rx, y - obj.y0, obj.Ry ))' + obj.C;
%            
%            % calculate vectors using calculated Rx,Ry
%            xFit2 = (Fermi2D( OD0Theretical, qTh, x - obj.x0, Rx./pixSz, 0*obj.y0, Ry./pixSz ))' + obj.C; %#ok<PROPLC>
%            yFit2 = (Fermi2D( OD0Theretical, qTh, 0*obj.x0 , Rx./pixSz, y - obj.y0, Ry./pixSz ))' + obj.C; %#ok<PROPLC>
%            
%            xFit = [xFit; xFit2];
%            yFit = [yFit; yFit2];
           
           
           
           %%%%%%%%%%%%%
% %            compare result without the f(x) function, and using simply sigma_x,y
% %            Conclusion - gives the same result, as expected. So f(e^q) is just for convenience of scaling R_i, and sigma_i is what really appears in the argument
%            arg1 = -exp(qTh - ( (x- obj.x0).^2./(sigmax./pixSz)^2 )); %+ (y- obj.y0).^2./(sy./pixSz)^2
%            arg2 = -exp(qTh);
%            xFit3 = OD0Theretical *real(polylog(2, arg1 ))./real(polylog(2, arg2 )) + obj.C;
%            
%            arg1 = -exp(qTh - ( (y- obj.y0).^2./(sigmay./pixSz)^2 )); %+ (y- obj.y0).^2./(sy./pixSz)^2
%            arg2 = -exp(qTh);
%            yFit3 = OD0Theretical *real(polylog(2, arg1 ))./real(polylog(2, arg2 )) + obj.C;
%            
%            xFit = [xFit; xFit2; xFit3];
%            yFit = [yFit; yFit2; yFit3];
%            %Bottom line xFit2 is narrower than the fit, which sort of makes sense(?). Should calcualte the temperature from the Gaussian wings fit, and not the FD fit!
        end
        
        function  plotFitResults(obj, appData)  % plots the text
            chipStart = appData.data.camera.chipStart;
            cla(appData.ui.pFitResults)
            text( 0, 1, ['Atoms Num: ' addCommas(obj.atomsNo)], 'fontSize', 20, 'Parent', appData.ui.pFitResults);
            text( 0, 0.8, 'fit function = Fermi2D( n0, q, x, Rx, y, Ry ) + c', 'fontsize', 12, 'Parent', appData.ui.pFitResults);
            text( 0, 0.4, {['n_0 = ' num2str(obj.OD)], ... %                ['A_G = ' num2str(obj.rRes1DG.A)], ...
                ['q = ' num2str(obj.q) ' +/- ' num2str(obj.conf(2), '%.4f')],... 
                ['T/T_F = ' num2str(obj.ToverTF) '+/-' num2str(obj.ToverTFError)], ... %' +/- ' num2str(obj.ToverTFError)
                ['(x_0, y_0) = (' num2str( obj.x0 * appData.data.camera.xPixSz * 1000) '; '...
                num2str((obj.y0-chipStart) * appData.data.camera.yPixSz * 1000) ') mm'], ...
                ['(R_x, R_y) = ' num2str(obj.Rx * appData.data.camera.xPixSz * 1e6, '%.2f') ';   '...
                num2str(obj.Ry * appData.data.camera.yPixSz * 1e6, '%.2f') ' um'],...
                ['C = ' num2str(obj.C) '+/-' num2str(obj.conf(7))]},...
                'fontsize', 12, 'Parent', appData.ui.pFitResults);
        end
        
    end
end

