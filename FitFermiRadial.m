classdef FitFermiRadial < FitTypes
    %FIT1DGAUSSIAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties ( Constant = true )
        ID = 'FitFermiRadial';
    end
    properties (SetAccess = private )
        xRes1DG = []; %1DGaussian x fit results
        yRes1DG = []; %1DGaussian y fit results
        rRes1DG = []; %1DGaussian radial fit results, on the azimuthaly-averaged data
        res = [];
        gof = [];
        output = [];
        conf = [];
        ToverTF = [];
        ToverTFError = [];
        FDTheoryToFitRatio = [];
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
            
            [pic, x0, y0] = appData.data.plots{appData.consts.plotTypes.ROI}.getAnalysisPic(appData); % only ROI
%             pic = pic(round(fitObj.yCenter-y0-3*fitObj.yRes.sigmaY) : round(fitObj.yCenter-y0+3*fitObj.yRes.sigmaY),...
%                 round(fitObj.xCenter-x0-3*fitObj.xRes.sigmaX) : round(fitObj.xCenter-x0+3*fitObj.xRes.sigmaX)); %take only 3 sigma around the center, to reduce amount of data
            [h, w] = size(pic);
            [X, Y] = meshgrid( (1: w) +x0-fitObj.xCenter-1, (1 : h) +y0-fitObj.yCenter-1); %need to check if the grid is centered around the cloud center
            
            %1D guassian fit for the azimuthaly averaged data
            [r, GAvg, GSem] = AzimuthalAverage( pic, X, Y ); %calculate azimuthal average. r is the radius, Gavg is the azimuthaly averaged OD
            startPoint = [0.5*(fitObj.xRes.Ax+fitObj.yRes.Ay) abs(0.5*(fitObj.xRes.sigmaX + fitObj.yRes.sigmaY)) 0.5*(fitObj.xRes.Cx + fitObj.yRes.Cy)]; %0.5*(fitObj.xRes.Cx + fitObj.yRes.Cy)
            lower = [0.2*startPoint(1) 0.2*startPoint(2) -2];
            upper = [2*startPoint(1) 2*startPoint(2) 2];
            weights = 1./(GSem);
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint, 'Lower', lower, 'Upper', upper, 'Weights', weights);
            f = fittype('A*exp(-r^2/(2*sigma^2))+C', 'coefficients', {'A', 'sigma', 'C'}, 'independent', 'r', 'dependent', 'y', 'options', s);
            [obj.rRes1DG, ~, output1DG]= fit(r', GAvg', f);
            
            %%% do the Gaussian fit again, this time only to the wings:            
            weightsGauss = weights;
            lower = [0.2*startPoint(1) 0.2*startPoint(2) -2];
            upper = [5*startPoint(1) 2*startPoint(2) 2];
            weightsGauss( abs(r)< 1.75*obj.rRes1DG.sigma) = 0;
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint, 'Lower', lower, 'Upper', upper, 'Weights', weightsGauss);
            f = fittype('A*exp(-r^2/(2*sigma^2))+C', 'coefficients', {'A', 'sigma', 'C'}, 'independent', 'r', 'dependent', 'y', 'options', s);
            [rRes1DGWings, ~, ~ ]= fit(r', GAvg', f);
            
            % process data after Gaussian fit: remove anyhting after 4 sigma (dosen't affect the fit)
%             GAvgCenter = GAvg(abs(r)<4*obj.rRes1DG.sigma); %reduce data to 4 sigma of Gaussian from the center
%             rCenter = r(abs(r)<4*obj.rRes1DG.sigma);
%             weights = weights(abs(r)<4*obj.rRes1DG.sigma);
%             remove center of data (noisy) - This part has been omitted, since we noe use weighting of the data points
%             GAvgCenter = GAvgCenter(0.5*obj.rRes1DG.sigma < abs(rCenter)); %cut data to 0.5 sigma of center, usually bad
%             rCenter = rCenter(0.5*obj.rRes1DG.sigma < abs(rCenter));
            
            GAvgCenter = GAvg;
            rCenter = r;
            %remove only the center data:
%             GAvgCenter = GAvg(0.5*obj.rRes1DG.sigma < abs(r)); %cut data to 0.5 sigma of center, usually bad
%             rCenter = r(0.5*obj.rRes1DG.sigma < abs(r));
            
              % 1D fit function - PolyLog 5/2, seems to be wrong, as integrating over \phi dosen't change PolyLog 2-> PolyLog 5/2
%             startPoint = [obj.rRes1DG.A -5 obj.rRes1DG.sigma*1.8 ];
%             lower = [0.2*startPoint(1) -60 0.05*startPoint(3)];
%             upper = [2*startPoint(1) 6 5*startPoint(3)];
%             s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint, 'Lower', lower, 'Upper', upper,...
%                 'MaxFunEvals',1000, 'MaxIter',1000, 'Weights', weights);
%             f = fittype('Fermi1D( n0, q, x, R ) + C', 'coefficients', {'n0', 'q', 'R'}, 'independent', 'x', 'dependent', 'y', 'options', s, 'problem', {'C'}); %
%             [obj.res, obj.gof, obj.output] = fit(rCenter', GAvgCenter', f, 'problem', { obj.rRes1DG.C }); %constant not used as a fit parameter - not necessary for the polylog fit.
            
            % Since radial averaging dosen't deacrese the dimension (assuming Rx=Ry), we use a 2D fit function PolyLog 2, using x=r, y=0, Ry=1 (to avoid division by 0)
            % using C as a fit parameter doesn't seems to affect the error in q
            startPoint = [obj.rRes1DG.A 0 obj.rRes1DG.sigma*1.8 obj.rRes1DG.C];
            lower = [0.2*startPoint(1) -60 0.05*startPoint(3) -1];
            upper = [2*startPoint(1) 6 5*startPoint(3) 1];
            s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint, 'Lower', lower, 'Upper', upper,...
                'MaxFunEvals',1000, 'MaxIter',1000, 'Weights', weights);
            f = fittype('Fermi2D( n0, q, r , R, 0, 1 ) + C', 'coefficients', {'n0', 'q', 'R', 'C'},...
                'independent', {'r'}, 'dependent', 'y', 'options', s);
            [obj.res, obj.gof, obj.output] = fit(rCenter', GAvgCenter', f); %constant not used as a fit parameter - not necessary for the polylog fit.
            
%             disp([num2str(obj.res.n0/startPoint(1)) ' ; ' num2str(obj.res.q/startPoint(2)) ' ; ' num2str(obj.res.R/startPoint(3))])
            obj.conf = confint(obj.res);
            obj.conf = (obj.conf(2,:)-obj.conf(1,:))/2;
            
            obj.ToverTF = real((-6*polylog(3,-exp(obj.res.q)))^(-1/3));
            obj.ToverTFError = real((-6*polylog(3,-exp(obj.res.q - obj.conf(2))))^(-1/3) - (-6*polylog(3,-exp(obj.res.q + obj.conf(2))))^(-1/3))./2;           
            
%             % display temperature, single arm ODT values
%             omega_r = (2*pi*25e3);
%             omega_z = (2*pi*600);
%             f = @(x) (1+x)./x*log(1+x);
%             TOF = 1.3e-3; %TOF in sec
%             R = obj.res.R*appData.data.camera.yPixSz; % R value in meters
%             MLi6 = 9.9883414e-27; 
%             Tr = 0.5*MLi6*(omega_r)^2*R^2./(1+omega_r^2*TOF^2)./f(exp(obj.res.q))/1.38e-23*1e6;
%             Tz = 0.5*MLi6*(omega_z)^2*R^2./(1+omega_z^2*TOF^2)./f(exp(obj.res.q))/1.38e-23*1e6;
%             disp(['(T_r, T_z) = (' num2str(Tr) ', ' num2str(Tz) ') uK'])
            
            %set fit params
            obj.xCenter = round(fitObj.xRes.x0); % should be indexes (integers)
            obj.yCenter= round(fitObj.yRes.y0);
            obj.xUnitSize = fitObj.xRes.sigmaX;
            obj.yUnitSize = fitObj.yRes.sigmaY;
            obj.maxVal = obj.res.n0;
            
            % calc ROI size (use ROIUnits.m) - MUST be after fit
            [obj.ROILeft obj.ROITop obj.ROIRight obj.ROIBottom] = appData.data.ROITypes{appData.data.ROIUnits}.getROICoords(appData, fitObj);
            
%             obj.atomsNo = fitObj.atomsNo;
            obj.atomsNo = appData.options.calcs{appData.options.calcAtomsNo}.calcAtomsNo(appData, obj, pic, ...
                [obj.ROILeft : obj.ROIRight] - x0+1, [obj.ROITop : obj.ROIBottom] - y0+1);
            
            obj.xData = fitObj.xData;
%             obj.xStart = fitObj.x0;
            obj.yData = fitObj.yData;
%             obj.yStart = fitObj.y0;

            %%%%%%%%%%%%%% plot fit results in a seperate figure
            fig = figure('Position', [463 50 1034 679]);
            ax = subplot(2,1,1);
            plot(r, GAvg, '.', 'markerSize', 10); hold on;
            plot(rCenter, obj.res(rCenter),'-r')
%             plot(rCenter(rCenter>0), obj.res(rCenter(rCenter>0)),'-r') %right  part of fit
%             plot(rCenter(rCenter<0), obj.res(rCenter(rCenter<0)),'-r','HandleVisibility','off') %left part of fit
            plot(r, obj.rRes1DG(r),'--black')
            plot(r, rRes1DGWings(r),':blue')
            xlabel('Radius [pixels]');
            ylabel('Optical Density');
            legend('Data', 'PolyLog fit', 'Gaussian fit', 'Gaussian wings fit')
            xlim(round(5*obj.rRes1DG.sigma)*[-1 +1])
            
            subplot(2,1,2);
            plot(rCenter, obj.output.residuals, '-r'); hold on;
            plot(r, output1DG.residuals, '--black')
            legend( 'PolyLog fit residuals', 'Gaussian fit residuals')
            grid on;
            sp = get(fig, 'Children');
            set(sp(2), 'Position', [0.078 0.0722 0.9 0.25]);
            set(sp(4), 'Position', [0.078 0.4359 0.9 0.54]);
            xlim(round(5*obj.rRes1DG.sigma)*[-1 +1])
            
            %%%%%%%%%%%%%%%%% % calculate the expected profile, from atom number constraint and wings temperature fit
            % this is done in order to compare the p-wave fit fugacity q to the 'theoretical' value of the fugacity qTh (i.e. calculated from trap frequencies and number of atoms).
%             kB = 1.3806504e-23; %J/K
%             m = 6/6.022e26;
%             pixSz = 6.45e-6;
%             sigmaAbs = 3*670.977338e-9^2/2/pi;
%             TOF = 0.3e-3;
%             out = calcXODTDensity(obj.atomsNo*1.8, 3.5, TOF*1e3, obj.xRes1DG.sigmaX * pixSz * 1000, obj.yRes1DG.sigmaY * pixSz * 1000);
%             
%             % Since we are assuming azimuthal symmetry, we calcualte only a single temperature
%             Twings = (rRes1DGWings.sigma * pixSz./TOF).^2 *m ./ kB; % Temperature from wings of the Gaussian fit
%             TFD = (obj.res.R * pixSz./TOF).^2 *m ./ (2*kB) .* exp(obj.res.q)./( (1+exp(obj.res.q)).*log(1+exp(obj.res.q)) ) ; % Temperature from wings of the FD fit, this gives different results - seems like the error in determining the temprature is too big.
%             % calculate FD temperature error:
% %             bounds = calcBounds(obj.res); % 'coefficients', {'n0', 'q', 'R', 'C'}
% %             [TFD, TFDErr] = propError('(R*6.45e-6/T)^2*(6/6.022e26)/(2*1.3806504e-23) * exp(q)./( (1+exp(q)).*log(1+exp(q)) )', {'R', 'T', 'q'}, {obj.res.R, TOF, obj.res.q}, {bounds(3), 0, bounds(2)});
%             T = Twings;
%             f = @(x) (1+x)./x.*log(1+x);
%             n2D0 = obj.atomsNo .* polylog(2,-exp(obj.res.q)) .* f(exp(obj.res.q))./( pi* obj.res.R.^2 * pixSz.^2 * polylog(3,-exp(obj.res.q)) ) ; %2D peak column density [m^2]
%             OD0 = sigmaAbs * real(n2D0); % sanity check - this returns exactly the experimental OD value, as it should. The reduced absorption cross section (by factor ~1.8) can go also be accounted for by using the uncalibrated number of atoms, so it's the same
%             
%             %%% expected value for peak OD from Temperature, Fermi energy, nmber of atoms, and trap frequencies:
%             % fugacity from Sommerfeld expansion - not accurate enough at high temperatures, 20% error at T/TF=0.3
%             mu = kB*out.Tf*( 1 - pi^2/3*(T./out.Tf)^2 ); %Chemical potential, Sommerfeld expansion
%             qTh = mu./(kB*T); %theoretical value of the logarighmic fugacity. Using T here from the fit seems like bootstrapping, but the result for the temperature agrees with the Gaussian wings fit (should double check);
%             % Fugacity from polylog (Eq. 24 in Varenna notes):
%             syms z
%             fugacity = double(vpasolve(polylog(3,-z)+1/6/(T./out.Tf)^3));
%             qTh = log(fugacity);
%             
%             sigma = sqrt(2*kB*T./m) * TOF; % 'Gaussian' cloud size (i.e. only from temperature)
%             R = sigma * sqrt(f(exp(qTh))); % Cloud size in x [m]
%             OD0Theretical = sigmaAbs * obj.atomsNo * polylog(2,-exp(qTh))* f(exp(qTh))./( pi* R^2  * polylog(3,-exp(qTh)) ) ; %2D peak column density [m^2], calcualted using the calibrated atoms number, and the expected Fermi temprature from the number of atoms
%             OD0Theretical = real(OD0Theretical);
%             
%             % Save ratio for plotting:
%             obj.FDTheoryToFitRatio = OD0Theretical./OD0;
%             
%             % calcualte vectors using measured Rx,Ry
%             %             xFit2 = (Fermi2D( OD0Theretical, obj.q, x - obj.x0, obj.Rx, 0*obj.y0, obj.Ry ))' + obj.C;
%             
%             % Plot theoretical line
%             GTheory = (Fermi2D( OD0Theretical, qTh, r, R./pixSz, 0, 1 ))' + obj.res.C;
%             axes(ax(1));
%             plot(r, GTheory, ':Green', 'LineWidth', 2);
%             legend('Data', 'PolyLog fit', 'Gaussian fit', 'Gaussian wings fit', 'PolyLog with calculated fugacity')
            
            %%%%%%%%%%%%%%%%%
            
            % last
            appData.data.fits{appData.consts.fitTypes.FermiRadial} = obj;
            appData = appData.data.plots{appData.consts.plotTypes.ROI}.setPic(appData, pic);
            
            
        end
        
       function normalizedROI = getNormalizedROI(obj, pic, x, y) % return normalized ROI (to the fitting constant)
           normalizedROI = pic(y, x) - obj.res.C;
       end
       
       function normalizedROI = getTheoreticalROI(obj, pic, x, y)
           [X, Y] = meshgrid(x, y);
           normalizedROI = mean([obj.xRes1DG.Ax obj.yRes1DG.Ay]) * exp( -0.5 * ((X-obj.xRes1DG.x0).^2 / obj.xRes1DG.sigmaX^2 + (Y-obj.yRes1DG.y0).^2 / obj.yRes1DG.sigmaY^2) );
       end
       
       function normalizedPic = normalizePic(obj, pic)
           normalizedPic = (pic - obj.res.C)/obj.maxVal;
       end
        
        function [xFit, yFit] = getXYFitVectors(obj, x, y)
            xFit = obj.xRes1DG.Ax * exp( -0.5 * (x-obj.xRes1DG.x0).^2 / obj.xRes1DG.sigmaX^2 ) + obj.xRes1DG.Cx;
            yFit = obj.yRes1DG.Ay * exp( -0.5 * (y-obj.yRes1DG.y0).^2 / obj.yRes1DG.sigmaY^2 ) + obj.yRes1DG.Cy;
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
            text( 0, 0.8, 'fit function = Fermi1D( n0, q, x, R ) + c', 'fontsize', 12, 'Parent', appData.ui.pFitResults);
            text( 0, 0.3, {['n_0 = ' num2str(obj.res.n0, '%.4f')], ...
                ['A_G = ' num2str(obj.rRes1DG.A, '%.4f')], ...
                ['(x_0, y_0) = (' num2str((appData.data.fits{appData.consts.fitTypes.oneDGaussian}.xRes.x0) * appData.data.camera.xPixSz * 1000) '; '...
                num2str((appData.data.fits{appData.consts.fitTypes.oneDGaussian}.yRes.y0-chipStart) * appData.data.camera.yPixSz * 1000) ') mm'], ...
                ['R = ' num2str(obj.res.R * appData.data.camera.xPixSz * 1000, '%.4f') ' +/- ' num2str(obj.conf(3) * appData.data.camera.xPixSz * 1000, '%.4f') ' mm'], ...
                ['q = ' num2str(obj.res.q, '%.4f') ' +/- ' num2str(obj.conf(2), '%.4f')], ...
                ['T/T_F = ' num2str(obj.ToverTF, '%.4f') ' +/- ' num2str(obj.ToverTFError, '%.4f') ], ...
                ['C = ' num2str(obj.res.C, '%.4f') ' +/- ' num2str(obj.conf(4), '%.4f')]}, ...
                'fontsize', 12, 'Parent', appData.ui.pFitResults);
            %['C_x = ' num2str(obj.xRes.Cx)], ...
            %['R^2 = ' num2str(obj.xGof.rsquare)]}, ...
%             text( 200, 80, {['A_y = ' num2str(obj.yRes.Ay)], ...
%                 ['y_0 = ' num2str((obj.yRes.y0-chipStart) * appData.data.camera.yPixSz * 1000) ' mm (from the chip)'], ...
%                 ['\sigma_y = ' num2str(obj.yRes.sigmaY * appData.data.camera.yPixSz * 1000) ' mm'], ...
%                 ['C_y = ' num2str(obj.yRes.Cy)], ...
%                 ['R^2 = ' num2str(obj.yGof.rsquare)]}, ...
%                 'fontsize', 12);
        end
        
    end
end


% function ret = Fermi1D( n0, q, x, R )
% % calcualte the 1D Fermi profile
% 
% persistent polylog5halfSpline
% %spline of the polylog function to speedup computation, in the range x = [-2000:10:-1000 -990:5:-100 -100:1:-10 -10:0.5:-1 -1:0.1:1];
% if isempty(polylog5halfSpline)
%     polylog5halfSpline = load('polylog5halfSpline.mat', 'polylog5halfSpline');
% end
% arg1 = -exp(q - x.^2./R^2.*((1+exp(q))./exp(q).*log(1+exp(q))) );
% arg2 = -exp(q);
% if any([arg1<-2000 ; arg1>1 ; arg2<-2000 ; arg2>1])
%    warning('PolyLog argument out of spline approximation range (-2000 < x < 1)')
%     return
% end
% 
% ret = n0*polylog5halfSpline.polylog5halfSpline( arg1 )./polylog5halfSpline.polylog5halfSpline( arg2 );
% 
% % ret = polylog(5/2, exp(q - x.^2./R^2.*f(exp(q))))./polylog(5/2, exp(q));
% % ret = n0*polylog(5/2, exp(q - x.^2./R^2.*((1+exp(q))/exp(q)*log(1+exp(q))) ))./polylog(5/2, exp(q));
% 
% end

% function [uniqueDistances, GAvg] = AzimuthalAveragev1( G, X, Y )
% % G is the picture of the Fermi cloud, with cloud center at (xc,yc)
% % D is a matrix who's values represent the distance from the center
% % [X,Y] = meshgrid( 1:size(G,2), 1:size(G,1) );
% D = sqrt( X.^2+Y.^2 );
% 
% % Analysis
% D = reshape(D,numel(D),1); %Reshape D to a column vector
% G = reshape(G,numel(G),1); %Reshape G to a column vector
% [D, IX] = sort(D); %[D,IX]=[sorted distances, sorted indices]
% G = G(IX); %Sort G according to the distance sorting
% 
% % average over values with identical distance from center
% 
% uniqueDistances = unique(D);
% GAvg = zeros(length(uniqueDistances),1);
% % halfBandWidth = 3;
% for i = 1 : length(uniqueDistances)
%     GAvg(i) = mean(G(D==uniqueDistances(i)));
% end
% 
% figure;
% plot(uniqueDistances, GAvg)
% % N = 20;
% % GAvg = nanmean(reshape([GAvg(:); nan(mod(-numel(GAvg),N),1)],N,[]));
% % uniqueDistances = nanmean(reshape([uniqueDistances(:); nan(mod(-numel(uniqueDistances),N),1)],N,[]));
% % hold on;
% % plot(uniqueDistances, GAvg,'-or')
% 
% end

% function [uniqueDistances, GAvg] = AzimuthalAveragev2( G, X, Y )
% % AzimuthalAverage - v2
% % G is the picture of the Fermi cloud, with cloud center at (xc,yc)
% % D is a matrix who's values represent the distance from the center
% % [X,Y] = meshgrid( 1:size(G,2), 1:size(G,1) );
% % v2 does not average over values with identical distance from center, so
% % the weight data won't be lost when doing a low-pass filter
% D = sqrt( X.^2+Y.^2 );
% 
% % Analysis
% D = reshape(D,numel(D),1); %Reshape D to a column vector
% G = reshape(G,numel(G),1); %Reshape G to a column vector
% [D, IX] = sort(D); %[D,IX]=[sorted distances, sorted indices]
% G = G(IX); %Sort G according to the distance sorting
% 
% %replicate data left to right (besides center-zero point)
% G = [flipud(G(2:end)); G];
% D = [flipud(-D(2:end)); D];
% % smooth the data:
% GAvg = smooth(D, G, 800,'sgolay');
% [uniqueDistances, indx] = unique(D);
% GAvg = GAvg(indx); %omit redundent data - points that have same value after smoothing
% 
% % figure;
% % % plot(D, G, 'o')
% % % hold on;
% % plot(uniqueDistances,GAvg,'or')
% % hold on;
% %% extra redution of data points by moving average
% N = 10;
% uniqueDistances = uniqueDistances(1:end-mod(length(uniqueDistances), N)); % omit last few points
% GAvg = GAvg(1:end-mod(length(GAvg), N));
% uniqueDistances = (mean(reshape(uniqueDistances,N,[])))';
% GAvg = (mean(reshape(GAvg,N,[])))';
% 
% %% plot before and after
% % figure;
% % plot(D, G, 'o')
% % hold on;
% % plot(uniqueDistances,GAvg,'oblack')
% 
% 
% end

function [r, GAvg, GSem] = AzimuthalAverage( G, X, Y )
% AzimuthalAverage - v3
% G is the picture of the Fermi cloud, with cloud center at (xc,yc)
% D is a matrix who's values represent the distance from the center
% [X,Y] = meshgrid( 1:size(G,2), 1:size(G,1) );
% v3 averages over bands with some distance from the origin, to avoid
% over-wighting area with larger r values (as they have more points).
D = sqrt( X.^2+Y.^2 );

% Analysis
D = reshape(D,numel(D),1); %Reshape D to a column vector
G = reshape(G,numel(G),1); %Reshape G to a column vector
[D, IX] = sort(D); %[D,IX]=[sorted distances, sorted indices]
G = G(IX); %Sort G according to the distance sorting

badnWidth = 3; %pixel band width to average over. Preferably between 2 to 5
maxR = round(max(D));
r = (0 : badnWidth : maxR);
GAvg = zeros(size(r));
GSem = zeros(size(r));
NPoints = zeros(size(r));

GAvg(1) = mean([G(1) ; G(r(1)<D & D<(r(2))) ]); %insert center point into the first averaging.
GSem(1) = std([G(1) ; G(r(1)<D & D<(r(2))) ]);
NPoints(1) = 1 + sum(r(1)<D & D<(r(2)));

for i = 2 : length(r)-1
	GAvg(i) = mean(G(r(i)<D & D<(r(i+1))));
    NPoints(i) = sum(r(i)<D & D<(r(i+1)));
    GSem(i) = std(G(r(i)<D & D<(r(i+1)))) ./ sqrt(NPoints(i)); %standard error of the mean for this point
    
end
r = r + badnWidth/2; %add zero radius point (which is not averaged), and shift radius to proper position
% 
% r = [0 r + badnWidth/2]; %add zero radius point (which is not averaged), and shift radius to proper position
% GAvg = [G(1) GAvg]; %add zero radius point (which is not averaged)

% plot raw and averaged data
% figure;
% plot(D, G, '.')
% hold on;
% plot(r, GAvg, '-', 'linewidth', 2)

% plot number of averaged points per radius band
% figure;
% plot(r, NPoints,'-')
% xlabel('Radius [Pixels]')
% ylabel('Number of averaged points')


rCut = round(0.8*length(r));
r = r( 1 : rCut); %cut data at 80% from max radius - data might become noisy there
GAvg = GAvg( 1: rCut );
GSem = GSem( 1: rCut );

% rCutLow = 20;
% GAvg = GAvg( rCutLow<r );
% GSem = GSem( rCutLow<r );
% r = r( rCutLow<r ); %cut data at 80% from max radius - data might become noisy there


%replicate data left to right (includes center-zero point, as it it averaged with the other points in the first band)
GAvg = [fliplr(GAvg) GAvg];
GSem = [fliplr(GSem) GSem];
r = [fliplr(-r) r];


end



