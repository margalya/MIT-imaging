
classdef CalcTheoretical < CalcAtomsNo
    %CALCREAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function atomsNo = calcAtomsNo(obj, appData, fitObj, pic, x, y)
            if strcmp(fitObj.ID,'Fit2DGaussian')
                hbar = 1.054571628e-34; %Planck's constant, [J*sec]
                f = 3.8422811521e+14; %F=2->F'=3 transition frequency, [Hz]
                Tim = 100e-6; %usual imaging pulse duration, [sec]
                G = 4095/5903; %approximate prosilica camera gain = (2^12-1)/{full well capacity}
                QE = 0.127; %approximate prosilica camera quantum efficiency at 780nm
                T = 0.884; %optical tranmission between chamber and CCD
                % define area of 2 sigma around the cloud's position
                %                 leftIndex = round(fitObj.yCenter-2*fitObj.sigmaY);
                %                 rightIndex = round(fitObj.yCenter+2*fitObj.sigmaY);
                %                 topIndex = round(fitObj.xCenter-2*fitObj.sigmaX);
                %                 bottomIndex = round(fitObj.xCenter+2*fitObj.sigmaX);
                % calculate mean photon number, 2 sigma around the cloud, from
                % 'withoutAtoms' picture
                %                 NpcMean = mean2( appData.data.plots{appData.consts.plotTypes.withoutAtoms}.pic( leftIndex:rightIndex, topIndex: bottomIndex) );
                % test command for the image area:
                % figure;imagesc( appData.data.plots{appData.consts.plotTypes.withAtoms}.pic( leftIndex:rightIndex, topIndex: bottomIndex) )
                % figure;imagesc( log(double(appData.data.plots{appData.consts.plotTypes.withoutAtoms}.pic( leftIndex:rightIndex, topIndex: bottomIndex)./appData.data.plots{appData.consts.plotTypes.withAtoms}.pic( leftIndex:rightIndex, topIndex: bottomIndex) ))); axis equal
                
                normalizedROI = fitObj.getNormalizedROI(pic, x, y);
                
                Npc = appData.data.plots{appData.consts.plotTypes.withoutAtoms}.pic( fitObj.ROITop : fitObj.ROIBottom, fitObj.ROILeft : fitObj.ROIRight); %number of photon counts per pixel, from 'withoutAtoms' image
                %                 figure;imagesc( appData.data.plots{appData.consts.plotTypes.withAtoms}.pic( fitObj.ROITop : fitObj.ROIBottom, fitObj.ROILeft : fitObj.ROIRight) )
                I = Npc.*(appData.data.camera.magnification^2 / T ) * hbar*2*pi*f*G / (appData.data.camera.xPixSz*appData.data.camera.yPixSz*Tim*QE); % evaluated intensity per pixel on the atoms
                scatcross = appData.consts.scatcross0 ./( 1 + 4*(4.28*1e6/appData.consts.linew)^2 + I/appData.consts.Isat);
                
                [X,Y] = meshgrid( 1:size(scatcross,2), 1:size(scatcross,1) );
                filt = ((Y-(fitObj.y0-fitObj.ROITop))/(fitObj.sigmaY)).^2+((X-(fitObj.x0-fitObj.ROILeft))/(fitObj.sigmaX)).^2 > 3^2;
%                 absorption = calcAbsoprtion(appData.data.plots{appData.consts.plotTypes.withAtoms}.pic( fitObj.ROITop : fitObj.ROIBottom, fitObj.ROILeft : fitObj.ROIRight),...
%                     appData.data.plots{appData.consts.plotTypes.withoutAtoms}.pic( fitObj.ROITop : fitObj.ROIBottom, fitObj.ROILeft : fitObj.ROIRight));
%                 figure;imagesc(absorption.*filt)
%                 absorption(~filt) = NaN;
                scatcross(filt) = appData.consts.scatcross0; %all values with distance of 3 sigma from cloud center in each direction are set to the usual value (to avoid affecting thew background calculation - 'C' parameter)
                
                atomsNo = round(appData.data.camera.xPixSz * appData.data.camera.yPixSz * sum(sum(normalizedROI ./ scatcross)));
            else
                atomsN  o = 0;
            end
        end
    end
end
