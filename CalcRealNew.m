
classdef CalcReal < CalcAtomsNo
    %CALCREAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function atomsNo = calcAtomsNo(obj, appData, fitObj, pic, x, y)
            if ~strcmp(fitObj.ID, 'Fit1DGaussian')
                normalizedROI = fitObj.getNormalizedROI(pic, x, y);
                scatcross = appData.consts.scatcross0{appData.options.atomType} * 1/(1+(appData.options.detuning*1e6*2/appData.consts.linew{appData.options.atomType})^2);
                atomsNo = round(appData.data.camera.xPixSz * appData.data.camera.yPixSz * sum(sum(normalizedROI)) / scatcross);
                if isfield(appData.consts.cameras{appData.options.cameraType}, 'photonPerADU') %If camera is setup for fluorescence imaging
                    photonCount = round(sum(normalizedROI(:)) * appData.consts.cameras{appData.options.cameraType}.photonPerADU); %Fluorescence imaging photon count
                end
            else
                tempPic = pic(y,x);
                nsigma = 3.2; %number of cloud sigmas to delete in order to calculate the mean noise.
                [X,Y] = meshgrid( fitObj.ROILeft:fitObj.ROIRight, fitObj.ROITop:fitObj.ROIBottom );
                mask = ((X-fitObj.xCenter)/(fitObj.xRes.sigmaX*nsigma)).^2 + ((Y-fitObj.yCenter)/(fitObj.yRes.sigmaY*nsigma)).^2 < 1;
                tempPic(mask) = NaN;
%                 figure;imagesc(mask) %show mask
%                 figure;imagesc(tempPic); text(40,40,['Mean = ' num2str(nanmean(nanmean(tempPic)))]) %show picture after applying cloud removal mask
                % normalizedROI = pic(y, x) - nanmean(nanmean(tempPic)); %count all pixels in pic
                
                tempPic2 = pic(y,x);
                tempPic2(~mask) = NaN;
                % figure;imagesc(tempPic2) %show atoms after mask, no backround
                normalizedROI = tempPic2 - nanmean(nanmean(tempPic)); %only count pixels nsigma away from cloud center
                % display(nanmean(nanmean(tempPic)));
                %% Algorithm to obtain number of atoms at convergence vs nsigma
                % tic
                atomsPerCount = 193.5363; % atoms per digital count in picture = pixelSize^2 / scatcross.
                nsigmaInitial = 1;
                nsigma = nsigmaInitial; %initial value of nsigma
                nsigmaStep = 0.05;
                
                for i = 1 : 5 % : ((nsigmaMax - nsigmaMin)/nsigmaStep)
                    tempPic = pic(y,x);
                    mask = ((X-fitObj.xCenter)/(fitObj.xRes.sigmaX*nsigma)).^2 + ((Y-fitObj.yCenter)/(fitObj.yRes.sigmaY*nsigma)).^2 < 1;
                    tempPic(mask) = NaN;
                    tempPic2 = pic(y,x);
                    tempPic2(~mask) = NaN;
                    normalizedROItemp = tempPic2 - nanmean(nanmean(tempPic));
                    atomsNo(i) = round( atomsPerCount * nansum(nansum(normalizedROItemp)) );
                    nsigma = nsigma + nsigmaStep;
                end
                
                atomsNoDiff = mean(diff(atomsNo)) / nsigmaStep; % initial value of relative change in the number of atoms
                flag = 1;
                if atomsNoDiff < 1
                    errordlg('atomsNoDiff < 1');
                    return
                end
                while atomsNoDiff > 1 && flag < 100 %continue running while atomsNoDiff = d(AtomNo)/d(nsigma)<50, and number of iterations < 100
                    tempPic = pic(y,x);
                    mask = ((X-fitObj.xCenter)/(fitObj.xRes.sigmaX*nsigma)).^2 + ((Y-fitObj.yCenter)/(fitObj.yRes.sigmaY*nsigma)).^2 < 1;
                    tempPic(mask) = NaN;
                    tempPic2 = pic(y,x);
                    tempPic2(~mask) = NaN;
                    normalizedROItemp = tempPic2 - nanmean(nanmean(tempPic));
                    atomsNo = [ atomsNo round( atomsPerCount * nansum(nansum(normalizedROItemp)) ) ]; %append last value to the vector
                    
                    %calculate d(AtomNo)/d(nsigma) - the change of atoms number vs nsigma. The value is averaged over the last 5 values to reduce noise.
                    atomsNoDiff = mean(diff(atomsNo(end-5:end))) / nsigmaStep;
                    
                    nsigma = nsigma + nsigmaStep; %increase nsigma, for the next iteration.
                    if ~exist('atomsNoDiffV','var') %record the difference, for checking convergence
                        atomsNoDiffV = atomsNoDiff;
                    else
                        atomsNoDiffV = [atomsNoDiffV atomsNoDiff];
                    end
                    flag = flag+1; %increase exit flag
                    
                end
%                 figure; imagesc(tempPic2)
%                 figure; plot(atomsNoDiffV,'o'); title('atomsNoDiffV');
%                 
%                 figure; plot(nsigmaInitial : nsigmaStep : (nsigmaInitial + (length(atomsNo)-1)*nsigmaStep), atomsNo,'o')
%                 xlabel('Number of sigma''s');
%                 ylabel('Number of atoms');
%                 text(5.5,1.1e4,{['Atoms Num = ' num2str(round(mean(atomsNo(end-3:end)))) '+/-' num2str(round(std(atomsNo(end-10:end)))) ] ; ['X convergence at = ' num2str(round(fitObj.xRes.sigmaX*nsigma)) ' pixels'] ; ['Y convergence at = ' num2str(round(fitObj.yRes.sigmaY*nsigma)) ' pixels'] })
%                 title('Ellipse counting area');
                % display(toc);
                atomsNo = atomsNo(end); %return last value as calculated number of atoms
                %% Check for (in)convergence for different values of C
                % scatcross = 2.904896021313274e-13;
                % % c = nanmean(nanmean(tempPic));
                % for j=1:10; %run for different values of c
                %     c = nanmean(nanmean(tempPic))+(j-5)*1e-5;
                %     leg{j} = ['c = ' num2str(c)]; %legend values
                %     for i=1:25
                %         nsigma = i;
                %         mask = ((X-fitObj.xCenter)/(fitObj.xRes.sigmaX*nsigma)).^2 + ((Y-fitObj.yCenter)/(fitObj.yRes.sigmaY*nsigma)).^2 < 1;
                %         tempPic2 = pic(y,x);
                %         tempPic2(~mask) = NaN;
                %         %     figure;imagesc(tempPic2)
                %         normalizedROItemp = tempPic2 - c;
                %         atomsNo(i,j) = round(2.3000e-06 * 2.3000e-06 * nansum(nansum(normalizedROItemp)) / scatcross);
                % %         if i==5 && j==1
                % %             figure;imagesc(normalizedROItemp)
                % %         end
                %     end
                % end
                %
                % % figure; plot(0.1:0.1:10,c,'o')
                % figure; plot(atomsNo)
                % legend(leg);
                % % figure; plot(diff(atomsNo))
                
                %% check for the effect of different values of sigma on the value of c
                %             scatcross = 2.904896021313274e-13;
                %
                %             imin = 1;
                %             imax = 100;
                %             factor = 0.15;
                %
                %             for i = imin : imax
                %                 nsigma = i*factor;
                %                 tempPic = pic(y,x);
                %                 mask = ((X-fitObj.xCenter)/(fitObj.xRes.sigmaX*nsigma)).^2 + ((Y-fitObj.yCenter)/(fitObj.yRes.sigmaY*nsigma)).^2 < 1;
                %                 tempPic(mask) = NaN;
                %                 c(i) = nanmean(nanmean(tempPic));
                %                 tempPic2 = pic(y,x);
                %                 tempPic2(~mask) = NaN;
                %                 %     figure;imagesc(tempPic2)
                %                 normalizedROItemp = tempPic2 - nanmean(nanmean(tempPic));
                %                 atomsNo(i) = round(2.3000e-06 * 2.3000e-06 * nansum(nansum(normalizedROItemp)) / scatcross);
                %                 %     if i==100
                %                 %         figure;imagesc(tempPic2)
                %                 %     end
                %             end
                %
                %             figure; plot(linspace(imin*factor,imax*factor,length(imin:imax)),c,'o')
                %             xlabel('Number of sigma''s');
                %             ylabel('Value of c');
                %             title('Ellipse counting area');
                %             figure; plot(linspace(imin*factor,imax*factor,length(imin:imax)),atomsNo,'o')
                %             xlabel('Number of sigma''s');
                %             ylabel('Number of atoms');
                %             title('Ellipse counting area');
                %% Plot single value of c, for various nsigma
                scatcross = 3*(670.977338e-9)^2/2/pi;
%                 c = 0.0354; %value from 2DFermi Fit 
                c = 0.014929; %value from 2DFermi Fit 351
%                 c = -0.021271; %value from 2DFermi Fit 317
%                 c = 0.0189; %347
                imin = 1;
                imax = 40*5;
                factor = 0.05;
                
                for i = imin : imax
                    nsigma = i*factor;
                    mask = ((X-fitObj.xCenter)/(fitObj.xRes.sigmaX*nsigma)).^2 + ((Y-fitObj.yCenter)/(fitObj.yRes.sigmaY*nsigma)).^2 < 1;
                    tempPic2 = pic(y,x);
                    tempPic2(~mask) = NaN;
                    %     figure;imagesc(tempPic2)
                    normalizedROItemp = tempPic2 - c;
                    atomsNo(i) = round(6.45e-06 * 6.45e-06 * nansum(nansum(normalizedROItemp)) / scatcross);
                    if i==round(12.45/factor)
                        figure;imagesc(normalizedROItemp)
                    end
                end
                
                % figure; plot(0.1:0.1:10,c,'o')
                x = linspace(imin*factor,imax*factor,length(imin:imax));
                figure; plot(x, atomsNo,'o') %linspace(1/2,70/2,70),
                xlabel('N_\sigma - size of counting area, in number of cloud sigma''s');
                ylabel('Number of atoms');
                hold on;
                diffAtoms = diff(atomsNo);
                plot(x(diffAtoms==0), atomsNo(diffAtoms==0), '-black')
                plot(x(diffAtoms~=0), (erf(x(diffAtoms~=0)./sqrt(2))).^2*max(atomsNo),'-r')
                
                atomFrac = atomsNo(diffAtoms~=0)./max(atomsNo(diffAtoms~=0));
                figure; plot(x(diffAtoms~=0), atomFrac,'o')
                xlabel('N_\sigma - size of counting area, in number of cloud sigma''s');
                ylabel('Fraction from atoms counted from maximum')
%                 hold on; plot(x(diffAtoms~=0), (erf(x(diffAtoms~=0)./sqrt(2))).^2,'r'); % Cumulative distribution function of a 2D Gaussian is not clear... erf^2 seems only an approximation
%                 [~,crossing] = min(abs(atomFrac-0.995)); % does not work well - depends strongly on the maximum position, and hence also on the value of c
%                 text(0.5,0.5,['99.5% crossing =~ ' num2str(x(crossing)) 'sigmas'])
            end
            
        end
    end
end
