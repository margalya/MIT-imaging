
classdef CalcReal < CalcAtomsNo
%CALCREAL Summary of this class goes here
%   Detailed explanation goes here

   properties
   end

   methods 
       function [atomsNo, photonCount] = calcAtomsNo(obj, appData, fitObj, pic, x, y)
           normalizedROI = fitObj.getNormalizedROI(pic, x, y);
           scatcross = appData.consts.scatcross0{appData.options.atomType} * 1/(1+(appData.options.detuning*1e6*2/appData.consts.linew{appData.options.atomType})^2); 
           atomsNo = round(appData.data.camera.xPixSz * appData.data.camera.yPixSz * sum(sum(normalizedROI)) / scatcross);
           if isfield(appData.consts.cameras{appData.options.cameraType}, 'photonPerADU') %If camera is setup for fluorescence imaging
               photonCount = round(sum(normalizedROI(:)) * appData.consts.cameras{appData.options.cameraType}.photonPerADU); %Fluorescence imaging photon count
           end
       end
   end
end 
