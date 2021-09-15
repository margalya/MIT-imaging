
classdef CalcTheoreticalFull < CalcAtomsNo
%CALCREAL Summary of this class goes here
%   Detailed explanation goes here

   properties
   end

   methods 
    function atomsNo = calcAtomsNo(obj, appData, fitObj, pic, x, y)
        xSize = appData.data.ROISizeX;
        ySize = appData.data.ROISizeY;
        x0 = fitObj.xCenter - xSize*fitObj.xUnitSize;
        x1 = fitObj.xCenter + xSize*fitObj.xUnitSize;
        y0 = fitObj.yCenter - ySize*fitObj.yUnitSize;
        y1 = fitObj.yCenter + ySize*fitObj.yUnitSize;         
        normalizedROI = fitObj.getTheoreticalROI(pic, round(x0:x1), round(y0:y1));
        scatcross = appData.consts.scatcross0 * 1/(1+(appData.options.detuning*1e6*2/appData.consts.linew)^2); 
        atomsNo = round(appData.data.camera.xPixSz * appData.data.camera.yPixSz * sum(sum(normalizedROI)) / scatcross);   
    end
   end
end 