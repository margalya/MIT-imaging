classdef ROIUnits
%ROIUNITS Summary of this class goes here
%   Detailed explanation goes here

   properties %( SetAccess = private )
%        xSize = -1; % in units
%        ySize = -1; % in units
   end

   methods ( Abstract = true )
       [x0 y0 x1 y1] = getROICoords(obj, appData, fitObj)
       [x0 y0] = getCenter(obj, appData, fitObj)
   end
end 
