classdef FitTypes
%FITTYPES Summary of this class goes here
%   Detailed explanation goes here

   properties ( Abstract = true, Constant = true )
       ID
   end
   
   properties (SetAccess = public)
                   
       atomsNo = -1;
       maxVal = -1;
       
       % the values are relative to the top/left corner of the image (no
       % chip start or ROI)
       ROILeft = -1;
       ROIRight = -1;
       ROITop = -1;
       ROIBottom = -1;
       
       xCenter = -1; % should be indexes (integers)
       yCenter= -1;
       
       xUnitSize = -1;
       yUnitSize = -1;
       
       xData = [];
       xStart = -1;
       yData = [];
       yStart = -1;
       
   end

   methods ( Abstract = true )
       appData                  = analyze                  (obj, appData) % do the analysis
%        [pic, xFit, yFit]          = getPicData             (obj, appData) % return:pic is normalized; xFit,yFit vectors of the fit (same length the original pic)
       normalizedROI       = getNormalizedROI  (obj, pic, x, y) % return normalized ROI (to the fitting constant). x,y are relative to pic.
       normalizedROI       = getTheoreticalROI  (obj, pic, x, y)
       
       normelizedPic        = normalizePic          (obj, pic) % normelized to the maxVal
       
       [xFit yFit]                 = getXYFitVectors     (obj, x, y) % x,y are absolute values. relative to the original fit.
                                          plotFitResults         (obj, appData) % plots the text. 
   end
   
   methods       
%         function pic = createSquare(obj, pic)
%             pic(obj.ROIBottom, obj.ROILeft : obj.ROIRight) = ones(1, obj.ROIRight - obj.ROILeft+1) * obj.lineColor;
%             pic(obj.ROITop, obj.ROILeft : obj.ROIRight) = ones(1, obj.ROIRight - obj.ROILeft+1) * obj.lineColor;
%             pic(obj.ROIBottom : obj.ROITop, obj.ROILeft) = ones(obj.ROITop - obj.ROIBottom+1, 1) * obj.lineColor;
%             pic(obj.ROIBottom : obj.ROITop, obj.ROIRight) = ones(obj.ROITop - obj.ROIBottom+1, 1) * obj.lineColor;
%         end
   end
end 
