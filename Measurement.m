
classdef Measurement
    
    properties (SetAccess = public)
        LVData = [];
        baseFolder = '';
        isInitialized = 0;
        
        noIterations = -1;
        iterationsOrder = -1;
        noMeasurements = -1;
        position = 0;
        
        valueTypes = LVData.getTypes();
    end
    
    methods
        function obj = Measurement(LVData)
            obj.LVData = LVData;
        end
        function newVec = iterateVec(obj, vec, itarationsOrder)
            newVec = zeros([length(vec) obj.noIterations]);
%             if obj.noIterations == 1
%                 newVec = vec';
%             end
            for i = 1 : obj.noIterations
                newVec(:, i) = vec';
            end
           
            switch itarationsOrder
                case 1 %iterate measurement
                    newVec = newVec';
                    newVec = newVec(:);
                case 2 % iterate loop
%                     LT = LT';
                    newVec = newVec(:);
                case 3 % random iterations
                    newVec = newVec(:);
                    newVec = newVec(randperm(length(newVec)));
            end
        end
        function currMeas = getCurrMeas(obj) %return string 'current measurement / total measurements'
            currMeas = [ num2str(obj.position) '/' num2str(obj.noMeasurements)];
        end
    end
    
    
    methods ( Abstract = true, Static = true )
        o = create(appData);
    end
    
    methods ( Abstract = true )
        obj = initialize(obj, appData);
        obj = edit(obj, appData);
        [obj, LVData] = next(obj, appData);
        str = getMeasStr(obj, appData);
    end        
    
end