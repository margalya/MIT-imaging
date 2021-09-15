classdef SternGerlach < Measurement
    %SternGerlach Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        TOFTime = 0;
        SGTRTime = 0;
        RFTime = 0;
%         SGTrapReleaseTime = 0;
        SGLTTimes = [];
        SGTrapReleaseTimes = [];
        
        typeEvents = -1;
        typeRFTime = -1;
        typeRFCommand = -1;
        SGTREvent = -1;
        imagingEvent = -1
        channel = -1;
        ramp = -1;
        newVal = -1;
    end
    
    methods ( Static = true )
        function o = create(appData)
            o = SternGerlach(appData);
        end
    end
     
    methods
        function obj = SternGerlach(appData)
            % first line - MUST
            obj = obj@Measurement(LVData.readLabview(appData.consts.defaultStrLVFile_Save));%appData.data.LVData);
            
            obj.typeEvents = obj.valueTypes.eventStartTime;            
            obj.typeRFTime = obj.valueTypes.RFTime;    
            obj.SGTREvent = obj.LVData.getEventIndex('SG Trap Release');
            obj.imagingEvent = obj.LVData.getEventIndex('Imaging');
            obj.channel = -1;
            obj.ramp = appData.consts.loops.SG.RFRampNo;
            obj.newVal = -1;
            
            obj.baseFolder = appData.save.saveDir;
            obj.SGTRTime = obj.LVData.eventsData(obj.SGTREvent).Event_Start;
            obj.TOFTime = obj.LVData.eventsData(obj.imagingEvent).Event_Start - obj.SGTRTime;
            obj.RFTime = obj.LVData.RFRamps.Data(obj.ramp).Time_milli_sec - floor(obj.SGTRTime/1e3);
            
            obj = obj.initialize(appData);
        end
        
        function obj = initialize(obj, appData)           
            if obj.isInitialized == 0
                answer = measDlg({'Folder' 'Num Iterations' 'Waiting Time (relative to first image) [ms]' 'Iterations Order'}, 'Stern Gerlach Measurements', [1 1 1 1]', ...
                    {obj.baseFolder appData.consts.loops.SG.noIterations appData.consts.loops.SG.SGTimes appData.consts.loops.iterationsStr}, ...
                    appData.consts.loops.options, {'folder' 'edit' 'edit' 'popupmenu'}, {obj.baseFolder 0 0 1}); %#ok<*NBRAK>
                obj.isInitialized = 1;
            else
                answer = measDlg({'Folder' 'Num Iterations' 'Waiting Time (relative to first image) [ms]'  'Iterations Order'}, 'Stern Gerlach Measurements', [1 1 1 1]', ...
                    {obj.baseFolder obj.noIterations vec2str(obj.SGLTTimes) appData.consts.loops.iterationsStr}, ...
                    appData.consts.loops.options, {'folder' 'edit' 'edit' 'popupmenu'}, {obj.baseFolder 0 0 1}); %#ok<*NBRAK>
            end
            if isempty(answer)
%                 obj = [];
                return
            end
            obj.baseFolder = answer{1};
            obj.noIterations = str2double(answer{2});
            obj.SGLTTimes = eval([ '[' answer{3}  ']' ]);
            if isempty(obj.SGLTTimes)
                errordlg('Input must be an array','Error', 'modal');
            end
            obj.iterationsOrder = answer{4};
            
            SGLT = obj.iterateVec(obj.SGLTTimes, obj.iterationsOrder);
            
            obj.SGTrapReleaseTimes = obj.SGTRTime + SGLT*1e3;
%             obj.imagingTimes = obj.imagingTime + LT*1e3;
            obj.noMeasurements = length(SGLT);
            
        end
                 
        function obj = edit(obj, appData)
%             tmpObj = obj;
%             obj.enableEdit = 0;
            obj = obj.initialize(appData);
%             if obj.noIterations == -1
%                 obj = tmpObj;
%             end
        end
        
        function [obj, newLVData] = next(obj, appData) 
            if obj.position == obj.noMeasurements
                newLVData = [];
                return
            end
            
            obj.position = obj.position + 1;
            
            newLVData = obj.LVData.changeValue(obj.typeEvents, obj.SGTREvent, -1, -1, ...
                    obj.SGTrapReleaseTimes(obj.position) );
            newLVData = newLVData.changeValue(obj.typeEvents, obj.imagingEvent, -1, -1, ...
                    obj.SGTrapReleaseTimes(obj.position) +  obj.TOFTime);
            newLVData = newLVData.changeValue(obj.typeRFTime, -1, -1, obj.ramp, ...
                    floor(obj.SGTrapReleaseTimes(obj.position)/1e3) +  obj.RFTime);
            
            set(appData.ui.pmSaveParam, 'Value', appData.consts.saveParams.darkTime);            
            set(appData.ui.etParamVal, 'String', num2str(obj.SGTrapReleaseTimes(obj.position)/1e3 - obj.SGTRTime/1e3));
            
        end
        
        function str = getMeasStr(obj, appData)
            str = [ appData.consts.availableLoops.str{appData.consts.availableLoops.SG}  ' (' obj.getCurrMeas() ')'];
        end
    end
    
end

