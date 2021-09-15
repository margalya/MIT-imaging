
classdef LT < Measurement
    %LT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LTTimes = [];
        trapReleaseTimes = [];
        imagingTimes = [];
        trapReleaseTime = -1;
        imagingTime = -1;
                
        type = -1;
        trapReleaseEvent = -1;
        imagingEvent = -1;
        channel = -1;
        ramp = -1;
        newVal = -1;
    end
    
    
     methods ( Static = true )
        function o = create(appData)
            o = LT(appData);
        end
     end
     
    methods        
        function obj = LT(appData)
            % first line - MUST
            obj = obj@Measurement(LVData.readLabview(appData.consts.defaultStrLVFile_Save));%appData.data.LVData);
            
            obj.type = obj.valueTypes.eventStartTime;
            obj.trapReleaseEvent = obj.LVData.getEventIndex('Trap Release');
            obj.imagingEvent = obj.LVData.getEventIndex('Imaging');
            obj.channel = -1;
            obj.ramp = -1;
            obj.newVal = -1;
            
            obj.baseFolder = appData.save.saveDir;
            obj.trapReleaseTime = obj.LVData.eventsData(obj.trapReleaseEvent).Event_Start;
            obj.imagingTime = obj.LVData.eventsData(obj.imagingEvent).Event_Start;
            
            obj = obj.initialize(appData);
        end
        
        function obj = initialize(obj, appData)           
            if obj.isInitialized == 0
                answer = measDlg({'Folder' 'Num Iterations' 'Life Time (relative to first image) [ms]' 'Iterations Order'}, 'Life Time Measurements', [1 1 1 1]', ...
                    {obj.baseFolder appData.consts.loops.LT.noIterations appData.consts.loops.LT.LTTimes appData.consts.loops.iterationsStr}, ...
                    appData.consts.loops.options, {'folder' 'edit' 'edit' 'popupmenu'}, {obj.baseFolder 0 0 1}); %#ok<*NBRAK>
                obj.isInitialized = 1;
            else
                answer = measDlg({'Folder' 'Num Iterations' 'Life Time (relative to first image) [ms]'  'Iterations Order'}, 'Life Time Measurements', [1 1 1 1]', ...
                    {obj.baseFolder obj.noIterations vec2str(obj.LTTimes) appData.consts.loops.iterationsStr}, ...
                    appData.consts.loops.options, {'folder' 'edit' 'edit' 'popupmenu'}, {obj.baseFolder 0 0 1}); %#ok<*NBRAK>
            end
            if isempty(answer)
%                 obj = [];
                return
            end
            obj.baseFolder = answer{1};
            obj.noIterations = str2double(answer{2});
            obj.LTTimes = eval([ '[' answer{3}  ']' ]);
            if isempty(obj.LTTimes)
                errordlg('Input must be an array','Error', 'modal');
            end
            obj.iterationsOrder = answer{4};
            
            LT = obj.iterateVec(obj.LTTimes, obj.iterationsOrder);
            
            obj.trapReleaseTimes = obj.trapReleaseTime + LT*1e3;
            obj.imagingTimes = obj.imagingTime + LT*1e3;
            obj.noMeasurements = length(LT);
            
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
            
            newLVData = obj.LVData.changeValue(obj.type, obj.trapReleaseEvent, obj.channel, obj.ramp, ...
                    obj.trapReleaseTimes(obj.position) );
            newLVData = newLVData.changeValue(obj.type, obj.imagingEvent, obj.channel, obj.ramp, ...
                    obj.imagingTimes(obj.position) );
            set(appData.ui.pmSaveParam, 'Value', appData.consts.saveParams.darkTime);
            set(appData.ui.etParamVal, 'String', num2str(obj.trapReleaseTimes(obj.position) - obj.trapReleaseTime));
            
        end
        
        function str = getMeasStr(obj, appData)
            str = [ appData.consts.availableLoops.str{appData.consts.availableLoops.LT}  ' (' obj.getCurrMeas() ')'];
        end
    end
    
    
end

