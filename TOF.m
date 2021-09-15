
classdef TOF < Measurement
    %TOF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        TOFTimes = [];
        imagingTimes = [];
        trapReleaseTime = -1;
                
        type = -1;
        event = -1;
        channel = -1;
        ramp = -1;
        newVal = -1;
    end
    
    
     methods ( Static = true )
        function o = create(appData)
            o = TOF(appData);
        end
     end
     
    methods        
        function obj = TOF(appData)
            % first line - MUST
            obj = obj@Measurement(LVData.readLabview(appData.consts.defaultStrLVFile_Save));%appData.data.LVData);
            
            obj.type = obj.valueTypes.eventStartTime;
            obj.event = obj.LVData.getEventIndex('Imaging');
            obj.channel = -1;
            obj.ramp = -1;
            obj.newVal = -1;
            
            obj.baseFolder = appData.save.saveDir;
            obj.trapReleaseTime = obj.LVData.eventsData(obj.LVData.getEventIndex('Trap Release')).Event_Start;
            
            % last line - MUST
            obj = obj.initialize(appData);
        end
        
        function obj = initialize(obj, appData)           
            if obj.isInitialized == 0
                answer = measDlg({'Folder' 'Num Iterations' 'Imaging Time [ms]'  'Iterations Order'}, 'Time of Flight Measurements', [1 1 1 1]', ...
                    {obj.baseFolder appData.consts.loops.TOF.noIterations appData.consts.loops.TOF.TOFTimes appData.consts.loops.iterationsStr}, ...
                    appData.consts.loops.options, {'folder' 'edit' 'edit' 'popupmenu'}, {obj.baseFolder 0 0 1}); %#ok<*NBRAK>
                obj.isInitialized = 1;
            else
                answer = measDlg({'Folder' 'Num Iterations' 'Imaging Time [ms]'  'Iterations Order'}, 'Time of Flight Measurements', [1 1 1 1]', ...
                    {obj.baseFolder obj.noIterations vec2str(obj.TOFTimes) {'Iterate Measurement' 'Iterate Loop' 'Random Iterations'}}, ...
                    appData.consts.loops.options, {'folder' 'edit' 'edit' 'popupmenu'}, {obj.baseFolder 0 0 1}); %#ok<*NBRAK>
            end
            if isempty(answer)
%                 obj = [];
                return
            end
            obj.baseFolder = answer{1};
            obj.noIterations = str2double(answer{2});
            obj.TOFTimes = eval([ '[' answer{3}  ']' ]);
            if isempty(obj.TOFTimes)
                errordlg('Input must be an array','Error', 'modal');
            end
            obj.iterationsOrder = answer{4};
            
            TOF = obj.iterateVec(obj.TOFTimes, obj.iterationsOrder);
            
            obj.imagingTimes = obj.trapReleaseTime + TOF*1e3;
            obj.noMeasurements = length(TOF);
            
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
            
            newLVData = obj.LVData.changeValue(obj.type, obj.event, obj.channel, obj.ramp, ...
                    obj.imagingTimes(obj.position) );
            set(appData.ui.pmSaveParam, 'Value', appData.consts.saveParams.TOF);
            set(appData.ui.etParamVal, 'String', num2str((obj.imagingTimes(obj.position) - obj.trapReleaseTime)*1e-3));
            
        end
        
        function str = getMeasStr(obj, appData)
            str = [ appData.consts.availableLoops.str{appData.consts.availableLoops.TOF}  ' (' obj.getCurrMeas() ')'];
        end
    end
    
    
end

