
function LVData = changeValue(LVData, type, eventIndex, channelIndex, rampNo, newVal)
% return:
% LVData = []: failure

types = LVData.getTypes;

% newValue: a struct. 
% fields depends on the type, mandatory field: newVal.type
% type can be:
% eventStartTime,        fields: type, eventIndex, newVal
% analogStartTime,     fields: type, eventIndex, channelIndex, rampNo, newVal
% analogRampTime,        fields: type, eventIndex, channelIndex, rampNo, newVal
% anlaogValue,                fields: type, eventIndex, channelIndex, rampNo, newVal
% digitalTime,                fields: type, eventIndex, channelIndex, rampNo, newVal
% digitalValue,             fields: type, eventIndex, channelIndex, rampNo, newVal
% RFTime,                             fields: type, rampNo, newVal
% RFCommand,                     fields: type, rampNo, newVal

switch type
    case types.eventStartTime
        if length(LVData.Events) <= eventIndex
            LVData = [];
        else
            LVData.Events(eventIndex+1).Event_Start = newVal;
        end
        
    case types.analogStartTime        
        if ( length(LVData.Events) <= eventIndex || ...
                length(LVData.Events(eventIndex+1).Analog_Channels) <= channelIndex || ...
                length(LVData.Events(eventIndex+1).Analog_Channels{channelIndex+1}) <= rampNo )
            LVData = [];
        else
           LVData.Events(eventIndex+1).Analog_Channels{channelIndex+1}(rampNo+1).Start = newVal;
        end
        
    case types.analogRampTime      
        if ( length(LVData.Events) <= eventIndex || ...
                length(LVData.Events(eventIndex+1).Analog_Channels) <= channelIndex || ...
                length(LVData.Events(eventIndex+1).Analog_Channels{channelIndex+1}) <= rampNo )
            LVData = [];
        else
           LVData.Events(eventIndex+1).Analog_Channels{channelIndex+1}(rampNo+1).Rise = newVal;
        end
        
    case types.anlaogValue
        if ( length(LVData.Events) <= eventIndex || ...
                length(LVData.Events(eventIndex+1).Analog_Channels) <= channelIndex || ...
                length(LVData.Events(eventIndex+1).Analog_Channels{channelIndex+1}) <= rampNo )
            LVData = [];
        else
           LVData.Events(eventIndex+1).Analog_Channels{channelIndex+1}(rampNo+1).Value = newVal;
        end
        
    case types.digitalTime
        if ( length(LVData.Events) <= eventIndex || ...
                length(LVData.Events(eventIndex+1).Digital_Channels) <= channelIndex || ...
                length(LVData.Events(eventIndex+1).Digital_Channels{channelIndex+1}) <= rampNo )
            LVData = [];
        else
           LVData.Events(eventIndex+1).Digital_Channels{channelIndex+1}(rampNo+1).Start = newVal;
        end        
        
    case types.digitalValue
        if ( length(LVData.Events) <= eventIndex || ...
                length(LVData.Events(eventIndex+1).Digital_Channels) <= channelIndex || ...
                length(LVData.Events(eventIndex+1).Digital_Channels{channelIndex+1}) <= rampNo )
            LVData = [];
        else
           LVData.Events(eventIndex+1).Digital_Channels{channelIndex+1}(rampNo+1).Value = newVal;
        end        
        
    case types.RFTime
        if ( length(LVData.RF_Ramp.Data) <= rampNo )
            LVData = [];
        else
           LVData.RF_Ramp.Data(rampNo+1).Time_milli_sec = newVal;
        end        
        
    case types.RFCommand
        if ( length(LVData.RF_Ramp.Data) <= rampNo )
            LVData = [];
        else
           LVData.RF_Ramp.Data(rampNo+1).String = newVal;
        end        
        
    otherwise
        LVData = [];
end

