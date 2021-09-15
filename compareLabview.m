
function ret = compareLabview(LVData1, LVData2)
% ret = [] - files equals
% otherwise, ret = cell of differenses

ret = [];
k = 1;
% general checking
if length(LVData1.analogProperties) ~= length(LVData2.analogProperties)
    ret{k} = ['Not Equal No of Analog Properties: ' num2str(length(LVData1.analogProperties)) ' -|- ' num2str(length(LVData2.analogProperties))];
    k = k+1;
end
if length(LVData1.digitalProperties) ~= length(LVData2.digitalProperties)
    ret{k} = ['Not Equal No of Digital Properties: ' num2str(length(LVData1.digitalProperties)) ' -|- ' num2str(length(LVData2.digitalProperties))];
    k = k+1;
end
if length(LVData1.AOMProperties) ~= length(LVData2.AOMProperties)
    ret{k} = ['Not Equal No of AOM Properties: ' num2str(length(LVData1.AOMProperties)) ' -|- ' num2str(length(LVData2.AOMProperties))];
    k = k+1;
end
if length(LVData1.eventsData) ~= length(LVData2.eventsData)
    ret{k} = ['Not Equal No of eventsData: ' num2str(length(LVData1.eventsData)) ' -|- ' num2str(length(LVData2.eventsData))];
    k = k+1;
end
if length(LVData1.RFRamps.Data) ~= length(LVData2.RFRamps.Data)
    ret{k} = ['Not Equal No of RF Ramps: ' num2str(length(LVData1.RFRamps.Data)) ' -|- ' num2str(length(LVData2.RFRamps.Data))];
    k = k+1;
end
    
% file names
if  ~isequalwithequalnans (LVData1.filepath, LVData2.filepath)
    ret{k} = ['Directory Name: ' LVData1.filepath ' -|- ' LVData2.filepath];
    k = k+1;
end
if  ~isequalwithequalnans (LVData1.filename, LVData2.filename)
    ret{k} = ['File Name: ' LVData1.filename ' -|- ' LVData2.filename];
    k = k+1;
end

% analog properties
for i = 1: length(LVData1.analogProperties)
    if  ~isequalwithequalnans (LVData1.analogProperties(i), LVData2.analogProperties(i))
        ret{k} = ['Analog Properties No: ' num2str(i)];
        k = k+1;
    end
end

% digital properties
for i = 1: length(LVData1.digitalProperties)
    if  ~isequalwithequalnans (LVData1.digitalProperties(i), LVData2.digitalProperties(i))
        ret{k} = ['Digital Properties No: ' num2str(i)];
        k = k+1;
    end
end

% AOM properties
for i = 1: length(LVData1.AOMProperties)
    if  ~isequalwithequalnans (LVData1.AOMProperties(i), LVData2.AOMProperties(i))
        ret{k} = ['AOM Properties No: ' num2str(i)];
        k = k+1;
    end
end

% eventsData
for i = 1: length(LVData1.eventsData)
    if  ~isequalwithequalnans (LVData1.eventsData(i).index, LVData2.eventsData(i).index)
        ret{k} = ['Event No: ' num2str(i) ', different index: ' num2str(LVData1.eventsData(i).index) ' -|- ' num2str(LVData2.eventsData(i).index)];
        k = k+1;
    end
    if  ~isequalwithequalnans (LVData1.eventsData(i).Event_Start, LVData2.eventsData(i).Event_Start)
        ret{k} = ['Event No: ' num2str(i) ', different Event Start: ' num2str(LVData1.eventsData(i).Event_Start) ' -|- ' num2str(LVData2.eventsData(i).Event_Start)];
        k = k+1;
    end
    if  ~isequalwithequalnans (LVData1.eventsData(i).Execute_Event, LVData2.eventsData(i).Execute_Event)
        ret{k} = ['Event No: ' num2str(i) ', different Execute Event: ' num2str(LVData1.eventsData(i).Execute_Event) ' -|- ' num2str(LVData2.eventsData(i).Execute_Event)];
        k = k+1;
    end
    if  ~isequalwithequalnans (LVData1.eventsData(i).Name, LVData2.eventsData(i).Name)
        ret{k} = ['Event No: ' num2str(i) ', different Name: ' LVData1.eventsData(i).Name ' -|- ' LVData2.eventsData(i).Name];
        k = k+1;
    end
    
    % analog channels
    for j = 1 : length(LVData1.eventsData(i).Analog_Channels)
        for jj = 1 :  length(LVData1.eventsData(i).Analog_Channels{j})
            if  ~isequalwithequalnans (LVData1.eventsData(i).Analog_Channels{j}(jj).Execute, LVData2.eventsData(i).Analog_Channels{j}(jj).Execute)
                ret{k} = ['Event No: ' num2str(i) ', Analog Channel No: ' num2str(j) ', Ramp No: ' num2str(jj) ...
                    ', different Execute: ' num2str(LVData1.eventsData(i).Analog_Channels{j}(jj).Execute) ...
                    ' -|- ' num2str(LVData2.eventsData(i).Analog_Channels{j}(jj).Execute)];
                k = k+1;
            end
             if  ~isequalwithequalnans (LVData1.eventsData(i).Analog_Channels{j}(jj).Start, LVData2.eventsData(i).Analog_Channels{j}(jj).Start)
                ret{k} = ['Event No: ' num2str(i) ', Analog Channel No: ' num2str(j) ', Ramp No: ' num2str(jj) ...
                    ', different Start: ' num2str(LVData1.eventsData(i).Analog_Channels{j}(jj).Start) ...
                    ' -|- ' num2str(LVData2.eventsData(i).Analog_Channels{j}(jj).Start)];
                k = k+1;
             end
             if  ~isequalwithequalnans (LVData1.eventsData(i).Analog_Channels{j}(jj).Rise, LVData2.eventsData(i).Analog_Channels{j}(jj).Rise)
                ret{k} = ['Event No: ' num2str(i) ', Analog Channel No: ' num2str(j) ', Ramp No: ' num2str(jj) ...
                    ', different Rise: ' num2str(LVData1.eventsData(i).Analog_Channels{j}(jj).Rise) ...
                    ' -|- ' num2str(LVData2.eventsData(i).Analog_Channels{j}(jj).Rise)];
                k = k+1;
             end
             if  ~isequalwithequalnans (LVData1.eventsData(i).Analog_Channels{j}(jj).Value, LVData2.eventsData(i).Analog_Channels{j}(jj).Value)
                ret{k} = ['Event No: ' num2str(i) ', Analog Channel No: ' num2str(j) ', Ramp No: ' num2str(jj) ...
                    ', different Value: ' num2str(LVData1.eventsData(i).Analog_Channels{j}(jj).Value) ...
                    ' -|- ' num2str(LVData2.eventsData(i).Analog_Channels{j}(jj).Value)];
                k = k+1;
            end
        end
    end
    
    % digital channels
    for j = 1 : length(LVData1.eventsData(i).Digital_Channels)
        for jj = 1 :  length(LVData1.eventsData(i).Digital_Channels{j})
             if  ~isequalwithequalnans (LVData1.eventsData(i).Digital_Channels{j}(jj).Execute, LVData2.eventsData(i).Digital_Channels{j}(jj).Execute)
                ret{k} = ['Event No: ' num2str(i) ', Digital Channel No: ' num2str(j) ', Ramp No: ' num2str(jj) ...
                    ', different Execute: ' num2str(LVData1.eventsData(i).Digital_Channels{j}(jj).Execute) ...
                    ' -|- ' num2str(LVData2.eventsData(i).Digital_Channels{j}(jj).Execute)];
                k = k+1;
             end
            if  ~isequalwithequalnans (LVData1.eventsData(i).Digital_Channels{j}(jj).Start, LVData2.eventsData(i).Digital_Channels{j}(jj).Start)
                ret{k} = ['Event No: ' num2str(i) ', Digital Channel No: ' num2str(j) ', Ramp No: ' num2str(jj) ...
                    ', different Execute: ' num2str(LVData1.eventsData(i).Digital_Channels{j}(jj).Start) ...
                    ' -|- ' num2str(LVData2.eventsData(i).Digital_Channels{j}(jj).Start)];
                k = k+1;
            end
            if  ~isequalwithequalnans (LVData1.eventsData(i).Digital_Channels{j}(jj).Value, LVData2.eventsData(i).Digital_Channels{j}(jj).Value)
                ret{k} = ['Event No: ' num2str(i) ', Digital Channel No: ' num2str(j) ', Ramp No: ' num2str(jj) ...
                    ', different Execute: ' num2str(LVData1.eventsData(i).Digital_Channels{j}(jj).Value) ...
                    ' -|- ' num2str(LVData2.eventsData(i).Digital_Channels{j}(jj).Value)];
                k = k+1;
            end
        end
    end
end
    
% % RF properties
% if  ~isequalwithequalnans (LVData1.RFRamps.Ag33250, LVData2.RFRamps.Ag33250)
%     ret{k} = ['RF Ramps, different Ag33250: ' num2str(LVData1.RFRamps.Ag33250) ' -|- ' num2str(LVData2.RFRamps.Ag33250)];
%     k = k+1;
% end
% if  ~isequalwithequalnans (LVData1.RFRamps.Amplitude_units, LVData2.RFRamps.Amplitude_units)
%     ret{k} = ['RF Ramps, different Amplitude units: ' num2str(LVData1.RFRamps.Amplitude_units) ' -|- ' num2str(LVData2.RFRamps.Amplitude_units)];
%     k = k+1;
% end
% if  ~isequalwithequalnans (LVData1.RFRamps.Offset, LVData2.RFRamps.Offset)
%     ret{k} = ['RF Ramps, different Offset: ' num2str(LVData1.RFRamps.Offset) ' -|- ' num2str(LVData2.RFRamps.Offset)];
%     k = k+1;
% end
% if  ~isequalwithequalnans (LVData1.RFRamps.Frequency_Deviation, LVData2.RFRamps.Frequency_Deviation)
%     ret{k} = ['RF Ramps, different Frequency Deviation: ' num2str(LVData1.RFRamps.Frequency_Deviation) ' -|- ' num2str(LVData2.RFRamps.Frequency_Deviation)];
%     k = k+1;
% end
% if  ~isequalwithequalnans (LVData1.RFRamps.Carrier_Frequency, LVData2.RFRamps.Carrier_Frequency)
%     ret{k} = ['RF Ramps, different Carrier Frequency: ' num2str(LVData1.RFRamps.Carrier_Frequency) ' -|- ' num2str(LVData2.RFRamps.Carrier_Frequency)];
%     k = k+1;
% end
% if  ~isequalwithequalnans (LVData1.RFRamps.Amplitude, LVData2.RFRamps.Amplitude)
%     ret{k} = ['RF Ramps, different Amplitude: ' num2str(LVData1.RFRamps.Amplitude) ' -|- ' num2str(LVData2.RFRamps.Amplitude)];
%     k = k+1;
% end
% if  ~isequalwithequalnans (LVData1.RFRamps.Shtut, LVData2.RFRamps.Shtut)
%     ret{k} = ['RF Ramps, different Shtut: ' num2str(LVData1.RFRamps.Shtut) ' -|- ' num2str(LVData2.RFRamps.Shtut)];
%     k = k+1;
% end
% if  ~isequalwithequalnans (LVData1.RFRamps.Start_Time, LVData2.RFRamps.Start_Time)
%     ret{k} = ['RF Ramps, different Start Time: ' num2str(LVData1.RFRamps.Start_Time) ' -|- ' num2str(LVData2.RFRamps.Start_Time)];
%     k = k+1;
% end
% 
% % RF ramps 
% for j = 1 : length(LVData1.RFRamps.Data)
%     if  ~isequalwithequalnans (LVData1.RFRamps.Data(j).Index, LVData2.RFRamps.Data(j).Index)
%         ret{k} = ['RF Ramps, Ramp No: ' num2str(j) ', different Index: ' num2str(LVData1.RFRamps.Data(j).Index) ' -|- ' num2str(LVData2.RFRamps.Data(j).Index)];
%         k = k+1;
%     end
%     if  ~isequalwithequalnans (LVData1.RFRamps.Data(j).Time_milli_sec, LVData2.RFRamps.Data(j).Time_milli_sec)
%         ret{k} = ['RF Ramps, Ramp No: ' num2str(j) ', different Time_milli_sec: ' num2str(LVData1.RFRamps.Data(j).Time_milli_sec) ' -|- ' num2str(LVData2.RFRamps.Data(j).Time_milli_sec)]; %#ok<*AGROW>
%         k = k+1;
%     end
%      if  ~isequalwithequalnans( LVData1.RFRamps.Data(j).String, LVData2.RFRamps.Data(j).String)
%         ret{k} = ['RF Ramps, Ramp No: ' num2str(j) ', different String: ' LVData1.RFRamps.Data(j).String ' -|- ' LVData2.RFRamps.Data(j).String];
%         k = k+1;
%     end
% end

compareDlg('PromptString','Differences in files:',...
                'SelectionMode','single',...
                'ListSize', [min(size(strvcat(ret{:}), 2)*6, 800)+25, min(length(ret)*20, 600)+25], ...
                'ListString',ret); %#ok<VCAT>
    
% listdlg('PromptString','Differences in files:',...
%                 'SelectionMode','single',...
%                 'ListSize', [min(size(strvcat(ret{:}), 2)*5, 800)+25, min(length(ret)*17, 600)+25], ...
%                 'ListString',ret); %#ok<VCAT>
            
            
end





