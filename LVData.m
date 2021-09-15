
classdef LVData
    
    properties (SetAccess = public)
        eventsData = [];
        analogProperties = [];
        digitalProperties = [];
        AOMProperties = []; 
        filepath = [];
        filename = [];
        RFRamps = [];
    end
    
    methods
          
        function i = getEventIndex(obj, eventName)
             for i = 1 : length(obj.eventsData)
                if strcmp(obj.eventsData(i).Name,eventName)
                    return
                end
            end
        end
        
        function i = getAnalogChIndex(obj, channelName)
             for i = 1 : length(obj.analogProperties)
                if strcmp(obj.analogProperties(i).Name,channelName)
                    return
                end
            end
        end
        
        function i = getDigitalChIndex(obj, channelName)
             for i = 1 : length(obj.digitalProperties)
                if strcmp(obj.digitalProperties(i).Name,channelName)
                    return
                end
            end
        end
        
        function str = getEventsNames(obj)
            str = cell(1, length(obj.eventsData));
            for i = 1 : length(str)
                str{i} = obj.eventsData(i).Name;
            end
        end
        function str = getAnalogNames(obj)
            str = cell(1, length(obj.analogProperties) + length(obj.AOMProperties));
            for i = 1 : length(obj.analogProperties)
                str{i} = obj.analogProperties(i).Name;
            end
            for i = 1 : length(obj.AOMProperties)
                str{length(obj.analogProperties) + i} = obj.AOMProperties(i).Name;
            end
        end 
        function str = getDigitalNames(obj)
            str = cell(1, length(obj.digitalProperties));
            for i = 1 : length(str)
                str{i} = obj.digitalProperties(i).Name;
            end
        end
        function str = getRFNames(obj)
            str = cell(1, length(obj.RFRamps));
            for i = 1 : length(str)
                str{i} = obj.RFRamps(i).GPIB;
            end
        end
        function num = getRampNum(obj, type, eventIndex, channelIndex)
            types = LVData.getTypes();
            num = -1;
            switch type
               case types.eventStartTime                    
               case {types.analogRampTime types.analogStartTime types.analogValue} 
                   num = length(obj.eventsData(eventIndex).Analog_Channels{channelIndex});
               case {types.digitalTime types.digitalValue }   
                   num = length(obj.eventsData(eventIndex).Digital_Channels{channelIndex});                   
               case {types.RFTime types.RFCommand}     
                   num = length(obj.RFRamps(1).Data);         
            end
        end
        
        function ret = writeLabview(obj,output_filename)

            output=fopen(output_filename,'wt+');
            if ( output == -1 )
                ret = 0;
                return;
            else
                ret = 1;
            end

            fprintf(output,['EN,     ' num2str(size(obj.eventsData,2)) ',\n']);
            fprintf(output,';;, Analog Channels Properties \n');
            fprintf(output,';;,    Pos, Number, Polling,  Scale ,  Shift , Start Value ,Max Value ,Min Value   , Name \n');

            for i=1:size(obj.analogProperties,2)
            fprintf(output,['AP,','%7.0f,','%7.0f,','%7.0f,',' %8.6f, ','%8.6f, ','%8.6f, ','%8.6f, ','%8.6f,','%s,\n'],...
                            obj.analogProperties(1,i).Channel,...
                            obj.analogProperties(1,i).Shtut,...
                            obj.analogProperties(1,i).Polling,...
                            obj.analogProperties(1,i).Shift,...
                            obj.analogProperties(1,i).Scale,...
                            obj.analogProperties(1,i).Start_Value,...
                            obj.analogProperties(1,i).Min,...
                            obj.analogProperties(1,i).Max,...
                            obj.analogProperties(1,i).Name);
            end

            for i=1:size(obj.AOMProperties,2)
                fprintf(output,['AO,','%7.0f,','%7.0f, ','%8.6f,','%10.6f,','%11.6f,','%s,\n'],...
                    obj.AOMProperties(1,i).Channel,...
                    obj.AOMProperties(1,i).Shtut,...
                    obj.AOMProperties(1,i).Start_Value,...
                    obj.AOMProperties(1,i).Min,...
                    obj.AOMProperties(1,i).Max,...
                    obj.AOMProperties(1,i).Name);
            end

            fprintf(output,';;, Digital Channels Properties \n');
            fprintf(output,';;,    Pos , Number,     Name \n');

            for i=1:size(obj.digitalProperties,2)
                fprintf(output,['DP,','%7.0f,','%7.0f,','%s,\n'],...
                    obj.digitalProperties(1,i).Channel,...
                    obj.digitalProperties(1,i).Shtut,...
                    obj.digitalProperties(1,i).Name);  
            end

            fprintf(output,';;, Events Properties \n');
            fprintf(output,';;,    Number,  Start Time, On/Off,    Name  \n');

            for i=1:size(obj.eventsData,2)
            fprintf(output,['EV,','%7.0f,','%16.0f,','%7.0f,','%s,\n'],...
                obj.eventsData(1,i).index,...
                obj.eventsData(1,i).Event_Start,...
                obj.eventsData(1,i).Execute_Event,...
                obj.eventsData(1,i).Name);
            end

            fprintf(output,';;,   Event,Channel,  Ramp,  On/Off   Start ,         Rise,     Value  \n');

            for i=1:size(obj.eventsData,2)
                    for j=1:size(obj.eventsData(1,i).Analog_Channels,2)
                        if ~isempty(obj.eventsData(1,i).Analog_Channels{1,j})
                            fprintf(output,['AR,','%7.0f,','%7.0f,','%7.0f,','%7.0f\n'],...
                                obj.eventsData(1,i).index,...
                                j-1,...
                                size(obj.eventsData(1,i).Analog_Channels{1,j},2),...
                                obj.eventsData(1,i).Analog_Channels{1,j}(1).Execute);

                            for k=1:size(obj.eventsData(1,i).Analog_Channels{1,j},2)
                                fprintf(output,['AD,','%7.0f,','%7.0f,','%7.0f,','%15.0f,','%15.0f, ','%8.6f, \n'],...
                                obj.eventsData(1,i).index,...
                                j-1,...
                                k-1,...
                                obj.eventsData(1,i).Analog_Channels{1,j}(k).Start,...
                                obj.eventsData(1,i).Analog_Channels{1,j}(k).Rise,...
                                obj.eventsData(1,i).Analog_Channels{1,j}(k).Value);
                            end
                        else
                            fprintf(output,['AR,','%7.0f,','%7.0f,','%7.0f,','%7.0f\n'],...
                                obj.eventsData(1,i).index,...
                                j-1,...
                                0,...
                                0);
                        end
                    end

                    for j=1:size(obj.eventsData(1,i).Digital_Channels,2)    
                        if ~isempty(obj.eventsData(1,i).Digital_Channels{1,j})
                            fprintf(output,['DR,','%7.0f,','%7.0f,','%7.0f,','%7.0f\n'],...
                                obj.eventsData(1,i).index,...
                                j-1,...
                                size(obj.eventsData(1,i).Digital_Channels{1,j},2),...
                                obj.eventsData(1,i).Digital_Channels{1,j}(1).Execute);

                            for k=1:size(obj.eventsData(1,i).Digital_Channels{1,j},2)
                                fprintf(output,['DD,','%7.0f,','%7.0f,','%7.0f,','%15.0f,','%7.0f\n'],...
                                obj.eventsData(1,i).index,...
                                j-1,...
                                k-1,...
                                obj.eventsData(1,i).Digital_Channels{1,j}(k).Start,...
                                obj.eventsData(1,i).Digital_Channels{1,j}(k).Value);
                            end
                        else
                            fprintf(output,['DR,','%7.0f,','%7.0f,','%7.0f,','%7.0f\n'],...
                                obj.eventsData(1,i).index,...
                                j-1,...
                                0,...
                                0);
                        end
                    end
            end
            
            for j = 1 : length(obj.RFRamps)
                fprintf(output,['TS,','%7.0f,','%7.0f,','%7.0f,','%7.0f,',  '%s,','%s\n'],...
                    obj.RFRamps(j).blabla,...
                    obj.RFRamps(j).ramp_index,...
                    obj.RFRamps(j).ramp_num,...
                    obj.RFRamps(j).blabla2,...
                    obj.RFRamps(j).GPIB,...
                    obj.RFRamps(j).Command);
            end
    %             fprintf(output,['TS,','%7.0f,','%7.0f,','%7.0f,','%9.0f,','%9.0f,',' %9.6f,','%7.0f\n'],...
    %                 obj.RFRamps.Ag33250,...
    %                 obj.RFRamps.Amplitude_units,...
    %                 obj.RFRamps.Offset,...
    %                 obj.RFRamps.Frequency_Deviation,...
    %                 obj.RFRamps.Carrier_Frequency,...
    %                 obj.RFRamps.Amplitude,...
    %                 obj.RFRamps.Shtut);

            for i=1:size(obj.RFRamps(1).Data,2)
%                 fprintf(output, [obj.RFRamps.Data(1,i).String '\n']);
                for j = 1 : length(obj.RFRamps)
                    fprintf(output,['TT,','%7.0f, ','%6.0f,','%7.0f,','%s\n'],...
                        obj.RFRamps(j).Data(1,i).blabla,...
                        obj.RFRamps(j).Data(1,i).Index,...
                        obj.RFRamps(j).Data(1,i).Time_milli_sec,...
                        obj.RFRamps(j).Data(1,i).String);
                end
%                 fprintf(output,['TT,','%7.0f, ','%7.0f, ','%s\n'],...
%                     obj.RFRamps.Data(1,i).Index,...
%                     obj.RFRamps.Data(1,i).Time_milli_sec,...
%                     obj.RFRamps.Data(1,i).String);
            end
            fclose(output);
        end
            
        function ret = changeValue(obj, type, eventIndex, channelIndex, rampNo, newVal)
            % return:
            % ret = []: failure

            types = LVData.getTypes();

            % newValue: a struct. 
            % fields depends on the type, mandatory field: newVal.type
            % type can be:
            % eventStartTime,        fields: type, eventIndex, newVal
            % analogStartTime,     fields: type, eventIndex, channelIndex, rampNo, newVal
            % analogRampTime,        fields: type, eventIndex, channelIndex, rampNo, newVal
            % analogValue,                fields: type, eventIndex, channelIndex, rampNo, newVal
            % digitalTime,                fields: type, eventIndex, channelIndex, rampNo, newVal
            % digitalValue,             fields: type, eventIndex, channelIndex, rampNo, newVal
            % RFTime,                             fields: type, rampNo, newVal
            % RFCommand,                     fields: type, rampNo, newVal

            switch type
                case types.eventStartTime
                    if length(obj.eventsData) < eventIndex
                        ret = [];
                        return
                    else
                        obj.eventsData(eventIndex).Event_Start = newVal;
                    end

                case types.analogStartTime        
                    if ( length(obj.eventsData) < eventIndex || ...
                            length(obj.eventsData(eventIndex).Analog_Channels) < channelIndex || ...
                            length(obj.eventsData(eventIndex).Analog_Channels{channelIndex}) < rampNo )
                        ret = [];
                        return
                    else
                       obj.eventsData(eventIndex).Analog_Channels{channelIndex}(rampNo).Start = newVal;
                    end

                case types.analogRampTime      
                    if ( length(obj.eventsData) < eventIndex || ...
                            length(obj.eventsData(eventIndex).Analog_Channels) < channelIndex || ...
                            length(obj.eventsData(eventIndex).Analog_Channels{channelIndex}) < rampNo )
                        ret = [];
                        return
                    else
                       obj.eventsData(eventIndex).Analog_Channels{channelIndex}(rampNo).Rise = newVal;
                    end

                case types.analogValue
                    if ( length(obj.eventsData) < eventIndex || ...
                            length(obj.eventsData(eventIndex).Analog_Channels) < channelIndex || ...
                            length(obj.eventsData(eventIndex).Analog_Channels{channelIndex}) < rampNo )
                        ret = [];
                        return
                    else
                       obj.eventsData(eventIndex).Analog_Channels{channelIndex}(rampNo).Value = newVal;
                    end

                case types.digitalTime
                    if ( length(obj.eventsData) < eventIndex || ...
                            length(obj.eventsData(eventIndex).Digital_Channels) < channelIndex || ...
                            length(obj.eventsData(eventIndex).Digital_Channels{channelIndex}) < rampNo )
                        ret = [];
                        return
                    else
                       obj.eventsData(eventIndex).Digital_Channels{channelIndex}(rampNo).Start = newVal;
                    end        

                case types.digitalValue
                    if ( length(obj.eventsData) < eventIndex || ...
                            length(obj.eventsData(eventIndex).Digital_Channels) < channelIndex || ...
                            length(obj.eventsData(eventIndex).Digital_Channels{channelIndex}) < rampNo )
                        ret = [];
                        return
                    else
                       obj.eventsData(eventIndex).Digital_Channels{channelIndex}(rampNo).Value = newVal;
                    end        

                case types.RFTime
                    if ( length(obj.RFRamps) < channelIndex || ...
                            length(obj.RFRamps(channelIndex).Data) < rampNo )
                        ret = [];
                        return
                    else
                       obj.RFRamps(channelIndex).Data(rampNo).Time_milli_sec = newVal;
                    end        

                case types.RFCommand
                    if ( length(obj.RFRamps) < channelIndex || ...
                            length(obj.RFRamps(channelIndex).Data) < rampNo )
                        ret = [];
                        return
                    else
                       obj.RFRamps(channelIndex).Data(rampNo).String = newVal;
                    end        

                otherwise
                    ret = [];
                    return
            end
            
            ret = obj;
        end
        
    end
    
    methods (Static = true)
      
        function newLVData=readLabview(filename)
            
            newLVData = [];

            %filename='C:\Ramon\RawData\Labview_data_to_load.txt';
            fileID=fopen(filename,'r');
            if fileID == -1
%                 errordlg(['No such file: ' filename],'Error', 'modal');
                return
            end
            i=1;
            RFC=1; % RF Counter

            while feof(fileID)==0  
                rawData{i}=fgetl(fileID); %#ok<*AGROW>
                i=i+1;
            end
            fclose(fileID);
            i = 1;
            while i <= length(rawData)
                if isempty(rawData{i})
                    i = i+1;
                    continue;
                end
                line=textscan(rawData{i}, '%s', 'delimiter',',');
                Data = [];
                for j=1:size(line{1,1},1)
                    Data{j}=line{1,1}{j,1};
                end
                switch Data{1}
                    case {'EN'}
            %             'Event Num.';
                        Event_Num=str2double(Data{2});
                        Events=struct('index',num2cell(0:Event_Num-1),'Event_Start',{'NaN'},'Execute_Event','NaN',...
                                'Name','NaN','Analog_Channels',cell(1),'Digital_Channels',cell(1));
                    case {'AP'}
            %             'Analog Prop.';
                        Analog_Properties(str2double(Data{2})+1).Shtut=str2double(Data{3});
                        Analog_Properties(str2double(Data{2})+1).Channel=str2double(Data{2});
                        Analog_Properties(str2double(Data{2})+1).Polling=str2double(Data{4});
                        Analog_Properties(str2double(Data{2})+1).Shift=str2double(Data{5});
                        Analog_Properties(str2double(Data{2})+1).Scale=str2double(Data{6});
                        Analog_Properties(str2double(Data{2})+1).Start_Value=str2double(Data{7});
                        Analog_Properties(str2double(Data{2})+1).Min=str2double(Data{8});
                        Analog_Properties(str2double(Data{2})+1).Max=str2double(Data{9});
                        Analog_Properties(str2double(Data{2})+1).Name=Data{10};
                    case {'DP'}
            %             'Digital Prop.';
                        Digital_Properties(str2double(Data{2})+1).Channel=str2double(Data(2));
                        Digital_Properties(str2double(Data{2})+1).Shtut=str2double(Data(3));
                        Digital_Properties(str2double(Data{2})+1).Name=Data{4};
                    case {'EV'}
            %             'Events Prop.';
                        Events(str2double(Data{2})+1).Event_Start=str2double(Data{3});
                        Events(str2double(Data{2})+1).Execute_Event=str2double(Data{4});
                        Events(str2double(Data{2})+1).Name=Data{5};
                    case {'AR'}
            %             'Analog Ramp';
                        A_Ramp = [];
                        exec = str2double(Data{5});
                        for j=1:str2double(Data{4})
                            line=textscan(rawData{i+j}, '%s', 'delimiter',',');
                            Data = [];
                            for k=1:size(line{1,1},1)
                                Data{k}=line{1,1}{k,1};
                            end
                            A_Ramp(j).Execute = exec;
                            A_Ramp(j).Start = str2double(Data{5});
                            if length(Data)<6
                                Data
                            end
                            A_Ramp(j).Rise = str2double(Data{6});
                            A_Ramp(j).Value = str2double(Data{7});
            %                 A_Ramp(j).String = rawData{i+j};
                        end % Ramp=(execute,Start,Rise,Value)

                        Events(str2double(Data{2})+1).Analog_Channels{str2double(Data{3})+1}=A_Ramp;
                         if ~isempty(j)
                            i=i+j;
                        end
                    case {'DR'}
            %             'Digital Ramp';
%                         if strcmp(Data{2}, '5')
%                             Data
%                         end
                        D_Ramp = [];
                        exec = str2double(Data{5});
                        for j=1:str2double(Data{4})                 
                            line=textscan(rawData{i+j}, '%s', 'delimiter',',');
                            Data = [];
                            for k=1:size(line{1,1},1)
                                Data{k}=line{1,1}{k,1};
                            end
                            D_Ramp (j).Execute = exec;
                            D_Ramp(j).Start = str2double(Data{5});
                            D_Ramp(j).Value = str2double(Data{6});
            %                 D_Ramp(j).String = rawData{i+j};
                        end % Ramp=(execute,Start,Value)
                        Events(str2double(Data{2})+1).Digital_Channels{str2double(Data{3})+1}=D_Ramp;
                        if ~isempty(j)
                            i=i+j;
                        end
                    case {'TS'}
            %             'RF Ramp';
                        RF_Ramp(str2double(Data{2})+1).blabla=str2double(Data{2});
                        RF_Ramp(str2double(Data{2})+1).ramp_index=str2double(Data{3});
                        RF_Ramp(str2double(Data{2})+1).ramp_num=str2double(Data{4});
                        RF_Ramp(str2double(Data{2})+1).blabla2=str2double(Data{5});
                        RF_Ramp(str2double(Data{2})+1).GPIB=Data{6};
                        RF_Ramp(str2double(Data{2})+1).Command=Data{7};
%                         RF_Ramp.Ag33250=str2double(Data{2});
%                         RF_Ramp.Amplitude_units=str2double(Data{3});
%                         RF_Ramp.Offset=str2double(Data{4});
%                         RF_Ramp.Frequency_Deviation=str2double(Data{5});
%                         RF_Ramp.Carrier_Frequency=str2double(Data{6});
%                         RF_Ramp.Amplitude=str2double(Data{7});
%                         RF_Ramp.Shtut=str2double(Data{8});
%                         RF_Ramp.Start_Time = -1;
                    case {'TT'}
                         RF_Ramp(str2double(Data{2})+1).Data(str2double(Data{3})+1).blabla = ...
                             str2double(Data{2});
                         RF_Ramp(str2double(Data{2})+1).Data(str2double(Data{3})+1).Index = ...
                             str2double(Data{3});
                         RF_Ramp(str2double(Data{2})+1).Data(str2double(Data{3})+1).Time_milli_sec=...
                             str2double(Data{4});
                         if length(Data) == 5 %check if GPIB command is not empty. Apparently the behavior of textscan changed between different matlab versions, so this line is necessary for the newer versions.
                             RF_Ramp(str2double(Data{2})+1).Data(str2double(Data{3})+1).String = ...
                                 Data{5};
                         else %if command is empty, put empty string
                             RF_Ramp(str2double(Data{2})+1).Data(str2double(Data{3})+1).String = ...
                                 '';
                         end
%                         RF_Ramp.Data(RFC).Index = str2double(Data{2});
%                         RF_Ramp.Data(RFC).Time_milli_sec=str2double(Data{3});
%                         RF_Ramp.Data(RFC).String = Data{4};
            %                  RF_Ramp.Data(RFC).String = rawData{i};

%                              if ( RF_Ramp.Start_Time == -1 )
%                                 RF_Ramp.Start_Time = RF_Ramp.Data(1).Time_milli_sec;
%                              end
                            RFC=RFC+1;
                    case {'AO'}
            %             'AOM Prop.';
                        AOM_Properties(str2double(Data{3})+1).Channel=str2double(Data{3});
                        AOM_Properties(str2double(Data{3})+1).Shtut=str2double(Data{2});
                        AOM_Properties(str2double(Data{3})+1).Start_Value=str2double(Data{4});
                        AOM_Properties(str2double(Data{3})+1).Min=str2double(Data{5});
                        AOM_Properties(str2double(Data{3})+1).Max=str2double(Data{6});
                        AOM_Properties(str2double(Data{3})+1).Name=Data{7};
                    case {'DD','AD',';;'}%'TT',
                        'Other';
                    case ''
                        1;
                    otherwise
                        'Otherwise';
                end
                i = i+1;
            end
            % i
            [path, name, ext] = fileparts(filename);
            newLVData = LVData();
            newLVData.eventsData=Events;
            newLVData.analogProperties=Analog_Properties;
            newLVData.digitalProperties=Digital_Properties;
            newLVData.AOMProperties=AOM_Properties;
            newLVData.filepath=path;
            newLVData.filename=[name ext];
            newLVData.RFRamps=RF_Ramp;
        end
            
        function types = getTypes
            % types struct (constants)
            types.str = {'Event Start Time' 'Analog Start Time' 'Analog Ramp Time' 'Analog Value' ...
                'Digital Time' 'Digital Value' 'RF Time' 'RF Command'};
            types.eventStartTime = 1;
            types.analogStartTime = 2;
            types.analogRampTime = 3;
            types.analogValue = 4;
            types.digitalTime = 5;
            types.digitalValue = 6;
            types.RFTime = 7;
            types.RFCommand = 8;
        end


    end
end