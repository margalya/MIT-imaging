classdef GeneralLoop < Measurement
    %GENERALLOOP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % consts
%         GUI = 0; %% MUST be first
%         tmpObj = [];
        
%         appData = [];
        folderIcon = [];
        fontSize = -1;
        componentHeight = -1;
        fPause = 1;
        baseBaseFolder = '';
        saveFolder = '';
        
%         tmpSaveName = 'tmpGeneralLoop.mat';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ui controls
        win = -1;
        pbOpenReadDir = -1;
        etReadDir = -1;
        etNumIterations = -1;
        pmIterationsOrder = -1;
        pbOK = -1;
        pbCancel = -1;
        pbSave = -1;
        pbLoad = -1;
        tbEnableEdit = -1;
        pbAddLoop = -1;
        pbRemoveLastLoop = -1;
        
        pmASaveParam = [];
        etASaveParamVals = [];
        sANoElements = [];
        
        pALoops = [];
        pbAAddChannel = [];
        pbARemoveLastChannel = [];
        pmMChangeTypes = {}; %2D cell array
        pmMEventName = {};
        pmMChannelName = {};     
        pmMRumpNo = {};
        etMVector = {};
        sMNoElements = {};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % data
        runs = [];
        
        enableEdit = 1;
        noLoops = 1;
        pALoopsStartHeight = 0;
        noChannels = 1; % an array at the length of noLoops
        
        MChangeTypes = {}; %2D cell array
        MEventName = {};
        MChannelName = {};     
        MRumpNo = {};
        MVector = {};
        MVectorRFCommands = {};
        MVectorStr = {};
        MNoElements = {};
        
        saveParam = -1;
        saveParamsStr = '';
        saveOtherParamStr = '';
        
        changeTypes = {};
        ASaveParam = [];
        ASaveParamsStr = {};
        ASaveParamVals = {};
        ASaveParamValsStr = {};
        ASaveOtherParamStr = {''};
        ANoElements = {};
        
        eventNames = '';
        analogNames = '';
        digitalNames = '';
        RFNames = '';
        
    end
    
    
     methods ( Static = true )
        function o = create(appData)
            o = GeneralLoop(appData);
        end
     end
     
    methods        
        function obj = GeneralLoop(appData)
            % first line - MUST
            obj = obj@Measurement(LVData.readLabview(appData.consts.defaultStrLVFile_Save));%appData.data.LVData);
            
%             obj.appData = appData;
            
            obj.LVData = appData.data.LVData;
            obj.baseBaseFolder = appData.save.saveDir;
            
            obj.noIterations = 1;
            obj.iterationsOrder = 1;            
            
            obj.fontSize = appData.consts.fontSize;
            obj.componentHeight = appData.consts.componentHeight;
            obj.folderIcon = imread('folder28.bmp');
            obj.saveFolder = appData.consts.loops.GenLoop.saveFolder;
            
            obj.changeTypes = LVData.getTypes();
            obj.saveParam = appData.consts.saveParams.default;
            obj.saveParamsStr = appData.consts.saveParams.str;
            obj.saveOtherParamStr = appData.consts.saveOtherParamStr;
            
            obj.eventNames = obj.LVData.getEventsNames();
            obj.analogNames = obj.LVData.getAnalogNames();
            obj.digitalNames = obj.LVData.getDigitalNames();
            obj.RFNames = obj.LVData.getRFNames();
           
            % last line - MUST
            obj = obj.initialize(appData, 1);
%             obj.tmpObj = obj;
        end
        
        function obj = initialize(obj, appData, doWait)  %#ok<INUSL>
%             flag = 0;
            if nargin < 3
                doWait = 1;
%                 flag = 1;
            end
%             if doWait == -1
%                 flag = 1;
%             end
            if obj.enableEdit
                enableStr = 'on';
            else
                enableStr = 'off';
            end
            if obj.win == -1 || ~ishandle(obj.win)
                obj.win = figure('Visible', 'off', ...
                    'Name', 'General Loop', ...
                    'Units', 'Pixels', ...
                    'Position', [100 100 730 500], ...
                    'Resize', 'off', ...
                    'MenuBar', 'None', ...
                    'Toolbar', 'None', ...
                    'NumberTitle', 'off' , ...  
                    'WindowStyle'      ,'normal', ...
                    'HandleVisibility', 'callback');
            end            
            
            obj.pbOpenReadDir = uicontrol(obj.win, ...
                'Style', 'pushbutton', ...
                'String', '', ...
                'Units', 'pixels', ...
                'Enable', enableStr, ...
                'Position', [5 5 30 30], ...
                'BackgroundColor', [0.8 0.8 0.8], ...
                'CData', obj.folderIcon);      
            obj.etReadDir = uicontrol(obj.win, ...
                'Style', 'edit', ...
                'String', obj.baseBaseFolder, ...
                'Units', 'pixels', ...
                'Enable', enableStr, ...
                'Position', [40 5 355 obj.componentHeight+5], ...
                'BackgroundColor', 'white', ...
                'HorizontalAlignment', 'left', ...
                'FontSize', obj.fontSize); 
            obj.tbEnableEdit = uicontrol(obj.win, ...
                'Style', 'togglebutton', ...
                'String', 'Enable Edit', ...
                'Units', 'pixels', ...
                'Value', obj.enableEdit, ...
                'Position', [400 5 100 obj.componentHeight+5], ...
                'BackgroundColor', [0.8 0.8 0.8], ...
                'FontSize', obj.fontSize); 
            
%             st1 = uicontrol(obj.win, ...
%                 'Style', 'text', ...
%                 'String', 'Num Iterations:', ...
%                 'Units', 'pixels', ...
%                 'Enable', enableStr, ...
%                 'Position', [400 5 150 obj.componentHeight], ...
%                 'BackgroundColor', [0.8 0.8 0.8], ...
%                 'HorizontalAlignment', 'left', ...
%                 'FontSize', obj.fontSize);  %#ok<NASGU>
%             obj.etNumIterations = uicontrol(obj.win, ...
%                 'Style', 'edit', ...
%                 'String', '1', ... % num2str(obj.noIterations), ...
%                 'Units', 'pixels', ...
%                 'Enable', enableStr, ...
%                 'Position', [500 5 20 obj.componentHeight+5], ...
%                 'BackgroundColor', 'white', ...
%                 'HorizontalAlignment', 'left', ...
%                 'FontSize', obj.fontSize);            
%             obj.pmIterationsOrder = uicontrol(obj.win, ...
%                 'Style', 'popupmenu', ...
%                 'String', {'Iterate Measurement' 'Iterate Loop' 'Random Iterations'}, ...
%                 'Value', obj.iterationsOrder, ...
%                 'Units', 'pixels', ...
%                 'Enable', enableStr, ...
%                 'Position', [525 12  150 20], ...
%                 'BackgroundColor', 'white', ...
%                 'HorizontalAlignment', 'left', ...
%                 'FontSize', obj.fontSize); 
%             obj.pbOK = uicontrol(obj.win, ...
%                 'Style', 'pushbutton', ...
%                 'String', 'OK', ...
%                 'Units', 'pixels', ...
%                 'Enable', enableStr, ...
%                 'Position', [750 40 100 obj.componentHeight+5], ...
%                 'BackgroundColor', [0.8 0.8 0.8], ...
%                 'FontSize', obj.fontSize); 
            
            
            obj.pbOK = uicontrol(obj.win, ...
                'Style', 'pushbutton', ...
                'String', 'OK', ...
                'Units', 'pixels', ...
                'Enable', enableStr, ...
                'Position', [615 5 100 obj.componentHeight+5], ...
                'BackgroundColor', [0.8 0.8 0.8], ...
                'FontSize', obj.fontSize); 
            obj.pbCancel = uicontrol(obj.win, ...
                'Style', 'pushbutton', ...
                'String', 'Cancel', ...
                'Units', 'pixels', ...
                'Position', [505 5 100 obj.componentHeight+5], ...
                'BackgroundColor', [0.8 0.8 0.8], ...
                'FontSize', obj.fontSize); 
            
            obj.pbSave = uicontrol(obj.win, ...
                'Style', 'pushbutton', ...
                'String', 'Save Loop', ...
                'Units', 'pixels', ...
                'Position', [5 40 100 obj.componentHeight+5], ...
                'BackgroundColor', [0.8 0.8 0.8], ...
                'FontSize', obj.fontSize); 
            obj.pbLoad = uicontrol(obj.win, ...
                'Style', 'pushbutton', ...
                'String', 'Load Loop', ...
                'Units', 'pixels', ...
                'Enable', enableStr, ...
                'Position', [110 40 100 obj.componentHeight+5], ...
                'BackgroundColor', [0.8 0.8 0.8], ...
                'FontSize', obj.fontSize); 
%             obj.tbEnableEdit = uicontrol(obj.win, ...
%                 'Style', 'togglebutton', ...
%                 'String', 'Enable Edit', ...
%                 'Units', 'pixels', ...
%                 'Value', obj.enableEdit, ...
%                 'Position', [215 40 100 obj.componentHeight+5], ...
%                 'BackgroundColor', [0.8 0.8 0.8], ...
%                 'FontSize', obj.fontSize); 
%                 
%             st1 = uicontrol(obj.win, ...
%                 'Style', 'text', ...
%                 'String', 'Save Param:', ...
%                 'Units', 'pixels', ...
%                 'Enable', enableStr, ...
%                 'Position', [320 40 90 obj.componentHeight], ...
%                 'BackgroundColor', [0.8 0.8 0.8], ...
%                 'HorizontalAlignment', 'left', ...
%                 'FontSize', obj.fontSize);          %#ok<NASGU>
%             obj.pmSaveParam = uicontrol(obj.win, ...
%                 'Style', 'popupmenu', ...
%                 'String', obj.saveParamsStr, ...
%                 'Value', obj.saveParam, ...
%                 'Units', 'pixels', ...
%                 'Enable', enableStr, ...
%                 'Position', [405 47  85 20], ...
%                 'BackgroundColor', 'white', ...
%                 'HorizontalAlignment', 'left', ...
%                 'FontSize', obj.fontSize); 
%             st1 = uicontrol(obj.win, ...
%                 'Style', 'text', ...
%                 'String', 'Param Vals:', ...
%                 'Units', 'pixels', ...
%                 'Enable', enableStr, ...
%                 'Position', [495 40 100 obj.componentHeight], ...
%                 'BackgroundColor', [0.8 0.8 0.8], ...
%                 'HorizontalAlignment', 'left', ...
%                 'FontSize', obj.fontSize);          %#ok<NASGU>
%             obj.etSaveParamVals = uicontrol(obj.win, ...
%                 'Style', 'edit', ...
%                 'String', '0', ... 
%                 'Units', 'pixels', ...
%                 'Enable', enableStr, ...
%                 'Position', [575 40 80 obj.componentHeight+5], ...
%                 'BackgroundColor', 'white', ...
%                 'HorizontalAlignment', 'left', ...
%                 'FontSize', obj.fontSize);  
%             obj.sNoElements = uicontrol(obj.win, ...
%                 'Style', 'text', ...
%                 'String', obj.noElements, ...
%                 'Units', 'pixels', ...
%                 'Enable', enableStr, ...
%                 'Position', [660 40 85 obj.componentHeight], ...
%                 'BackgroundColor', [0.8 0.8 0.8], ...
%                 'HorizontalAlignment', 'left', ...
%                 'FontSize', obj.fontSize);  

            st1 = uicontrol(obj.win, ...
                'Style', 'text', ...
                'String', 'Num Iterations:', ...
                'Units', 'pixels', ...
                'Enable', enableStr, ...
                'Position', [215 40 100 obj.componentHeight], ...
                'BackgroundColor', [0.8 0.8 0.8], ...
                'HorizontalAlignment', 'left', ...
                'FontSize', obj.fontSize);  %#ok<NASGU>
            obj.etNumIterations = uicontrol(obj.win, ...
                'Style', 'edit', ...
                'String', '1', ... % num2str(obj.noIterations), ...
                'Units', 'pixels', ...
                'Enable', enableStr, ...
                'Position', [315 40 20 obj.componentHeight+5], ...
                'BackgroundColor', 'white', ...
                'HorizontalAlignment', 'left', ...
                'FontSize', obj.fontSize);            
            obj.pmIterationsOrder = uicontrol(obj.win, ...
                'Style', 'popupmenu', ...
                'String', {'Iterate Measurement' 'Iterate Loop' 'Random Iterations'}, ...
                'Value', obj.iterationsOrder, ...
                'Units', 'pixels', ...
                'Enable', enableStr, ...
                'Position', [340 47  150 20], ...
                'BackgroundColor', 'white', ...
                'HorizontalAlignment', 'left', ...
                'FontSize', obj.fontSize); 
            
            obj.pbAddLoop = uicontrol(obj.win, ...
                'Style', 'pushbutton', ...
                'String', 'Add Loop', ...
                'Units', 'pixels', ...
                'Enable', enableStr, ...
                'Position', [615 40 100 obj.componentHeight+5], ...
                'BackgroundColor', [0.8 0.8 0.8], ...
                'FontSize', obj.fontSize); 
%             obj.pbRemoveLastLoop = uicontrol(obj.win, ...
%                 'Style', 'pushbutton', ...
%                 'String', 'Remove Last', ...
%                 'Units', 'pixels', ...
%                 'Enable', enableStr, ...
%                 'Position', [860 5 100 obj.componentHeight+5], ...
%                 'BackgroundColor', [0.8 0.8 0.8], ...
%                 'FontSize', obj.fontSize); 
            obj.pbRemoveLastLoop = uicontrol(obj.win, ...
                'Style', 'pushbutton', ...
                'String', 'Remove Last', ...
                'Units', 'pixels', ...
                'Enable', enableStr, ...
                'Position', [505 40 100 obj.componentHeight+5], ...
                'BackgroundColor', [0.8 0.8 0.8], ...
                'FontSize', obj.fontSize); 
            
            obj.pALoopsStartHeight = 0;
            for i = 1 : obj.noLoops
                panelTitle = '';
                if i == 1
                    panelTitle = 'Inner Loop';
                else
                    for j = obj.noLoops-i+2 : obj.noLoops
                        panelTitle = [panelTitle 'Outer ']; %#ok<AGROW>
                    end
                    panelTitle = [panelTitle 'Loop']; %#ok<AGROW>
                end
                loopHeaight = max( 2*(obj.componentHeight+5)+(2+4)*5, ... %components and spaces
                    (obj.noChannels(i)+1)*(obj.componentHeight+5)+(obj.noChannels(i)+4)*5);
                obj.pALoops(i) = uipanel('Parent', obj.win, ...
                    'Title', panelTitle, ...
                    'TitlePosition', 'lefttop', ...
                    'Units', 'pixels', ...
                    'Position', [5 70+obj.pALoopsStartHeight  720 loopHeaight], ... 
                    'BackgroundColor', [0.8 0.8 0.8], ...
                    'FontSize', obj.fontSize);
                obj.pALoopsStartHeight = obj.pALoopsStartHeight +loopHeaight;
                            
%                 obj.pbAAddChannel(i) = uicontrol(obj.pALoops(i), ...
%                     'Style', 'pushbutton', ...
%                     'String', 'Add Channel', ...
%                     'Units', 'pixels', ...
%                     'Enable', enableStr, ...
%                     'Position', [855 obj.componentHeight+15 100 obj.componentHeight+5], ...
%                     'BackgroundColor', [0.8 0.8 0.8], ...
%                     'FontSize', obj.fontSize ); 

                st1 = uicontrol(obj.pALoops(i), ...
                    'Style', 'text', ...
                    'String', 'Save Param:', ...
                    'Units', 'pixels', ...
                    'Enable', enableStr, ...
                    'Position', [5 5 90 obj.componentHeight], ...
                    'BackgroundColor', [0.8 0.8 0.8], ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize', obj.fontSize);          %#ok<NASGU>
                obj.pmASaveParam(i) = uicontrol(obj.pALoops(i), ...
                    'Style', 'popupmenu', ...
                    'String', obj.saveParamsStr, ...
                    'Value', obj.saveParam, ...
                    'Units', 'pixels', ...
                    'Enable', enableStr, ...
                    'Position', [90 12  85 20], ...
                    'BackgroundColor', 'white', ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize', obj.fontSize);
%                 obj.ASaveParamsStr{i} = obj.saveParamsStr;
%                 obj.ASaveParam(i) = obj.saveParam;
                st1 = uicontrol(obj.pALoops(i), ...
                    'Style', 'text', ...
                    'String', 'Param Vals:', ...
                    'Units', 'pixels', ...
                    'Enable', enableStr, ...
                    'Position', [180 5 100 obj.componentHeight], ...
                    'BackgroundColor', [0.8 0.8 0.8], ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize', obj.fontSize);          %#ok<NASGU>
                obj.etASaveParamVals(i) = uicontrol(obj.pALoops(i), ...
                    'Style', 'edit', ...
                    'String', '0', ...
                    'Units', 'pixels', ...
                    'Enable', enableStr, ...
                    'Position', [260 5 80 obj.componentHeight+5], ...
                    'BackgroundColor', 'white', ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize', obj.fontSize);
%                 obj.ASaveParamVals{i} = [];
                obj.sANoElements(i) = uicontrol(obj.pALoops(i), ...
                    'Style', 'text', ...
                    'String', '0 Elements', ...
                    'Units', 'pixels', ...
                    'Enable', enableStr, ...
                    'Position', [345 5 100 obj.componentHeight], ...
                    'BackgroundColor', [0.8 0.8 0.8], ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize', obj.fontSize);
%                 obj.ANoElements{i} = '0 Elements';
                
                obj.pbAAddChannel(i) = uicontrol(obj.pALoops(i), ...
                    'Style', 'pushbutton', ...
                    'String', 'Add Channel', ...
                    'Units', 'pixels', ...
                    'Enable', enableStr, ...
                    'Position', [610 5 100 obj.componentHeight+5], ...
                    'BackgroundColor', [0.8 0.8 0.8], ...
                    'FontSize', obj.fontSize ); 
                obj.pbARemoveLastChannel(i) = uicontrol(obj.pALoops(i), ...
                    'Style', 'pushbutton', ...
                    'String', 'Remove Last', ...
                    'Units', 'pixels', ...
                    'Enable', enableStr, ...
                    'Position', [500 5 100 obj.componentHeight+5], ...
                    'BackgroundColor', [0.8 0.8 0.8], ...
                    'FontSize', obj.fontSize); 
                
                for j = 1 : obj.noChannels(i)
                    %%%%%%%%%%%%%%%%%%%%%%%
                    % add a line (change channel)  
                    channelHeight = 42 +(j-1)*30;
                    obj.pmMChangeTypes{i}{j} = uicontrol(obj.pALoops(i), ...
                        'Style', 'popupmenu', ...
                        'String', obj.changeTypes.str, ...
                        'Value', 1, ...
                        'Units', 'pixels', ...
                        'Enable', enableStr, ...
                        'Position', [5 channelHeight  120 20], ...
                        'BackgroundColor', 'white', ...
                        'HorizontalAlignment', 'left', ...
                        'FontSize', obj.fontSize);
%                     obj.MChangeTypes{i}{j} = 1;
                    obj.pmMEventName{i}{j} = uicontrol(obj.pALoops(i), ...
                        'Style', 'popupmenu', ...
                        'String', 'Events Names', ...
                        'Value', 1, ...
                        'Units', 'pixels', ...
                        'Enable', enableStr, ...
                        'Position', [130 channelHeight  120 20], ...
                        'BackgroundColor', 'white', ...
                        'HorizontalAlignment', 'left', ...
                        'FontSize', obj.fontSize);
%                     obj.MEventName{i}{j} = 1;
                    obj.pmMChannelName{i}{j} = uicontrol(obj.pALoops(i), ...
                        'Style', 'popupmenu', ...
                        'String', 'Channel Names', ...
                        'Value', 1, ...
                        'Units', 'pixels', ...
                        'Enable', enableStr, ...
                        'Position', [255 channelHeight  120 20], ...
                        'BackgroundColor', 'white', ...
                        'HorizontalAlignment', 'left', ...
                        'FontSize', obj.fontSize);  
%                     obj.MChannelName{i}{j} = 1;
                    obj.pmMRumpNo{i}{j} = uicontrol(obj.pALoops(i), ...
                        'Style', 'popupmenu', ...
                        'String', 'Rump Num', ...
                        'Value', 1, ...
                        'Units', 'pixels', ...
                        'Enable', enableStr, ...
                        'Position', [380 channelHeight  50 20], ...
                        'BackgroundColor', 'white', ...
                        'HorizontalAlignment', 'left', ...
                        'FontSize', obj.fontSize); 
%                     obj.MRumpNo{i}{j} = 1;
                    obj.etMVector{i}{j} = uicontrol(obj.pALoops(i), ...
                        'Style', 'edit', ...
                        'String', '0', ...
                        'Units', 'pixels', ...
                        'Enable', enableStr, ...
                        'Position', [435 channelHeight-8 195 obj.componentHeight+6], ...
                        'BackgroundColor', 'white', ...
                        'HorizontalAlignment', 'left', ...
                        'FontSize', obj.fontSize);
%                     obj.MVector{i}{j} = [];
                    obj.sMNoElements{i}{j} = uicontrol(obj.pALoops(i), ...
                        'Style', 'text', ...
                        'String', '0 elements', ...
                        'Units', 'pixels', ...
                        'Enable', enableStr, ...
                        'Position', [630 channelHeight-5 85 obj.componentHeight], ...
                        'BackgroundColor', [0.8 0.8 0.8], ...
                        'HorizontalAlignment', 'left', ...
                        'FontSize', obj.fontSize);   
%                     obj.MNoElements{i}{j} = '0 elements';

                    set(obj.pmMChangeTypes{i}{j}, 'Callback', {@obj.pmMChangeTypes_Callback});
                    set(obj.pmMEventName{i}{j}, 'Callback', {@obj.pmMEventName_Callback});
                    set(obj.pmMChannelName{i}{j}, 'Callback', {@obj.pmMChannelName_Callback});  
                    set(obj.etMVector{i}{j}, 'Callback', {@obj.etMVector_Callback});
                    set(obj.pmMRumpNo{i}{j}, 'Callback', {@obj.pmMRumpNo_Callback}); 
                end
                
            set(obj.pmASaveParam(i), 'Callback', {@obj.pmASaveParam_Callback});     
            set(obj.etASaveParamVals(i), 'Callback', {@obj.etASaveParamVals_Callback}); 
            set(obj.pbAAddChannel(i), 'Callback', {@obj.pbAAddChannel_Callback});
            set(obj.pbARemoveLastChannel(i), 'Callback', {@obj.pbARemoveLastChannel_Callback});
            
            end
            
            set(obj.pbOpenReadDir, 'Callback', {@obj.pbOpenReadDir_Callback});             
            set(obj.etReadDir, 'Callback', {@obj.etReadDir_Callback});          
            set(obj.etNumIterations, 'Callback', {@obj.etNumIterations_Callback});
            set(obj.pmIterationsOrder, 'Callback', {@obj.pmIterationsOrder_Callback});
            set(obj.pbOK, 'Callback', {@obj.pbOK_Callback});
            set(obj.pbCancel, 'Callback', {@obj.pbCancel_Callback});
            set(obj.pbSave, 'Callback', {@obj.pbSave_Callback});
            set(obj.pbLoad, 'Callback', {@obj.pbLoad_Callback});
            set(obj.tbEnableEdit, 'Callback', {@obj.tbEnableEdit_Callback});
            set(obj.pbAddLoop, 'Callback', {@obj.pbAddLoop_Callback});            
            set(obj.pbRemoveLastLoop, 'Callback', {@obj.pbRemoveLastLoop_Callback});     
%             set(obj.pmASaveParam, 'Callback', {@obj.pmASaveParam_Callback});     
%             set(obj.etASaveParamVals, 'Callback', {@obj.etASaveParamVals_Callback}); 
            
            minHeight = max(3*(2*(obj.componentHeight+5)+(2+4)*5), obj.pALoopsStartHeight); %min height of 3 panels
            set(obj.win, 'Position', [500 100 730 70+minHeight], 'MenuBar', 'None','Toolbar', 'None');
            guidata(obj.win, obj);
%             if obj.noIterations ~= -1 
                obj.updateGUI();
%             end
            set(obj.win, 'Visible', 'on');
          
            if  doWait == 1 && ishandle(obj.win)

                % Go into uiwait if the figure handle is still valid.
                % This is mostly the case during regular use.
                
%                 tmpObj = obj;                
%                 st = dbstack;

%                 uiwait(obj.win);
 
                obj.fPause = 1;
                guidata(obj.win, obj);
                while obj.fPause == 1
                    pause(1);
                    if  ishandle(obj.win)
                        obj = guidata(obj.win);
                    else
                        obj.fPause = 0;
                    end
                end
                    
                
                
                if  ishandle(obj.win)
                    obj = guidata(obj.win);
                    close(obj.win);
                    
                    if obj.noIterations ~= -1
                        obj.runs = 1:obj.noMeasurements;
                        obj.runs = obj.iterateVec(obj.runs, obj.iterationsOrder);
                        
                        obj.makeFolders(obj.baseBaseFolder, obj.noLoops);
                        
                    end
                else
                    obj.noIterations = -1;
                end
%                 if obj.noIterations == -1 && strcmp(st(2).name, 'imaging/lbMeasurementsList_KeyPressFcn')
%                     obj = tmpObj;
%                 end
                obj.win = -1;
            end
                       
        end
        
        function makeFolders(obj, baseFolder, loopInd)
            if loopInd == 1
                return
            end
            for i = 1 : length(obj.ASaveParamVals{loopInd})
                folders{i} = [baseFolder '\' obj.ASaveParamsStr{loopInd}{obj.ASaveParam(loopInd)} ' - ' ...
                    num2str(obj.ASaveParamVals{loopInd}(i))]; %#ok<AGROW>
            end
            for i = 1 : length(obj.ASaveParamVals{loopInd})
                [s,mess,messid] = mkdir(folders{i}); %#ok<NASGU,ASGLU>
                obj.makeFolders(folders{i}, loopInd-1)
            end
        end
        
        function obj = edit(obj, appData)
%             obj = guidata(obj.win);
            
            tmpObj = obj;
            obj.enableEdit = 0;
            obj = obj.initialize(appData, 1);
            if obj.noIterations == -1
                obj = tmpObj;
            end
            
%             guidata(obj.win, obj);
        end
                
        function [obj, newLVData] = next(obj, appData) 
%             obj = guidata(obj.win);
            if obj.position == obj.noMeasurements
                newLVData = [];
                return
            end
            
            obj.position = obj.position + 1;
            
            runInd = obj.runs(obj.position)-1;
            ind = zeros(1, obj.noLoops);
            for i = 1 : obj.noLoops
                totMul = 1;
                for j = 1 : obj.noLoops-i
                    totMul = totMul * length(obj.MVector{j}{1});
                end
                ind(i) = floor(runInd/totMul);
                runInd = rem(runInd, totMul);
            end
            ind = fliplr(ind)+1;
            
            newLVData = obj.LVData;
            obj.baseFolder = obj.baseBaseFolder;
            for i = obj.noLoops : -1 : 1
                if i > 1
                    obj.baseFolder = [obj.baseFolder '\' obj.ASaveParamsStr{i}{obj.ASaveParam(i)} ' - ' ...
                        num2str(obj.ASaveParamVals{i}(ind(i)))];
                end
                for j = 1 : obj.noChannels(i)
                    switch obj.MChangeTypes{i}{j}
                        case obj.valueTypes.eventStartTime
                            type = obj.valueTypes.eventStartTime;
                            event = obj.MEventName{i}{j};
                            channel = -1;
                            ramp = -1;
                            value = obj.MVector{i}{j}(ind(i));
                        case obj.valueTypes.analogRampTime
                            type = obj.valueTypes.analogRampTime;
                            event = obj.MEventName{i}{j};
                            channel = obj.MChannelName{i}{j};
                            ramp = obj.MRumpNo{i}{j};
                            value = obj.MVector{i}{j}(ind(i));
                        case obj.valueTypes.analogStartTime
                            type = obj.valueTypes.analogStartTime;
                            event = obj.MEventName{i}{j};
                            channel = obj.MChannelName{i}{j};
                            ramp = obj.MRumpNo{i}{j};
                            value = obj.MVector{i}{j}(ind(i));
                        case obj.valueTypes.analogValue
                            type = obj.valueTypes.analogValue;
                            event = obj.MEventName{i}{j};
                            channel = obj.MChannelName{i}{j};
                            ramp = obj.MRumpNo{i}{j};
                            value = obj.MVector{i}{j}(ind(i));
                        case obj.valueTypes.digitalTime
                            type = obj.valueTypes.digitalTime;
                            event = obj.MEventName{i}{j};
                            channel = obj.MChannelName{i}{j};
                            ramp = obj.MRumpNo{i}{j};
                            value = obj.MVector{i}{j}(ind(i));
                        case obj.valueTypes.digitalValue
                            type = obj.valueTypes.digitalValue;
                            event = obj.MEventName{i}{j};
                            channel = obj.MChannelName{i}{j};
                            ramp = obj.MRumpNo{i}{j};
                            value = obj.MVector{i}{j}(ind(i));
                        case obj.valueTypes.RFTime
                            type = obj.valueTypes.RFTime;
                            event = -1;
                            channel = obj.MChannelName{i}{j};
                            ramp = obj.MRumpNo{i}{j};
                            value = obj.MVector{i}{j}(ind(i));                        
                        case obj.valueTypes.RFCommand
                            type = obj.valueTypes.RFCommand;
                            event = -1;
                            channel = obj.MChannelName{i}{j};
                            ramp = obj.MRumpNo{i}{j};
                            value = [obj.MVectorRFCommands{i}{j} num2str(obj.MVector{i}{j}(ind(i)))];
                    end
                    newLVData = newLVData.changeValue(type, event, channel, ramp, value );
                end
            end

            % change save param, param val and folder
            appData.save.saveOtherParamStr = obj.ASaveOtherParamStr{1};
            set(appData.ui.pmSaveParam, 'String', obj.ASaveParamsStr{1});
            set(appData.ui.pmSaveParam, 'Value', obj.ASaveParam(1));
            set(appData.ui.etParamVal, 'String', num2str(obj.ASaveParamVals{1}(ind(1)))); 
%             set(appData.ui.etSaveDir, 'String', folder);
            

%             guidata(obj.win, obj);
        end
        
        function str = getMeasStr(obj, appData)
%             obj = guidata(obj.win);
%             str = 'aa';
            str = [appData.consts.availableLoops.str{appData.consts.availableLoops.generalLoop}  ' (' obj.getCurrMeas() ')'];
%             guidata(obj.win, obj);
        end
    
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        % callbacks
        function pbOpenReadDir_Callback(obj, object, eventdata)   %#ok<INUSL>
            dirName = uigetdir(get(obj.etReadDir, 'String'));
            if ( dirName ~= 0 )
                set(obj.etReadDir, 'String', dirName);
                obj.etReadDir_Callback(obj.etReadDir, eventdata);
            end
        end    
        
        function etReadDir_Callback(obj, object, eventdata)
            obj = guidata(obj.win);
            obj.baseBaseFolder = get(object, 'String');
            guidata(obj.win, obj);
        end    
        
        function etNumIterations_Callback(obj, object, eventdata)
            obj = guidata(obj.win);
            val = str2double(get(object, 'String'));
            if ( isnan(val) || val <= 0 || floor(val) ~= val )
                set(object, 'String', '1');
                errordlg('Input must be positive integer','Error', 'modal');
            else
                obj.noIterations = val;
            end
            guidata(obj.win, obj);
        end
        
        function pmIterationsOrder_Callback(obj, object, eventdata) 
            obj = guidata(obj.win);
            obj.iterationsOrder = get(object, 'Value');
            guidata(obj.win, obj);
        end           
        
        function pmASaveParam_Callback(obj, object, eventdata) 
            obj = guidata(obj.win);
            i = obj.getPIndex(object, obj.pmASaveParam);
            obj.ASaveParam(i) = get(object, 'Value');%obj.pmASaveParam, 'Value');
            obj.ASaveParamsStr{i} = get(object, 'String');
            if ( obj.ASaveParam(i) == length(obj.ASaveParamsStr{i}) )
                param = inputdlg('Enter param name:', 'Other param input');
                if isempty(param)
                    return
                end
                obj.ASaveOtherParamStr(i) = param;
                obj.ASaveParamsStr{i} = {['O.P. - ' param{1}] obj.ASaveParamsStr{i}{2:end}};
                set(object, 'String', obj.ASaveParamsStr{i});
                obj.ASaveParam(i) = 1;
                set(object, 'Value', 1);
            end
            guidata(obj.win, obj);
        end  
        
        function etASaveParamVals_Callback(obj, object, eventdata) 
            obj = guidata(obj.win);
            i = obj.getPIndex(object, obj.etASaveParamVals);
            try
                vec = eval(get(object, 'String'));
                obj.ASaveParamVals{i} = vec;
                obj.ASaveParamValsStr{i} = get(object, 'String');
%                 obj.noMeasurements = length(vec);
                obj.ANoElements{i} = [num2str(length(vec)) ' elements'];
                set(obj.sANoElements(i), 'String', obj.ANoElements{i});
                guidata(obj.win, obj);
            catch ex %#ok<NASGU>
                set(object, 'String', '');
                errordlg('Input must be an array','Error', 'modal');
            end
        end    
        
        function pbSave_Callback(obj, object, eventdata) 
%             obj = guidata(obj.win);
            obj = obj.updateObj();
%             guidata(obj.win, obj);
%              [file path] = uigetfile([obj.baseFolder '\*.mat']);
            [file path] = uiputfile([obj.saveFolder '\*.mat']);
            if file == 0
                return
            end
            o = obj; %#ok<NASGU>
%             o = obj.createSaveObj();             %#ok<NASGU>
            save([path file], 'o');
        end
        
        function pbLoad_Callback(obj, object, eventdata) 
            obj = guidata(obj.win);
%              [file path] = uigetfile([obj.baseFolder '\*.mat']);
             [file path] = uigetfile([obj.saveFolder '\*.mat']);
             if file == 0
                 return
             end
%             o = obj.createLoadObj([path file]);
            load([path file]);
%             o.tmpObj = tmpO.o;
            
            o.win = obj.win;
            o.baseBaseFolder = obj.baseBaseFolder;
            obj = o;
            guidata(obj.win, obj);
            obj.clearFigure(obj.win);
%             obj = obj.initialize([], 0);
            obj.updateGUI(o);
%             guidata(obj.win, obj);
        end     
        
        function pbOK_Callback(obj, object, eventdata)
            obj = obj.updateObj();
            %             obj = guidata(obj.win);
            
            totElements = 1;
            for i = 1 : obj.noLoops
                for j = 1 : obj.noChannels(i)
%                     if length(obj.saveParamVals) ~= length(obj.MVector{i}{j})
                    if length(obj.ASaveParamVals{i}) ~= length(obj.MVector{i}{j})
                        errordlg('All vectors mast be at the same length','Error', 'modal');
                        return
                    end
                end
%                 if length(obj.MVector{1}{1}) ~= length(obj.MVector{i}{1})
%                     errordlg('All vectors mast be at the same length','Error', 'modal');
%                     return
%                 end
                totElements = totElements * length(obj.ASaveParamVals{i});
            end
%             if length(obj.saveParamVals) ~= totElements
%                 errordlg('All vectors mast be at the same length','Error', 'modal');
%                 return
%             end
            obj.noMeasurements = totElements*obj.noIterations;
            
%             obj.noIterations = -2;
%             guidata(obj.win, obj);
%             uiresume(obj.win);     
            obj.fPause = 0;
            guidata(obj.win, obj);
%             delete(obj.win);
%             obj.win = -1;    
%             o = obj; %#ok<NASGU>
%             save(obj.tmpSaveName, 'o');
        end
%         function closeWin_Callback(obj, object, eventdata)
%             obj.noIterations = -1;
%             guidata(obj.win, obj);
%         end

        function pbCancel_Callback(obj, object, eventdata)
            obj = guidata(obj.win);
            obj.noIterations = -1; 
            obj.fPause = 0;
            guidata(obj.win, obj);
%             uiresume(obj.win); 
%             delete(obj.win);
%             obj.win = -1;  
%             obj.tmpObj.noIterations = -1;
%             obj.tmpObj.win = -1;  
%             o = obj.tmpObj;
%             if isempty(o)
%                 o = obj; %#ok<NASGU>
%             end
%             save(obj.tmpSaveName, 'o');
        end       
        
        function tbEnableEdit_Callback(obj, object, eventdata) 
            obj = guidata(obj.win);
            obj.enableEdit = get(object, 'Value');
%             obj.clearFigure(obj.win);
            guidata(obj.win, obj);
            obj = obj.initialize([], 0); %#ok<NASGU>
        end     
        
        function pbAddLoop_Callback(obj, object, eventdata) 
            obj = guidata(obj.win);
            obj = obj.updateObj();
            obj.noLoops = obj.noLoops + 1;
            obj.noChannels = [obj.noChannels 1];
            guidata(obj.win, obj);
            obj.clearFigure(obj.win);
%             guidata(obj.win, obj);
%             obj = obj.initialize([], 0); %#ok<NASGU>
        end      
        
        function pbRemoveLastLoop_Callback(obj, object, eventdata)
            obj = guidata(obj.win);
            if obj.noLoops == 1
                return
            end
            obj.noLoops = obj.noLoops - 1;            
            obj.noChannels = obj.noChannels(1:end-1);
            guidata(obj.win, obj);
%             obj = obj.updateObj();
%             guidata(obj.win, obj);
            obj.clearFigure(obj.win);
%             obj = obj.initialize([], 0); %#ok<NASGU>
        end     
        
        function pbAAddChannel_Callback(obj,  object, eventdata)
            obj = guidata(obj.win);
            pIndex = obj.getPIndex(object, obj.pbAAddChannel);
%             for i = 1 : length(obj.pbAAddChannel)
%                 if object == obj.pbAAddChannel(i)
%                     pIndex = i;
%                     break;
%                 end
%             end
            obj = obj.updateObj();
            obj.noChannels(pIndex) = obj.noChannels(pIndex) + 1;
            guidata(obj.win, obj);
            obj.clearFigure(obj.win);
%             obj = obj.initialize([], 0); %#ok<NASGU>
        end    
        
        function pbARemoveLastChannel_Callback(obj,  object, eventdata) %#ok<*INUSD>
            obj = guidata(obj.win);
            pIndex = obj.getPIndex(object, obj.pbARemoveLastChannel);
%             for i = 1 : length(obj.pbARemoveLastChannel)
%                 if object == obj.pbARemoveLastChannel(i)
%                     pIndex = i;
%                     break;
%                 end
%             end
            if obj.noChannels(pIndex) == 1
                return
            end
            obj.noChannels(pIndex) = obj.noChannels(pIndex) - 1;
            guidata(obj.win, obj);
            obj = obj.updateObj();
            guidata(obj.win, obj);
            obj.clearFigure(obj.win);
%             obj = obj.initialize([], 0);
        end
        
        function pmMChangeTypes_Callback(obj,  object, eventdata)
            obj = guidata(obj.win);
            [i j] = getPCIndex(obj, object, obj.pmMChangeTypes);
            obj.MChangeTypes{i}{j} = get(object, 'Value');
            
            switch obj.MChangeTypes{i}{j}
                case obj.valueTypes.eventStartTime
                    set(obj.pmMEventName{i}{j}, 'String', obj.eventNames);
                    set(obj.pmMChannelName{i}{j}, 'String', 'Channel Names');
                    set(obj.pmMRumpNo{i}{j}, 'String', 'Rump Num');
                    set(obj.pmMEventName{i}{j}, 'Value', 1);
                    set(obj.pmMChannelName{i}{j}, 'Value', 1);
                    set(obj.pmMRumpNo{i}{j}, 'Value', 1);
                    set(obj.etMVector{i}{j}, 'String', 'Values');
                    set(obj.sMNoElements{i}{j}, 'String','0 elements');
%                    obj.MEventName{i}{j} = 1;
%                    obj.MEventName{i}{j} = 1;
%                    obj.MRumpNo{i}{j} = 1;                 
%                     obj.MVector{i}{j} = '';
%                     obj.MVector{i}{j} = 0;               
               case {obj.valueTypes.analogRampTime obj.valueTypes.analogStartTime obj.valueTypes.analogValue}                   
                   set(obj.pmMEventName{i}{j}, 'String', obj.eventNames);                 
                   set(obj.pmMChannelName{i}{j}, 'String', obj.analogNames);                  
                   set(obj.pmMEventName{i}{j}, 'Value', 1);                 
                   set(obj.pmMChannelName{i}{j}, 'Value', 1);
                   str = cell(1, obj.LVData.getRampNum(get(obj.pmMChangeTypes{i}{j}, 'Value'), 1, 1));
                   for k = 1 : length(str)
                       str{k} = num2str(k);
                   end
                   if isempty(str)
                       set(obj.pmMRumpNo{i}{j}, 'String', 'Rump Num');
                   else
                       set(obj.pmMRumpNo{i}{j}, 'String', str);
                   end
                   set(obj.pmMRumpNo{i}{j}, 'Value', 1);
                   set(obj.etMVector{i}{j}, 'String', 'Values');
                   set(obj.sMNoElements{i}{j}, 'String','0 elements');
%                    obj.MEventName{i}{j} = 1;
%                    obj.MEventName{i}{j} = 1;
%                    obj.MRumpNo{i}{j} = 1;                 
               case {obj.valueTypes.digitalTime obj.valueTypes.digitalValue }                 
                   set(obj.pmMEventName{i}{j}, 'String', obj.eventNames);                 
                   set(obj.pmMChannelName{i}{j}, 'String', obj.digitalNames);
                   set(obj.pmMEventName{i}{j}, 'Value', 1);
                   set(obj.pmMChannelName{i}{j}, 'Value', 1);
                   str = cell(1, obj.LVData.getRampNum(get(obj.pmMChangeTypes{i}{j}, 'Value'), 1, 1));
                   for k = 1 : length(str)
                       str{k} = num2str(k);
                   end
                   if isempty(str)
                       set(obj.pmMRumpNo{i}{j}, 'String', 'Rump Num');
                   else
                       set(obj.pmMRumpNo{i}{j}, 'String', str);
                   end
                   set(obj.pmMRumpNo{i}{j}, 'Value', 1);
                   set(obj.etMVector{i}{j}, 'String', 'Values');
                   set(obj.sMNoElements{i}{j}, 'String','0 elements');
%                    obj.MEventName{i}{j} = 1;
%                    obj.MEventName{i}{j} = 1;
%                    obj.MRumpNo{i}{j} = 1;                 
                case {obj.valueTypes.RFTime obj.valueTypes.RFCommand}
                    set(obj.pmMEventName{i}{j}, 'String', 'Events Names');
                    set(obj.pmMChannelName{i}{j}, 'String', obj.RFNames);
                    set(obj.pmMEventName{i}{j}, 'Value', 1);
                    set(obj.pmMChannelName{i}{j}, 'Value', 1);
                    str = cell(1, obj.LVData.getRampNum(get(obj.pmMChangeTypes{i}{j}, 'Value'), -1, -1));
                    for k = 1 : length(str)
                        str{k} = num2str(k);
                    end
                    if isempty(str)
                        set(obj.pmMRumpNo{i}{j}, 'String', 'Rump Num');
                    else
                        set(obj.pmMRumpNo{i}{j}, 'String', str);
                    end
                    set(obj.pmMRumpNo{i}{j}, 'Value', 1);
                    set(obj.etMVector{i}{j}, 'String', 'Values');
                    set(obj.sMNoElements{i}{j}, 'String','0 elements');
                    %                    obj.MEventName{i}{j} = 1;
                    %                    obj.MEventName{i}{j} = 1;
                    %                    obj.MRumpNo{i}{j} = 1;
            end
%             obj.etMVector_Callback(obj.etMVector{i}{j}, [])
            % get event, chanel and ramp data
            obj.MEventName{i}{j} = get(obj.pmMEventName{i}{j}, 'Value');
            obj.MChannelName{i}{j} = get(obj.pmMChannelName{i}{j}, 'Value');
            obj.MRumpNo{i}{j} = get(obj.pmMRumpNo{i}{j}, 'Value');
            guidata(obj.win, obj);
        end
        
        function pmMEventName_Callback(obj,  object, eventdata)
            obj = guidata(obj.win);
            [i j] = getPCIndex(obj, object, obj.pmMEventName);
            obj.MEventName{i}{j} = get(obj.pmMEventName{i}{j}, 'Value');
            
%             str = cell(1, obj.LVData.getRampNum(obj.MChangeTypes{i}{j}, ...
%                 obj.MEventName{i}{j}, obj.MChannelName{i}{j}));            
            str = cell(1, obj.LVData.getRampNum(get(obj.pmMChangeTypes{i}{j}, 'Value'), ...
                get(obj.pmMEventName{i}{j}, 'Value'), get(obj.pmMChannelName{i}{j}, 'Value')));
            for k = 1 : length(str)
                str{k} = num2str(k);
            end
            if isempty(str)
                set(obj.pmMRumpNo{i}{j}, 'String', 'Rump Num');
            else
                set(obj.pmMRumpNo{i}{j}, 'String', str);
            end
            guidata(obj.win, obj);
        end
        
        function pmMChannelName_Callback(obj,  object, eventdata)
            obj = guidata(obj.win);
            [i j] = getPCIndex(obj, object, obj.pmMChannelName);
            obj.MChannelName{i}{j} = get(obj.pmMChannelName{i}{j}, 'Value');
            
%             str = cell(1, obj.LVData.getRampNum(obj.MChangeTypes{i}{j}, ...
%                 obj.MEventName{i}{j}, obj.MChannelName{i}{j}));
            
            str = cell(1, obj.LVData.getRampNum(get(obj.pmMChangeTypes{i}{j}, 'Value'), ...
                get(obj.pmMEventName{i}{j}, 'Value'), get(obj.pmMChannelName{i}{j}, 'Value')));
            for k = 1 : length(str)
                str{k} = num2str(k);
            end
            if isempty(str)
                set(obj.pmMRumpNo{i}{j}, 'String', 'Rump Num');
            else
                set(obj.pmMRumpNo{i}{j}, 'String', str);
            end
            guidata(obj.win, obj);
        end
        
        function pmMRumpNo_Callback(obj,  object, eventdata)
            obj = guidata(obj.win);
            [i j] = getPCIndex(obj, object, obj.pmMRumpNo);
            obj.MRumpNo{i}{j} = get(object, 'Value');
            guidata(obj.win, obj);
        end
        
        function etMVector_Callback(obj,  object, eventdata)
            obj = guidata(obj.win);
            [i j] = getPCIndex(obj, object, obj.etMVector);
            [obj.MVectorRFCommands{i}{j} obj.MVector{i}{j}] = obj.evalStr(object, get(object, 'String'), i, j);
            obj.MVectorStr{i}{j} = get(object, 'String');
%             try
%                 vec = eval(get(object, 'String'));
%             catch ex
%                 if get(obj.pmMChangeTypes{i}{j}, 'Value') == obj.valueTypes.RFCommand
%                     s=regexp( get(object, 'String'), ' ', 'split');
%                     flag = 0;
%                     for k = 1 : length(s)
%                         try
%                             vec = eval(s{k});
%                             flag = 1;
%                             break;
%                         catch ex
%                         end
%                     end
%                     if flag == 0
%                         set(object, 'String', '');
%                         set(obj.sMNoElements{i}{j}, 'String','0 elements');
%                         errordlg('Input must be an array','Error', 'modal');
%                         return
%                     end
%                 else
%                     set(object, 'String', '');
%                     set(obj.sMNoElements{i}{j}, 'String','0 elements');
%                     errordlg('Input must be an array','Error', 'modal');
%                     return
%                 end
%             end
%             obj.MVector{i}{j} = vec;
%             obj.MVectorRFCommands{i}{j} = RFCommand;
            obj.MNoElements{i}{j} = [num2str(length(obj.MVector{i}{j})) ' elements'];
            set(obj.sMNoElements{i}{j}, 'String', obj.MNoElements{i}{j});
            guidata(obj.win, obj);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = updateObj(obj) %createSaveObj(obj)
            %             o = obj;
            obj = guidata(obj.win);
            
            % clear obj
            obj.ASaveParam = [];
            obj.ASaveParamsStr = {};
            obj.ASaveParamVals = {};
            obj.ANoElements = {};
            obj.MChangeTypes = {};
            obj.MEventName = {};
            obj.MChannelName = {};
            obj.MRumpNo = {};
            obj.MVectorRFCommands = {};
            obj.MVector = {};
            obj.MNoElements = {};
            
            obj.baseBaseFolder =  get(obj.etReadDir, 'String');
            obj.noIterations = str2double(get(obj.etNumIterations, 'String'));
            obj.enableEdit = get(obj.tbEnableEdit, 'Value');
            obj.iterationsOrder = get(obj.pmIterationsOrder, 'Value');
%             obj.saveParam = get(obj.pmASaveParam, 'Value');
%             obj.saveParamVals = eval(get(obj.etASaveParamVals, 'String'));
%             obj.noElements = get(obj.sANoElements, 'String');
            for i = 1 : obj.noLoops
                obj.ASaveParam(i) = get(obj.pmASaveParam(i), 'Value');
                obj.ASaveParamsStr{i} = get(obj.pmASaveParam(i), 'String');
                obj.ASaveParamVals{i} = eval(get(obj.etASaveParamVals(i), 'String'));
                obj.ASaveParamValsStr{i} = get(obj.etASaveParamVals(i), 'String');
                obj.ANoElements{i} = get(obj.sANoElements(i), 'String');
                for j = 1 : obj.noChannels(i)
                    obj.MChangeTypes{i}{j} = get(obj.pmMChangeTypes{i}{j}, 'Value'); %2D cell array
                    obj.MEventName{i}{j} = get(obj.pmMEventName{i}{j}, 'Value');
                    obj.MChannelName{i}{j} = get(obj.pmMChannelName{i}{j}, 'Value');
                    obj.MRumpNo{i}{j} = get(obj.pmMRumpNo{i}{j}, 'Value');
%                     obj.MVector{i}{j} = obj.evalStr(get(obj.etMVector{i}{j}, 'String'), i, j);
                    [obj.MVectorRFCommands{i}{j} obj.MVector{i}{j}] = obj.evalStr(obj.etMVector{i}{j}, get(obj.etMVector{i}{j}, 'String'), i, j);
                    obj.MVectorStr{i}{j} = get(obj.etMVector{i}{j}, 'String');
                    obj.MNoElements{i}{j} = get(obj.sMNoElements{i}{j}, 'String');
                end
            end
            guidata(obj.win, obj);
        end
        
            
        function updateGUI(obj, o)
            obj = guidata(obj.win);
            if nargin < 2
                o = obj;
            end
            set(obj.etReadDir, 'String', num2str(o.baseBaseFolder));
            set(obj.etNumIterations, 'String', num2str(o.noIterations));
            set(obj.tbEnableEdit, 'Value', o.enableEdit);
            set(obj.pmIterationsOrder, 'Value', o.iterationsOrder);
%             set(obj.pmASaveParam, 'Value', o.saveParam);
%             set(obj.etASaveParamVals, 'String', vec2str(o.saveParamVals));
%             set(obj.sANoElements, 'String', o.noElements);
            for i = 1 : length(o.MChangeTypes)%o.noLoops
                set(obj.pmASaveParam(i), 'Value', o.ASaveParam(i)); 
                set(obj.pmASaveParam(i), 'String', o.ASaveParamsStr{i});
                set(obj.etASaveParamVals(i), 'String', o.ASaveParamValsStr{i});%vec2str(o.ASaveParamVals{i}));
                set(obj.sANoElements(i), 'String', o.ANoElements{i});
                for j = 1 : length(o.MChangeTypes{i})%o.noChannels(i)
                    set(obj.pmMChangeTypes{i}{j}, 'Value', o.MChangeTypes{i}{j}); %2D cell array
                    obj.pmMChangeTypes_Callback( obj.pmMChangeTypes{i}{j}, []);
                    set(obj.pmMEventName{i}{j}, 'Value', o.MEventName{i}{j});
                    set(obj.pmMChannelName{i}{j}, 'Value', o.MChannelName{i}{j});
                    obj.pmMChannelName_Callback( obj.pmMChannelName{i}{j},  []);
                    set(obj.pmMRumpNo{i}{j}, 'Value', o.MRumpNo{i}{j});
                    
                    set(obj.etMVector{i}{j}, 'String', o.MVectorStr{i}{j});%[obj.MVectorRFCommands{i}{j} o.MVectorStr{i}{j}]);%vec2str(o.MVector{i}{j})]);
                    set(obj.sMNoElements{i}{j}, 'String', o.MNoElements{i}{j});
                end
            end
            guidata(obj.win, obj);
        end
        
        function clearFigure(obj, fig)
%             obj = guidata(obj.win);
            obj.pALoops = [];
            obj.pmASaveParam = [];
            obj.etASaveParamVals = [];
            obj.sANoElements = [];
            
            obj.pbAAddChannel = [];
            obj.pbARemoveLastChannel = [];
            obj.pmMChangeTypes = {}; %2D cell array
            obj.pmMEventName = {};
            obj.pmMChannelName = {};     
            obj.pmMRumpNo = {};
            obj.etMVector = {};
            obj.sMNoElements = {};
            
            clf(fig, 'reset');
            guidata(obj.win, obj);
            obj = obj.initialize([], 0); %#ok<NASGU>
        end
        
        function [pIndex cIndex] = getPCIndex(obj, object, objectM) %#ok<MANU>
            pIndex = 0;
            cIndex = 0;
            for i = 1 : length(objectM)
                for j = 1 : length(objectM{i})
                    if object == objectM{i}{j}
                        pIndex = i;
                        cIndex = j;
                        break;
                    end
                end
            end
        end
        
        function [pIndex] = getPIndex(obj, object, objectA) %#ok<MANU>
            pIndex = 0;
            for i = 1 : length(objectA)
                if object == objectA(i)
                    pIndex = i;
                    break;
                end
            end
        end
        
        function [RFCommand vec] = evalStr(obj, object, str, i, j)
            RFCommand = '';
            vec = [];
            try
                vec = eval(str);
            catch ex %#ok<NASGU>
                if i ~= 0 && get(obj.pmMChangeTypes{i}{j}, 'Value') == obj.valueTypes.RFCommand
                    s=regexp( str, ' ', 'split');
                    flag = 0;
                    for k = 1 : length(s)
                        try
                            vec = eval([s{k:end}]);
                            flag = 1;
                            break;
                        catch ex %#ok<NASGU>
                            RFCommand = [RFCommand s{k} ' ']; %#ok<AGROW>
                        end
                    end
                    if flag == 0
                        set(object, 'String', '');
                        set(obj.sMNoElements{i}{j}, 'String','0 elements');
                        errordlg('Input must be an array','Error', 'modal');
                        return
                    end
                else
                    set(object, 'String', '');
                    set(obj.sMNoElements{i}{j}, 'String','0 elements');
                    errordlg('Input must be an array','Error', 'modal');
                    return
                end
            end
        end
        
        
        
    end
    
    
    
end