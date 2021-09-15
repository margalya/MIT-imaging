function analyzeFoldersGUI( appData )
%add: control if parameter is spatial or not (to multiply by yPixSize)
%add: checkbox for checkGOF or not
height = 400;
width = 750;
h = figure('Visible', 'on', ...
    'Name', 'Analyze Folders', ...
    'Units', 'Pixels', ...
    'Position', [200 200 width height], ...
    'Resize', 'on', ...
    'MenuBar', 'None', ...
    'Toolbar', 'None');

Operator.str = {'mean +/- std', 'mean +/- sem', 'Norm. Vis.', 'MeanX +/- std'};
Operator.meanStd = 1;
Operator.meanSEM = 2;
Operator.normVis = 3;
Operator.meanXstd = 4;
Operator.value = 1;
pmOperator = uicontrol(h, ...
    'Style', 'popupmenu', ...
    'String', Operator.str, ...
    'Value', Operator.value, ...
    'Units', 'pixels', ...
    'Position', [90 height-30 120 27], ...
    'BackgroundColor', 'white', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', appData.consts.fontSize, ...
    'Callback', {@pmOperator_Callback});
% text: Operator
uicontrol(h, ...
    'Style', 'text', ...
    'String', 'Operator', ...
    'Units', 'pixels', ...
    'Position', [10 height-33 65 27], ...
    'BackgroundColor', [0.8 0.8 0.8], ...
    'HorizontalAlignment', 'left', ...
    'FontSize', appData.consts.fontSize);

ParamName.str = {'Visibility', 'lambda', 'sigmaY', 'phi', 'Y Position', 'mF1'};
ParamName.v = 1;
ParamName.lambda = 2;
ParamName.sigmaY = 3;
ParamName.phi = 4;
ParamName.yPos = 5;
ParamName.mF1 = 6;
ParamName.value = 1;
pmParamName = uicontrol(h, ...
    'Style', 'popupmenu', ...
    'String', ParamName.str, ...
    'Value', ParamName.value, ...
    'Units', 'pixels', ...
    'Position', [90 height-60 120 27], ...
    'BackgroundColor', 'white', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', appData.consts.fontSize, ...
    'Callback', {@pmParamName_Callback});
% text: Parameter
uicontrol(h, ...
    'Style', 'text', ...
    'String', 'Parameter', ...
    'Units', 'pixels', ...
    'Position', [10 height-63 70 27], ...
    'BackgroundColor', [0.8 0.8 0.8], ...
    'HorizontalAlignment', 'left', ...
    'FontSize', appData.consts.fontSize);

%push button analyze
uicontrol(h, ...
    'Style', 'pushbutton', ...
    'String', 'Analyze', ...
    'Units', 'pixels', ...
    'Position', [width-100 height-40 90 27], ...
    'BackgroundColor', [0.8 0.8 0.8], ...
    'FontSize', appData.consts.fontSize, ...
    'Callback', {@pbAnalyze_Callback});

%push button select folders
uicontrol(h, ...
    'Style', 'pushbutton', ...
    'String', 'Select Folders', ...
    'Units', 'pixels', ...
    'Position', [width-230 height-40 120 27], ...
    'BackgroundColor', [0.8 0.8 0.8], ...
    'FontSize', appData.consts.fontSize, ...
    'Callback', {@pbSelectFolders_Callback});

lbFoldersListDisplay = uicontrol(h, ...
    'Style', 'listbox', ...
    'String', '', ...
    'Units', 'pixels', ...
    'Position', [10 10 width-20 0.72*height], ...
    'BackgroundColor', 'white', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', appData.consts.fontSize);

%edit text x Values
etxValues = uicontrol(h, ...
    'Style', 'edit', ...
    'String', '', ...
    'Units', 'pixels', ...
    'Position', [90 height-90 300 27], ...
    'BackgroundColor', 'white', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', appData.consts.fontSize);
% text: x Values
uicontrol(h, ...
    'Style', 'text', ...
    'String', 'x Values', ...
    'Units', 'pixels', ...
    'Position', [10 height-93 70 27], ...
    'BackgroundColor', [0.8 0.8 0.8], ...
    'HorizontalAlignment', 'left', ...
    'FontSize', appData.consts.fontSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function pmParamName_Callback(object, eventdata) %#ok<INUSD>
        ParamName.value = get(object, 'Value');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function pmOperator_Callback(object, eventdata) %#ok<INUSD>
        Operator.value = get(object, 'Value');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function pbSelectFolders_Callback(object, eventdata) %#ok<INUSD>
        foldersList = uipickfiles('FilterSpec', appData.analyze.readDir);
        if isnumeric(foldersList) % selection is empty
            foldersList='';
        elseif length(foldersList)==1
            [~, ~, ext] = fileparts(foldersList{1});
            if strcmp(ext,'.m') %if input is m-file containing 'foldersList' variable
                xValues = [];
                xLabel = '';
                run(foldersList{1});
                if exist('xValues','var')
                    set(etxValues, 'String', strrep(num2str(xValues),'         ',' ') );
                end
            elseif strcmp(ext, '.mat')
                load(foldersList{1});
            end
        end
        set(lbFoldersListDisplay,'string',foldersList);
    end

    function pbAnalyze_Callback(object, eventdata) %#ok<INUSD>
        % loads fitting parameter 'paramName' (string) from data base, and return the mean param value and its SEM -
        % standard error of the mean. In case that the desired parameter is the
        % phase, the paramSEM is the minimum stndard deviation of the phase
        
        warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle')
        
        paramName = ParamName.str{ParamName.value};
        operator = Operator.str{Operator.value};
        LineSpec = appData.consts.pmLineSpec.str{appData.analyze.LineSpec};
        LineSpec = strrep(LineSpec,', mean',''); %remove ' mean' from LineSpec string
        
        foldersList = get(lbFoldersListDisplay,'string');
        if ~iscell(foldersList)
            return
        end
        xValues = get(etxValues,'string');
        if isempty(xValues)
            xValues = 1:length(foldersList);
        else
            xValues = str2num(xValues); %#ok<ST2NM>
        end
        
        if (length(xValues)~=length(foldersList))
            warndlg('x Values and folders list are not in the same length')
            return
        end
        
        % load('E:\Dropbox\BEC2\Arrow of time\Data analysis\Fig 4 optimized data\foldersListMax.mat');
        % foldersList = foldersListMax;
        % xValues = [4:2:16];
        % run('E:\Dropbox\BEC2\Arrow of time\Data analysis\Vn vs initial position\loadFoldersList_wo_50Hz.m');
        % xValues = 1:length(foldersList); %#ok<USENS>
        
        %initialize empty arrays
        param = zeros(1, length(foldersList) );
        paramStd = param;
        paramSEM = param;
        visPicMean = param;
        visPicMeanconf = param;
        paramX = [];
        paramXStd = [];
        paramXSEM = [];
        yPixSz = appData.consts.cameras{appData.options.cameraType}.yPixSz*1e6;
        chipStart = appData.data.camera.chipStart;
        
        %if parameter is the visibility of the average picture - then always read
        %visibility parameter
        if strcmp(operator,'Norm. Vis.');
            paramName = 'v';
        end
        progressbar(0, 0);
        for i = 1 : length(foldersList)
            
            fileList = dir([foldersList{i} '\*.mat']);
            fileList = fileList( ~arrayfun(@(x) strcmp(x.name,'data-1000.mat'),fileList) ); %remove data-1000.mat and 1001 from the list
            fileList = fileList( ~arrayfun(@(x) strcmp(x.name,'data-1001.mat'),fileList) );
            N = length(fileList);
            
            %initialize temporary parameter
            tempParam = zeros(1, N);
            sse = tempParam;
            rsquare = tempParam;
            rmse = tempParam;
            lambda = tempParam;
            saveParamVal = tempParam;
            %go through all files in folder and extract paramName from fit result
            for j = 1 : N
                load( [ foldersList{i} '\' fileList(j).name ], 'savedData');
                if strcmp(operator, 'MeanX +/- std') %not averaging param over entire folder
                    if strcmp(paramName,'mF1');
                        tempParam(j) = savedData.data.fits{ savedData.data.fitType }.mF1*100; %extract mF1, in [%]
                        saveParamVal(j) = savedData.save.saveParamVal;
                    end
                else % averaging param over entire folder
                    try
                        if strcmp(paramName,'Y Position'); %extract y position, taking into account chipstart, and also transform to mm
                            tempParam(j) = (savedData.data.fits{ savedData.data.fitType }.yCenter-savedData.data.camera.chipStart) ...
                                * savedData.data.camera.yPixSz *1000;
                        elseif strcmp(paramName,'lambda');
                            tempParam(j) = savedData.data.fits{ savedData.data.fitType }.res.(paramName) * savedData.data.camera.yPixSz *1000; %extract lambda, with proper scaling from pixels to mm
                        elseif strcmp(paramName,'sigmaY');
                            tempParam(j) = savedData.data.fits{ savedData.data.fitType }.res.w * savedData.data.camera.yPixSz *1000; %extract Gaussian width for fitfringesY, with proper scaling from pixels to mm
                        else
                            tempParam(j) = savedData.data.fits{ savedData.data.fitType }.res.(paramName);
                        end
                        sse(j) = savedData.data.fits{ savedData.data.fitType }.gof.sse;
                        rsquare(j) = savedData.data.fits{ savedData.data.fitType }.gof.rsquare;
                        rmse(j) = savedData.data.fits{ savedData.data.fitType }.gof.rmse;
                        lambda(j) = savedData.data.fits{ savedData.data.fitType }.res.('lambda');
                        %
                    catch
                        display([ 'error at: folder ' num2str(i) ', file ' num2str(j)] );
                    end
                end
                progressbar([], j/N);
            end
            
            checkGOF(sse, rsquare, rmse, lambda, i, fileList); %check for goodness of fit in current folder
            
            switch operator
                case 'mean +/- std' %mean over all single shots inside each folder
                    param(i) = mean(tempParam);
                    paramStd(i) = std(tempParam);
                    if strcmp(paramName,'phi');
                        paramStd(i) = minstd(tempParam); %for phase - calculate the minstd, overwrite previous value
                    end
                case 'mean +/- sem'
                    param(i) = mean(tempParam);
                    paramSEM(i) = std(tempParam) / sqrt(N);
                    if strcmp(paramName,'phi');
                        paramSEM(i) = minstd(tempParam) / sqrt(N); %for phase - calculate the minstd, overwrite previous value
                    end
                case 'Norm. Vis.'
                    param(i) = mean(tempParam); %visibility
                    paramSEM(i) = std(tempParam) / sqrt(N); %visibility SEM
                    
                    load( [ foldersList{i} '\data-1000.mat' ],'savedData');
                    visPicMean(i) = savedData.data.fits{ savedData.data.fitType }.res.v; %visPicMean
                    tmpconfint = confint( savedData.data.fits{ savedData.data.fitType }.res );
                    %         find index of the visibility parameter (changed in different versions of FitFringesY.m)
                    visIndx = find(strcmp(coeffnames( savedData.data.fits{ savedData.data.fitType }.res ),'v'));
                    visPicMeanconf(i) = (tmpconfint(2,visIndx)-tmpconfint(1,visIndx))/2; %visPicMeanconf
                case 'MeanX +/- std' %not averaging param over entire folder
                    clear x_unique yMean yStd
                    [ x_unique, yMean, yStd, ySEM ] = meanX( saveParamVal, tempParam);
                    paramX = [paramX yMean];
                    paramXStd = [paramXStd yStd];
                    paramXSEM = [paramXSEM ySEM];
                case 'none'
            end
            clear savedData
            
            progressbar(i/length(foldersList), 0);
        end
        
        figure( 'FileName', foldersList{1});
        switch operator
            case 'mean +/- std'
                %         errorbar(xValues, (param-chipStart) * appData.data.camera.yPixSz * 1000, paramStd * appData.data.camera.yPixSz * 1000, 'o'); %for Y pos - multiply by pixel size and subtract chipstart
                errorbar(xValues, param, paramStd, 'o'); %for regular
                %         xlabel('First kick time [\mus]');
            case 'mean +/- sem'
                errorbar(xValues, param, paramSEM, 'o');
                %         plot(xValues, exp(- (paramStd.^2) /2), LineSpec);
            case 'Norm. Vis.'
                Vn = visPicMean./param; %Normalized visibility
                VnErr = Vn .* sqrt( (visPicMeanconf./visPicMean).^2 + (paramSEM./param).^2);
                errorbar(xValues, Vn, VnErr, 'o');
                %         xlabel(xLabel);
                ylabel('Normalized visbility');
                legend('V_{Mean pic}/<V_{Single Shot}>');
            case 'MeanX +/- std'
                errorbar(1:length(paramX), paramX, paramXStd, 'o')
                xlabel('Point number')
                ylabel('param +/- std')
            case 'none'
        end
        
    end

    function checkGOF( sse, rsquare, rmse, lambda, i, fileList)
        % function checks if there is a suspicious image in which the fit might be
        % bad, by cheking the goodness of fit (GOF) using mean and std values
        
        if any( (sse-mean(sse)) > 3 * std(sse)) %alert on values > mean()+n*std()
            display( char( ['Outlier SSE value at: folder ' num2str(i) ', files '], char(fileList((sse-mean(sse)) > 3 * std(sse)).name)));
        end
        
        if any( rsquare < (mean(rsquare) - 3 * std(rsquare)) ) %alert on values < mean()+n*std()
            display( char( ['Outlier RSquare value at: folder ' num2str(i) ', files '], char(fileList( rsquare < (mean(rsquare) - 3 * std(rsquare)) ).name)));
        end
        
        if any( (rmse-mean(rmse)) > 3 * std(rmse)) %alert on values > mean()+n*std()
            display( char( ['Outlier RMSE value at: folder ' num2str(i) ', files '], char(fileList((rmse-mean(rmse)) > 3 * std(rmse)).name)));
        end
        
        if any( abs(lambda-mean(lambda)) > 3 * std(lambda)) %alert on both sides: mean()+/-n*std()
            display( char( ['Outlier wavelength value at: folder ' num2str(i) ', files '], char(fileList( abs(lambda-mean(lambda)) > 3 * std(lambda) ).name)));
        end
        
    end

end