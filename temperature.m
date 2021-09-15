
function [graph, Result] = temperature(totAppData, readDir, atomType)
if numel(totAppData) > 1
    len = length(totAppData);
    tof = zeros(1, len); %[sec]
    sx = zeros(1, len); %[meter]
    sy = zeros(1, len); %[meter]
    
    switch atomType
        case 1
            mass = totAppData{1}.consts.MNa;
        case 2
            mass = totAppData{1}.consts.MLi6;
    end
    
    fitType = totAppData{1}.data.fitType;
    
    for i = 1 : len
        %     if ( totAppData{i}.save.saveParam ~= totAppData{i}.consts.saveParams.TOF )
        %         warndlg({['The save parameter in data-' num2str( totAppData{i}.save.picNo) ' is not TOF (' num2str(totAppData{i}.consts.saveParams.TOF) ').']; ...
        %             ['It is:' num2str(totAppData{i}.save.saveParam) '.']} , 'Warning', 'modal');
        %     end
        tof(i) = totAppData{i}.save.saveParamVal;
        sx(i) = totAppData{i}.data.fits{ fitType }.xUnitSize;
        sy(i) = totAppData{i}.data.fits{ fitType }.yUnitSize;
    end
    
    tof = tof/(1000*1); % change to seconds
    sx = sx  * totAppData{i}.data.camera.xPixSz; % change to meters
    sy = sy  * totAppData{i}.data.camera.yPixSz; % change to meters
    
    tof2 = tof.^2;
    sx2 = sx.^2;
    sy2 = sy.^2;
    
    [resX gofX outX] = fit(tof2', sx2', 'poly1');
    if ( length(sx2) > 2 )
        confX = confint(resX);
        confX = (confX(2,:)-confX(1,:))/2;
    end
    [resY gofY outY] = fit(tof2', sy2', 'poly1');
    if ( length(sy2) > 2 )
        confY = confint(resY);
        confY = (confY(2,:)-confY(1,:))/2;
    end
    
    graph = figure( 'FileName', [readDir '_temp.fig']);
    plot( tof2, sx2, 'ob');
    hold on
    % plot(tof2, resX.p1*tof2+resX.p2, 'c');
    plot(resX, 'c');
    plot( tof2, sy2, 'or');
    % plot(tof2, resY.p1*tof2+resY.p2, 'm');
    plot(resY, 'm');
    hold off
    
    title('Temperature calculations');
    set(gca,'Ylabel',text('String', 'sigma^2 [m^2]'));
    set(gca,'Xlabel',text('String', 'TOF^2 [sec^2]'));
    set(gcf, 'Name', 'Temperature');
    legend('off')
    
    Tx = resX.p1 / (totAppData{1}.consts.Kb)*mass * 1e6; % we use: slope = sigma_v^2 = k_B*T/m => T = slope*m/k_B
    TxError = confX(1) / totAppData{1}.consts.Kb*mass * 1e6;
    Ty = resY.p1 / (totAppData{1}.consts.Kb)*mass * 1e6;
    TyError = confY(1) / totAppData{1}.consts.Kb*mass * 1e6;
    Result = [Tx Ty];
    
    %This piece of code saves x and y temperature data and their erros to main
    %workspace variables Tx, Txe, Ty, TyE
    %assign temperatures and errors to base workspace:
    % assignin('base','Txtemp',Tx)
    % assignin('base','TxEtemp',TxError)
    % assignin('base','Tytemp',Ty)
    % assignin('base','TyEtemp',TyError)
    %evaluate in main workspace: append data to existing vector; create vector if
    %it dosen't exist
    % evalin('base', ['if exist(''Tx'',''var'')' newline  'Tx=[Tx Txtemp]; TxE=[TxE TxEtemp]; Ty=[Ty Tytemp]; TyE=[TyE TyEtemp];' newline 'else' newline 'Tx=Txtemp;TxE=TxEtemp; Ty=Tytemp; TyE=TyEtemp;' newline 'end'])
    
    if numel(totAppData)>2
        text( 0.1, 0.5, {'fit function: ax + b', ...
            ['in x direction (blue): ' num2str(resX.p1) 'x + ' num2str(resX.p2) ', R^2 = ' num2str(gofX.rsquare) ] ...
            ['       T_x = ' num2str(Tx) ...
            ' +/- ' num2str(TxError) ' \muK'] ...
            ['       \sigma_{0x} = ' num2str(sqrt(resX.p2)*1000) ' +/- ' num2str(confX(2)*1000) ' mm' ] ...
            ['in y direction (red):  ' num2str(resY.p1) 'x + ' num2str(resY.p2) ', R^2 = ' num2str(gofY.rsquare) ] ...
            ['       T_y = ' num2str(Ty)  ...
            ' +/- ' num2str(TyError) ' \muK'] ...
            ['       \sigma_{0y} = ' num2str(sqrt(resY.p2)*1000) ' +/- ' num2str(confX(2)*1000) ' mm' ] },...
            'Units', 'Normalized'); %     ['average density = ' num2str(n*1e-6) '1/cm^3'] },...
    else
        text( 0.1, 0.5, {'fit function: ax + b', ...
            ['in x direction (blue): ' num2str(resX.p1) 'x + ' num2str(resX.p2) ', R^2 = ' num2str(gofX.rsquare) ] ...
            ['       T_x = ' num2str(TxError) ' \muK'] ...
            ['       \sigma_{0x} = ' num2str(sqrt(resX.p2)*1000)  ' mm' ] ...
            ['in y direction (red):  ' num2str(resY.p1) 'x + ' num2str(resY.p2) ', R^2 = ' num2str(gofY.rsquare) ] ...
            ['       T_y = ' num2str(TyError)  ' \muK'] ...
            ['       \sigma_{0y} = ' num2str(sqrt(resY.p2)*1000) ' mm' ] },...
            'Units', 'Normalized'); %     ['average density = ' num2str(n*1e-6) '1/cm^3'] },...
    end
    save( [readDir '_temp.fig'], 'Tx', 'TxError', 'Ty', 'TyError')
else
    graph =[];
    Result = [];
end


