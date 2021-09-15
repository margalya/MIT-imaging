
function [graph Result] = temperature(totAppData, readDir, atomType)
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

for ( i = 1 : len )
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

%average data with same TOF values
[ tof2_unique, sx2Mean, ~, sx2Sem ] = meanXAM( tof2, sx2);
[ ~, sy2Mean, ~, sy2Sem ] = meanXAM( tof2, sy2);

[resX, gofX, outX] = fit(tof2_unique', sx2Mean', 'poly1');  %#ok<ASGLU>
% if ( length(sx2) > 2 )
confX = confint(resX);
confX = (confX(2,:)-confX(1,:))/2;
% end
[resY, gofY, outY] = fit(tof2_unique', sy2Mean', 'poly1');  %#ok<ASGLU>
% if ( length(sy2) > 2 )
confY = confint(resY);
confY = (confY(2,:)-confY(1,:))/2;
% end

graph = figure( 'FileName', [readDir '_temp.fig']);
errorbar( tof2_unique, sx2Mean, sx2Sem, 'ob');
hold on
% plot(tof2, resX.p1*tof2+resX.p2, 'c');
plot(resX, 'c');
errorbar( tof2_unique, sy2Mean, sy2Sem, 'or');
% plot(tof2, resY.p1*tof2+resY.p2, 'm');
plot(resY, 'm');
hold off

title('Temperature calculations');
set(gca,'Ylabel',text('String', 'sigma^2 [m^2]'));
set(gca,'Xlabel',text('String', 'TOF^2 [sec^2]'));
set(gcf, 'Name', 'Temperature');
legend('off')

Tx = resX.p1 / (2*totAppData{1}.consts.Kb)*mass * 1e6; % we use: sigma_v = sqrt(2*k_B*T/m) => T = sigma_v^2*m/(2*k_B)
Ty = resY.p1 / (2*totAppData{1}.consts.Kb)*mass * 1e6;

Result = [Tx Ty];

%evaluate the density, elastic collision rate, and htree-body loss rate
% if all([sqrt(resX.p2) sqrt(resY.p2)])
%     n = Natoms./( mean([sqrt(resX.p2) sqrt(resY.p2)]) ); %density in 1/m^3
% end

% %This piece of code saves x and y temperature data and their erros to main
% %workspace variables Tx, Txe, Ty, TyE
% %assign temperatures and errors to base workspace:
% assignin('base','Txtemp',resX.p1 / totAppData{1}.consts.Kb*totAppData{1}.consts.MNa * 1e6) 
% assignin('base','TxEtemp',confX(1) / totAppData{1}.consts.Kb*totAppData{1}.consts.MNa * 1e6)
% assignin('base','Tytemp',resY.p1 / totAppData{1}.consts.Kb*totAppData{1}.consts.MNa * 1e6)
% assignin('base','TyEtemp',confY(1) / totAppData{1}.consts.Kb*totAppData{1}.consts.MNa * 1e6)
% %evaluate in main workspace: append data to existing vector; create vector if
% %it dosen't exist
% evalin('base', ['if exist(''Tx'',''var'')' sprintf('\n')  'Tx=[Tx Txtemp]; TxE=[TxE TxEtemp]; Ty=[Ty Tytemp]; TyE=[TyE TyEtemp];' sprintf('\n') 'else' sprintf('\n') 'Tx=Txtemp;TxE=TxEtemp; Ty=Tytemp; TyE=TyEtemp;' sprintf('\n') 'end'])

if numel(totAppData)>2
text( 0.1, 0.5, {'fit function: ax + b', ...
    ['in x direction (blue): ' num2str(resX.p1) 'x + ' num2str(resX.p2) ', R^2 = ' num2str(gofX.rsquare) ] ...
    ['       T_x = ' num2str(Tx) ...
            ' +/- ' num2str(confX(1) / totAppData{1}.consts.Kb*mass * 1e6/2) ' \muK'] ...
    ['       \sigma_{0x} = ' num2str(sqrt(resX.p2)*1000) ' +/- ' num2str(confX(2)*1000) ' mm' ] ...
    ['in y direction (red):  ' num2str(resY.p1) 'x + ' num2str(resY.p2) ', R^2 = ' num2str(gofY.rsquare) ] ...
    ['       T_y = ' num2str(Ty)  ...
            ' +/- ' num2str(confY(1) / totAppData{1}.consts.Kb*mass * 1e6/2) ' \muK'] ...
    ['       \sigma_{0y} = ' num2str(sqrt(resY.p2)*1000) ' +/- ' num2str(confX(2)*1000) ' mm' ] },...
    'Units', 'Normalized'); %     ['average density = ' num2str(n*1e-6) '1/cm^3'] },...
else
    text( 0.1, 0.5, {'fit function: ax + b', ...
    ['in x direction (blue): ' num2str(resX.p1) 'x + ' num2str(resX.p2) ', R^2 = ' num2str(gofX.rsquare) ] ...
    ['       T_x = ' num2str(resX.p1 / totAppData{1}.consts.Kb*mass * 1e6) ' \muK'] ...
    ['       \sigma_{0x} = ' num2str(sqrt(resX.p2)*1000)  ' mm' ] ...
    ['in y direction (red):  ' num2str(resY.p1) 'x + ' num2str(resY.p2) ', R^2 = ' num2str(gofY.rsquare) ] ...
    ['       T_y = ' num2str(resY.p1 / totAppData{1}.consts.Kb*mass * 1e6)  ' \muK'] ...
    ['       \sigma_{0y} = ' num2str(sqrt(resY.p2)*1000) ' mm' ] },...
    'Units', 'Normalized'); %     ['average density = ' num2str(n*1e-6) '1/cm^3'] },...
end
else
     graph =[];
     Result = [];
end

end

function [ x_unique, yMean, yStd, ySem ] = meanXAM( x, y)
% meanXAM calcualtes the mean y value and its standard deviation, per value of x
% AM is for analyzeMeasurement.m function (to avoid duplication with the
% external meanX.m)

[x_unique,m,n]=unique(x);
yMean=zeros(1,length(m));
yStd=zeros(1,length(m));
ySem=zeros(1,length(m));
for j=1:length(m)
    yMean(j)=mean( y(n==j ));
    yStd(j)=std( y(n==j));
    ySem(j)=std( y(n==j))/sqrt(sum(n==j));
%     display(sum(n==j));
end

end