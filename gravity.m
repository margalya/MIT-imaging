
function graph = gravity(totAppData)
if numel(totAppData) > 2
len = length(totAppData);
tof = zeros(1, len); %[sec]
y = zeros(1, len); %[meter]

fitType = totAppData{1}.data.fitType;

for ( i = 1 : len )
%     if ( totAppData{i}.save.saveParam ~= totAppData{i}.consts.saveParams.TOF )
%         warndlg({['The save parameter in data-' num2str( totAppData{i}.save.picNo) ' is not TOF (' num2str(totAppData{i}.consts.saveParams.TOF) ').']; ...
%             ['It is:' num2str(totAppData{i}.save.saveParam) '.']} , 'Warning', 'modal');
%     end
    tof(i) = totAppData{i}.save.saveParamVal;
    y(i) = totAppData{i}.data.fits{ fitType }.yCenter;
end

tof = tof/1000; % change to seconds
y = y  * totAppData{i}.data.camera.yPixSz; % change to meters

[res gof out] = fit(tof', y', 'poly2'); %#ok<NASGU>
conf = confint(res);
conf = (conf(2,:)-conf(1,:))/2;

s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [9.8 y(1)], 'Lower', [0 0], 'Upper', [20 y(end)]);
f = fittype('0.5*g*x^2+x0', 'coefficients', {'g', 'x0'}, 'independent', 'x', 'dependent', 'y', 'options', s);
[res1 gof1 out1] = fit(tof', y', f); %#ok<NASGU>
conf1 = confint(res1);
conf1 = (conf1(2,:)-conf1(1,:))/2;
s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', [0 y(1)], 'Lower', [0 0], 'Upper', [0.1 y(end)]);
f = fittype('4.9*x^2+v0*x+x0', 'coefficients', {'v0', 'x0'}, 'independent', 'x', 'dependent', 'y', 'options', s);
[res2 gof2 out2] = fit(tof', y', f); %#ok<NASGU>
conf2 = confint(res2);
conf2 = (conf2(2,:)-conf2(1,:))/2;

graph = figure
plot( tof, y, 'ob');
hold on
plot(res, 'c');
% plot(tof, res.p1*tof.^2+res.p2*tof+res.p3, 'c');
hold off

title('Gravity calculations');
set(gca,'Ylabel',text('String', 'y_0 [m]'));
set(gca,'Xlabel',text('String', 'TOF [sec]'));
set(gcf, 'Name', 'Gravity');

MRb = 1.44e-25; % kg
kB = 1.38065e-23; % J/K
if numel(totAppData) > 3
text( tof(1)+0.1*tof(2), (max(y)-min(y)) *0.9+min(y), {'fit function: at^2 + bt + c', ...
    [num2str(res.p1) 't^2 + ' num2str(res.p2) 't + ' num2str(res.p3) ', R^2 = ' num2str(gof.rsquare) ] ...
    ['g = ' num2str(res.p1*2) ' +/- ' num2str(conf(1)) ' m/sec^2' ], ...
    ['v_0 = ' num2str(res.p2*100) ' +/- ' num2str(100*conf(2)) ' cm/sec' ], ...
    ['v_0 = ' num2str(MRb*res.p2^2/kB*1e6) ' +/- ' num2str(MRb*res.p2/kB*1e6*conf(2)) ' \muK (T=mv^2/k_B)' ], ...
    ['y_0 = '  num2str(res.p3*1e3 - totAppData{1}.data.camera.chipStart * totAppData{1}.data.camera.yPixSz*1e3) ...
            ' +/- ' num2str(conf(3)*1e3) ' mm (from the chip)' ]});
end

graph = [graph figure]
plot( tof, y, 'ob');
hold on
plot(res1, 'c');
plot(res2, 'g');
% plot(tof, 0.5*res1.g*tof.^2+res1.x0, 'c');
% plot(tof, 4.9*tof.^2+res2.v0*tof+res2.x0, 'g');
hold off

title('Gravity calculations - fixed g / no v_0');
set(gca,'Ylabel',text('String', 'y_0 [m]'));
set(gca,'Xlabel',text('String', 'TOF [sec]'));
set(gcf, 'Name', 'Gravity');

if numel(totAppData) > 3
text( tof(1)+0.1*tof(2), (max(y)-min(y)) *0.7+min(y), {'fit function (cyan): 0.5gt^2 + y_0', ...
    ['0.5*' num2str(res1.g) 't^2 + ' num2str(res1.x0) ', R^2 = ' num2str(gof1.rsquare) ] ...
    ['g = ' num2str(res1.g) ' +/- ' num2str(conf1(1)) ' m/sec^2' ], ...
    ['y_0 = '  num2str(res1.x0*1e3 - totAppData{1}.data.camera.chipStart * totAppData{1}.data.camera.yPixSz*1e3) ...
         ' +/- ' num2str(conf1(2)*1e3) ' mm (from the chip)' ] ...
    [] ...
    'fit function (green): 4.9t^2 + v_0t + y_0', ...
    ['4.9t^2 + ' num2str(res2.v0) 't + ' num2str(res2.x0) ', R^2 = ' num2str(gof2.rsquare) ] ...
    ['v_0 = ' num2str(res2.v0*100) ' +/- ' num2str(conf2(1)*100) ' cm/sec' ], ...
    ['v_0 = ' num2str(MRb*res2.v0^2/kB*1e6) ' +/- ' num2str(MRb*res.p2/kB*1e6*conf2(1)) ' \muK (T=mv^2/k_B)' ], ...
    ['y_0 = ' num2str(res2.x0*1e3 - totAppData{1}.data.camera.chipStart * totAppData{1}.data.camera.yPixSz*1e3) ...
         ' +/- ' num2str(conf2(2)*1e3) ' mm (from the chip)' ]});
end
else 
    graph = [];
end