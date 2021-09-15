
function threeBodyLoss(totAppData)
% fit the loss rate to a three-body loss function, assuming size is constant
% v2: in this version, we use numerical <n^2> in the nbumerical fit, but it depends strongly on the interpolation
len = length(totAppData);
DT = zeros(1, len); %[sec]
N = zeros(1, len); %[meter]

fitType = totAppData{1}.data.fitType;

for i = 1 : len 
    DT(i) = totAppData{i}.save.saveParamVal; %dark time
    N(i) = totAppData{i}.data.fits{ fitType }.atomsNo; %number of atoms
end
[ DT_unique, NMean, ~, NSem ] = meanXAM( DT, N);

% user input
prompt = {'States 6 to 1 scaling factor:', 'Original time-scale unit (msec / sec) :', 'ODT control voltage [V]:', 'TOF [ms]:'};
dlgtitle = 'User input';
dims = [1 35];
definput = {'2.6036', 'msec', '2', '0.5'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
NMean = str2double(answer{1})*NMean; %calibrate state 1 count using state 6
if strcmp(answer{2}, 'msec')
    DT_unique = DT_unique*1e-3; %change to seconds (usually the relevant scale)
end

% fit - pure three-body loss
startPoint = [max(NMean) 1];
lower = [0 1e-10 ];
upper = [2*max(NMean) 200];
s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint, 'Lower', lower, 'Upper', upper);
f = fittype('(2*K3V2*1e-12*t + 1/N0^2)^(-1/2)', 'coefficients', {'N0', 'K3V2'}, 'independent', 't', 'dependent', 'y', 'options', s); %N(t) = (2*K3/V^2*t + 1/N0^2)^(-1/2)
[res, gof] = fit(DT_unique', NMean', f);

conf = confint(res);
conf = (conf(2,:)-conf(1,:))/2;

figure('Filename', [totAppData{1}.ui.etReadDir.String '_threeBodyLoss.fig']);
errorbar( DT_unique, NMean, NSem, 'ob');
hold on
plot(res, 'c');
xlim([min(DT_unique)-0.05*range(DT_unique) max(DT_unique)+0.05*range(DT_unique)])
ylim([min(NMean)-0.05*range(NMean) max(NMean)+0.05*range(NMean)])

title('Three-body loss fit');
set(gca,'Ylabel',text('String', 'No. of atoms +/- sem'));
set(gca,'Xlabel',text('String', 'Dark Time [sec]'));
set(gcf, 'Name', 'Three-body loss');

hText = text( 0.3, 0.8, {'fit function: (2*K3/V^2*t + 1/N0^2)^{-1/2}', ...
    ['R^2 = ' num2str(gof.rsquare) ] ...
    []...
    ['N_0 = (' num2str((res.N0)*1e-6) ' +/- ' num2str( conf(1)*1e-6) ')*10^6 atoms' ], ...
    ['K_3/V^2 = (' num2str(res.K3V2) ' +/- ' num2str(conf(2)) ')*10^{-12}/sec/atom^2'],...
    ['Loss rate at t=0 = K_3/V^2*N_0^2 = (' num2str(res.K3V2*1e-12*(res.N0)^2) ')/sec']},...
    'Units', 'Normalized');

save([totAppData{1}.ui.etReadDir.String '_threeBodyLoss.mat'], 'res', 'conf')

%%%%%%%%%%% fit - one-body + three-body loss, time dependent volume %%%%%%%%%%
if ~isempty(answer{3}) && ~isempty(answer{4})
% K1 = 0; %
K1 = 1/31; % data from 2020-02-06
K2 = 0; %assuming no two-body loss
% K3StartPoint = abs(mean(diff(NMean)./diff(DT_unique)./NMean(1:end-1).^3)*V0^2)*1e30;

% startPoint2 = [res.N0 1 8];
% lower2 = [0 0.001 0.1];
% upper2 = [4*res.N0 100 20];

if exist([totAppData{1}.ui.etReadDir.String '_sizeX.fig'],'File') && exist([totAppData{1}.ui.etReadDir.String '_sizeY.fig'],'File')
    % extract size data vs hold time
    [~, SigmaX ] = extractFigFileData([totAppData{1}.ui.etReadDir.String '_sizeX.fig']);
    [~, SigmaY ] = extractFigFileData([totAppData{1}.ui.etReadDir.String '_sizeY.fig']);
else
   warndlg(['Three-body loss, time-dependent volume:' newline 'Cannot find SizeX and SizeY fig files for current folder' newline 'Aborting second fit'])
   return;
end

% Vt = V0.*(sigma_r.^2 .* sigma_z)./(sigma_r(1).^2 .* sigma_z(1)); %calculate the time-dependent volume according to V(t)=V0*(sigma_r(t).^2 .* sigma_z(t))./(sigma_r(t=0).^2 .* sigma_z(t=0))
output = calcXODTDensity(NMean, str2double(answer{3}), str2double(answer{4}), SigmaX, SigmaY);
% n2meanThNumeric = output.n2meanThNumeric;
% n2meanPolyFit = polyfit(holdTime, n2meanThNumeric, 7); % use second-order polynomial for time-dependence of the K3 Volume
% figure; plot(holdTime, n2meanThNumeric, 'o'); hold on; plot( linspace(min(holdTime), max(holdTime)), polyval(n2meanPolyFit, linspace(min(holdTime), max(holdTime))), 'r')
n2meanInverseSpline = csaps(DT_unique, 1./output.n2meanThNumeric); %spline for the inverse of the numeric <n^2>, which is usually smoother than the regular numeric <n^2>
%%% plot n2meanInverseSpline result
tempFig = figure; plot(DT_unique, 1./output.n2meanThNumeric, 'o'); hold on; plot(linspace(min(DT_unique),max(DT_unique)), fnval(n2meanInverseSpline,linspace(min(DT_unique),max(DT_unique))))
xlabel('Hold time [sec]')
ylabel('1/<n^2>')
pause(3);
close(tempFig);

startPoint2 = [res.N0 res.K3V2*1e-12*output.K3Volume(1)*1e42]; %1/30
lower2 = [0 0];
upper2 = [4*res.N0 100];
s2 = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint2, 'Lower', lower2, 'Upper', upper2);
f2 = fittype('TwoThreeBodyLossTimeEvolution(t, N0, K1, K2, K3*1e-42, n2meanInverseSpline)',...
    'coefficients', {'N0', 'K3'}, 'independent', 't', 'dependent', 'y', 'options', s2, 'problem', {'K1', 'K2', 'n2meanInverseSpline'});
[res2, gof2] = fit(DT_unique', NMean', f2, 'problem', {K1, K2, n2meanInverseSpline});

conf2 = confint(res2);
conf2 = (conf2(2,:)-conf2(1,:))/2;
% figure('Filename', [totAppData{1}.ui.etReadDir.String '_one+threeBodyLoss.fig']);
% errorbar( DT_unique, NMean, NSem, 'ob');
plot(DT_unique, res2(DT_unique), 'r');

title('One-Body + three-body loss fit');
set(gca,'Ylabel',text('String', 'No. of atoms +/- sem'));
set(gca,'Xlabel',text('String', 'Hold time [sec]'));
set(gcf, 'Name', 'Three-body loss');

%%% calcualte K3 values in cm^3
K3 = res.K3V2 .* output.K3Volume(1);
K3Err = conf(2) .* output.K3Volume(1);
K3Numeric = res2.K3*1e-42*1e12;
K3NumericErr = conf2(2)*1e-42*1e12;

%%% update previous text regarding the analytic result of K3, using the
delete(hText);
text( 0.3, 0.8, {'fit function: (2*K3/V^2*t + 1/N0^2)^{-1/2}', ...
    ['R^2 = ' num2str(gof.rsquare) ] ...
    []...
    ['N_0 = (' num2str((res.N0)*1e-6) ' +/- ' num2str( conf(1)*1e-6) ')*10^6 atoms' ], ...
    ['K_3 = (' num2str(K3) ' +/- ' num2str(K3Err) ')*cm^6/sec'],...
    ['Loss rate at t=0 = K_3/V^2*N_0^2 = (' num2str(res.K3V2*1e-12*(res.N0)^2) ')/sec']},...
    'Units', 'Normalized');

text( 0.3, 0.4, {'fit function 2: Numerical 1+3 body loss', ...
    ['R^2 = ' num2str(gof2.rsquare) ] ...
    []...
    ['N_0 = (' num2str((res2.N0)*1e-6) ' +/- ' num2str( conf2(1)*1e-6) ')*10^6 atoms' ], ... %    ['\tau = 1/K_1 = ' num2str((1/res2.K1)) ' +/- ' num2str( conf2(2)) ' sec' ], ...
    ['K_3/V_0^2 = (' num2str(res2.K3*1e-42*1e12/output.K3Volume(1)) ' +/- ' num2str(conf2(2)*1e-42*1e12/output.K3Volume(1)) ')*10^{-12}/sec/atom^2'],...
    ['K_3 = (' num2str(K3Numeric) ' +/- ' num2str(K3NumericErr) ') cm^6/sec']},...
    'Units', 'Normalized');
save([totAppData{1}.ui.etReadDir.String '_threeBodyLoss.mat'], 'res2', 'conf2', 'K3', 'K3Err', 'K3Numeric', 'K3NumericErr', '-append' )

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