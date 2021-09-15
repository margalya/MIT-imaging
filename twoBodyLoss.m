
function twoBodyLoss(totAppData)
% fit the loss rate to a two-body loss function, assuming size is constant
len = length(totAppData);
DT = zeros(1, len); %[sec]
N = zeros(1, len); %[meter]

fitType = totAppData{1}.data.fitType;

for i = 1 : len 
    DT(i) = totAppData{i}.save.saveParamVal; %dark time
    N(i) = totAppData{i}.data.fits{ fitType }.atomsNo; %number of atoms
end
[ DT_unique, NMean, ~, NSem ] = meanXAM( DT, N);

% fit
K2VStartPoint = abs(mean(diff(NMean)./diff(DT_unique)./NMean(1:end-1).^2))*1e6;
startPoint = [max(NMean) K2VStartPoint];
lower = [0 0 ];
upper = [2*max(NMean) 10*K2VStartPoint];
s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint, 'Lower', lower, 'Upper', upper);
f = fittype('N0/(1 + K2V*1e-6*N0*t)', 'coefficients', {'N0', 'K2V'}, 'independent', 't', 'dependent', 'y', 'options', s); %N(t) = N0/(1 + k2*N0*t/V)
[res, gof] = fit(DT_unique', NMean', f);

conf = confint(res);
conf = (conf(2,:)-conf(1,:))/2;

figure('Filename', [totAppData{1}.ui.etReadDir.String '_twoBodyLoss.fig']);
errorbar( DT_unique, NMean, NSem, 'ob');
hold on
plot(res, 'c');
% tGrid = [1e-4 5e-4 1e-3 5e-3 1e-2 2.5e-2 5e-2 7.5e-2 1e-1 linspace(0.1, max(DT_unique))];
% plot(tGrid, res(tGrid), 'c');
xlim([min(DT_unique)-0.05*range(DT_unique) max(DT_unique)+0.05*range(DT_unique)])
ylim([min(NMean)-0.05*range(NMean) max(NMean)+0.05*range(NMean)])

title('Two-body loss fit');
set(gca,'Ylabel',text('String', 'No. of atoms +/- sem'));
set(gca,'Xlabel',text('String', 'Dark Time [sec]'));
set(gcf, 'Name', 'Two-body loss');

text( 0.3, 0.8, {'fit function: N_0/(1 + K_2/V N_0 t)', ...
    ['R^2 = ' num2str(gof.rsquare) ] ...
    []...
    ['N_0 = (' num2str((res.N0)*1e-6) ' +/- ' num2str( conf(1)*1e-6) ')*10^6 atoms' ], ...
    ['K_2/V = (' num2str(res.K2V) ' +/- ' num2str(conf(2)) ')*10^{-6}/sec/atom']},...
    'Units', 'Normalized');

save([totAppData{1}.ui.etReadDir.String '_twoBodyLoss.mat'], 'res', 'conf')

%%%%%%%%%%% fit - one-body + two-body loss, analytic %%%%%%%%%%
% looks like there's a problem - fit dosen't converge, should check
% fit
% K1 = 1/31;
% startPoint2 = [max(NMean) res.K2V];
% lower2 = [0 res.K2V/10 ];
% upper2 = [2*max(NMean) 10*res.K2V];
% s2 = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint2, 'Lower', lower2, 'Upper', upper2);
% f2 = fittype('K1/(exp(K1*t)*K1/N0+(exp(K1*t) - 1)*K2V*1e-6)', 'coefficients', {'N0', 'K2V'}, 'independent', 't', 'dependent', 'y', 'options', s2, 'problem', {'K1'}); %N(t) = N0/(1 + k2*N0*t/V)
% [res2, gof2] = fit(DT_unique', NMean', f2, 'problem', {K1});
% 
% conf2 = confint(res2);
% conf2 = (conf2(2,:)-conf2(1,:))/2;
% 
% % figure('Filename', [totAppData{1}.ui.etReadDir.String '_twoBodyLoss.fig']);
% % errorbar( DT_unique, NMean, NSem, 'ob');
% % hold on
% plot(res2, 'r');
% xlim([min(DT_unique)-0.05*range(DT_unique) max(DT_unique)+0.05*range(DT_unique)])
% ylim([min(NMean)-0.05*range(NMean) max(NMean)+0.05*range(NMean)])
% 
% title('Two-body loss fit');
% set(gca,'Ylabel',text('String', 'No. of atoms +/- sem'));
% set(gca,'Xlabel',text('String', 'Dark Time [sec]'));
% set(gcf, 'Name', 'Two-body loss');
% 
% text( 0.3, 0.4, {'fit function2 2: K1/(exp(K1*t)*K1/N0+(exp(K1*t) - 1)*K2/V)', ...
%     ['R^2 = ' num2str(gof2.rsquare) ] ...
%     []...
%     ['N_0 = (' num2str((res2.N0)*1e-6) ' +/- ' num2str( conf2(1)*1e-6) ')*10^6 atoms' ], ...
%     ['\tau = 1/K_1 = ' num2str((1/res2.K1)) ' sec (fixed value)' ], ...
%     ['K_2/V = (' num2str(res2.K2V) ' +/- ' num2str(conf2(2)) ')*10^{-6}/sec/atom']},...
%     'Units', 'Normalized');

%%%%%%%%%%% fit - one-body + two-body loss %%%%%%%%%%
% K3 = 0; %assuming only three body loss
% V = 2.9049e-16; %Volume: V = sigma_r*sigma_r*sigma_z, taken for 3V X ODT, and 300uK temperature.
% 
% % startPoint2 = [res.N0 1 8];
% % lower2 = [0 0.001 0.1];
% % upper2 = [4*res.N0 100 20];
% 
% startPoint2 = [res.N0 1/60 96];
% lower2 = [0 0 1];
% upper2 = [4*res.N0 10 300];
% [holdTime, SizeX ] = extractFigFileData([totAppData{1}.ui.etReadDir.String '_sizeX.fig']);
% [~, SizeY ] = extractFigFileData([totAppData{1}.ui.etReadDir.String '_sizeY.fig']);
% s2 = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint2, 'Lower', lower2, 'Upper', upper2);
% f2 = fittype('TwoThreeBodyLossTimeEvolution(t, N0, K1, K2, K3, V, sigma_r, sigma_z, tgrid)',...
%     'coefficients', {'N0', 'K1', 'K2'}, 'independent', 't', 'dependent', 'y', 'options', s2, 'problem', {'K3', 'V', 'sigma_r', 'sigma_z', 'tgrid'});  % , 'N0'
% [res2, gof2] = fit(DT_unique', NMean', f2, 'problem', {K3, V, SizeX, SizeY, holdTime}); % , res.N0
% 
% conf2 = confint(res2);
% conf2 = (conf2(2,:)-conf2(1,:))/2;
% % figure('Filename', [totAppData{1}.ui.etReadDir.String '_one+threeBodyLoss.fig']);
% % errorbar( DT_unique, NMean, NSem, 'ob');
% plot(DT_unique, res2(DT_unique), 'r');
% 
% title('One-Body + two-body loss fit');
% % set(gca,'Ylabel',text('String', 'No. of atoms +/- sem'));
% % set(gca,'Xlabel',text('String', 'Dark Time [sec]'));
% % set(gcf, 'Name', 'Three-body loss');
% 
% %note: we scale K2 and K3 inside 'TwoThreeBodyLossTimeEvolution' function the following way:
% % K2 = K2fit*V*1e-6;
% % K3 = K3fit*V^2*1e-12;
% text( 0.3, 0.4, {'fit function 2: Numerical 1+2 body loss', ...
%     ['R^2 = ' num2str(gof2.rsquare) ] ...
%     []...
%     ['N_0 = ' num2str((res2.N0)*1e-6) ' +/- ' num2str( conf2(1)*1e-6) '*10^6 atoms' ], ...
%     ['\tau = 1/K_1 = ' num2str((1/res2.K1)) ' +/- ' num2str( conf2(2)) ' sec' ], ...
%     ['K_2/V = ' num2str(res2.K2) ' +/- ' num2str(conf2(3)) ' 10^{-6}/sec/atom^2']},...
%     'Units', 'Normalized');

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