
function lifeTime2(totAppData)

len = length(totAppData);
DT = zeros(1, len); %[sec]
N = zeros(1, len); %[meter]

fitType = totAppData{1}.data.fitType;

for ( i = 1 : len )
%     if ( totAppData{i}.save.saveParam ~= totAppData{i}.consts.saveParams.darkTime )
%         warndlg({['The save parameter in data-' num2str( totAppData{i}.save.picNo) ' is not Dark Time (' num2str(totAppData{i}.consts.saveParams.darkTime) ').']; ...
%             ['It is:' num2str(totAppData{i}.save.saveParam) '.']} , 'Warning', 'modal');
%     end
    DT(i) = totAppData{i}.save.saveParamVal;
    N(i) = totAppData{i}.data.fits{ fitType }.atomsNo;
end
[ DT_unique, NMean, ~, NSem ] = meanXAM( DT, N);
[res gof] = fit(DT_unique', NMean', 'exp2');
conf = confint(res);
conf = (conf(2,:)-conf(1,:))/2;

figure('Filename', [totAppData{1}.ui.etReadDir.String '_lifetime2.fig']);
errorbar( DT_unique, NMean, NSem, 'ob');
hold on
% plot(DT, res.a*exp(res.b*DT) + res.c*exp(res.d*DT), 'c');
plot(res, 'c');
hold off
xlim([min(DT_unique) max(DT_unique)])
ylim([min(NMean) max(NMean)])

title('Life Time (2 exp) calculations');
set(gca,'Ylabel',text('String', 'No. of atoms +/- sem'));
set(gca,'Xlabel',text('String', 'Dark Time [sec]'));
set(gcf, 'Name', 'Life TIme (2 exp)');

text( DT(1)+0.1*DT(end), max(N(:)) *0.5, {'fit function: ae^{bt} + ce^{dt}', ...
    [num2str(res.a) 'e^{ ' num2str(res.b) 't} + ' num2str(res.c) 'e^{ ' num2str(res.d) 't}, R^2 = ' num2str(gof.rsquare) ] ...
    []...
    ['N_0 = ' num2str((res.a+res.c)*1e-6) ' +/- ' num2str( sqrt(conf(1)^2+conf(3)^2)*1e-6) '*10^6' ], ...
    ['A_1 = ' num2str(res.a*1e-6) ' +/- ' num2str(conf(1)*1e-6) '*10^6'], ...
    ['\tau_1 = ' num2str(-1/res.b) ' +/- ' num2str(conf(2)/res.b^2) ' sec' ], ...
    ['A_2 = ' num2str(res.c*1e-6) ' +/- ' num2str(conf(3)*1e-6) '*10^6'], ...
    ['\tau_2 = ' num2str(-1/res.d) ' +/- ' num2str(conf(4)/res.d^2) ' sec' ]});

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