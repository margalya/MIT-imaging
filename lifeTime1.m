
function lifeTime1(totAppData)

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
[res, gof, ~] = fit(DT_unique', NMean', 'exp1'); %#ok<NASGU>
conf = confint(res);
conf = (conf(2,:)-conf(1,:))/2;

figure('Filename', [totAppData{1}.ui.etReadDir.String '_lifetime.fig']);
errorbar( DT_unique, NMean, NSem, 'ob');
hold on
% plot(DT, res.a*exp(res.b*DT), 'c');
plot(res, 'c');
hold off

title('Life Time (1 exp) calculations');
set(gca,'Ylabel',text('String', 'No. of atoms +/- SEM'));
set(gca,'Xlabel',text('String', 'Dark Time [sec]'));
set(gcf, 'Name', 'Life TIme (1 exp)');

text( 0.5, 0.5, {'fit function: ae^{bt}', ...
    [num2str(res.a) 'e^{ ' num2str(res.b) 't}, R^2 = ' num2str(gof.rsquare) ] ...
    [] ...
    ['N_0 = ' num2str(res.a*1e-6) ' +/- ' num2str(conf(1)*1e-6) '*10^6' ], ...
    ['\tau = ' num2str(-1/res.b) ' +/- ' num2str(conf(2)/res.b^2) ' sec' ]},...
    'Units', 'Normalized');

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