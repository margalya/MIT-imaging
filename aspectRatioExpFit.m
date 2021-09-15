function aspectRatioExpFit(totAppData, readDir)

for j= 1 : length(totAppData)
    fitType = totAppData{j}.data.fitType;
    aspectRatio(j) = totAppData{j}.data.fits{ fitType }.xUnitSize ./ totAppData{j}.data.fits{ fitType }.yUnitSize; %#ok<AGROW>
%     tempDiff(j) = (totAppData{j}.data.fits{ fitType }.xUnitSize.^2 - totAppData{j}.data.fits{ fitType }.yUnitSize.^2 )./(totAppData{j}.data.fits{ fitType }.xUnitSize.^2 + totAppData{j}.data.fits{ fitType }.yUnitSize.^2 ); %#ok<AGROW>
%     tempDiff2(j) = (totAppData{j}.data.fits{ fitType }.xUnitSize.^2 - totAppData{j}.data.fits{ fitType }.yUnitSize.^2 )./(2 .* totAppData{j}.data.fits{ fitType }.xUnitSize.^2 + totAppData{j}.data.fits{ fitType }.yUnitSize.^2 ); %#ok<AGROW>
    val(j) = totAppData{j}.save.saveParamVal; %#ok<AGROW>
end


% plot
[ val_unique, aspectRatioMean, ~, aspectRatioSem ] = meanXAM( val, aspectRatio);

% aspectRatioMean = aspectRatioMean(val_unique<1000);
% aspectRatioSem = aspectRatioSem(val_unique<1000);
% val_unique = val_unique(val_unique<1000);

aspectRatioMean = aspectRatioMean(val_unique>0);
aspectRatioSem = aspectRatioSem(val_unique>0);
val_unique = val_unique(val_unique>0);

figure('Filename', [readDir '_aspectRatioExpFit.fig']);
errorbar( val_unique, aspectRatioMean, aspectRatioSem, 'ob');

% fit to decaying exponent
startPoint = [range(aspectRatio) mean(val) aspectRatio(end)];
lower = [-2*range(aspectRatio) 0 -1 ];
upper = [2*range(aspectRatio) 2*(max(val)-min(val)) 2*max(aspectRatio)];
s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint, 'Lower', lower, 'Upper', upper, 'Robust', 'on');
f = fittype('a*exp(-t/tau)+c', 'coefficients', {'a', 'tau', 'c'}, 'independent', 't', 'dependent', 'y', 'options', s);
[res, gof] = fit(val_unique', aspectRatioMean', f);

% % % [res, gof, ~] = fit(val_unique', aspectRatioMean', 'exp2');
conf = confint(res);
conf = (conf(2,:)-conf(1,:))/2;

x = linspace(0,max(val_unique),500);
hold on
plot(x, res(x), 'c');

title('Size Aspect Ratio');
xlabel('Time');
ylabel('Aspect Ratio \sigma_X / \sigma_Y +/- SEM');
set(gcf, 'Name', 'Size Aspect Ratio');

text( 0.5, 0.5, {'fit function: a*exp(-t/tau)+c', ...
    [num2str(res.a) 'e^{ -t/' num2str(res.tau) '} + ' num2str(res.c) ', R^2 = ' num2str(gof.rsquare) ] ...
    [] ...
    ['a = ' num2str(res.a) ' +/- ' num2str(conf(1)) ], ...
    ['\tau = ' num2str(res.tau) ' +/- ' num2str(conf(2)) '' ],...
    ['c = ' num2str(res.c) ' +/- ' num2str(conf(3)) ]},...
    'Units', 'Normalized');

save([readDir '_aspectRatio.mat'], 'res', 'conf')

% also analyze the 1D fit results, if fittype is 2DGaussian
% if totAppData{1}.data.fitType == totAppData{1}.consts.fitTypes.twoDGaussian
%     clear aspectRatio val
%     fitType = totAppData{j}.consts.fitTypes.oneDGaussian;
%     for j= 1 : length(totAppData)
%         aspectRatio(j) = totAppData{j}.data.fits{ fitType }.xUnitSize ./ totAppData{j}.data.fits{ fitType }.yUnitSize;
%         val(j) = totAppData{j}.save.saveParamVal;
%     end
%     
%     % plot
%     [ val_unique, aspectRatioMean, ~, aspectRatioSem ] = meanXAM( val, aspectRatio);
%     errorbar( val_unique, aspectRatioMean, aspectRatioSem, 'oblack');
%     
%     % fit to decaying exponent
%     startPoint = [range(aspectRatio) mean(val) aspectRatio(end)];
%     lower = [-2*range(aspectRatio) 0 -1 ];
%     upper = [2*range(aspectRatio) 2*(max(val)-min(val)) 2*max(aspectRatio)];
%     s = fitoptions('Method','NonlinearLeastSquares', 'Startpoint', startPoint, 'Lower', lower, 'Upper', upper);
%     f = fittype('a*exp(-t/tau)+c', 'coefficients', {'a', 'tau', 'c'}, 'independent', 't', 'dependent', 'y', 'options', s);
%     [res1D, gof1D] = fit(val_unique', aspectRatioMean', f);
%     
%     % % % [res, gof, ~] = fit(val_unique', aspectRatioMean', 'exp2');
%     conf = confint(res1D);
%     conf = (conf(2,:)-conf(1,:))/2;
%     
%     x = linspace(0,max(val_unique),50);
%     plot(x, res1D(x), 'r');
%     
%     text( 0.5, 0.1, {'1DG data: fit function: a*exp(-t/tau)+c', ...
%         [num2str(res1D.a) 'e^{ -t/' num2str(res1D.tau) '} + ' num2str(res1D.c) ', R^2 = ' num2str(gof1D.rsquare) ] ...
%         [] ...
%         ['a = ' num2str(res1D.a) ' +/- ' num2str(conf(1)) ], ...
%         ['\tau = ' num2str(res1D.tau) ' +/- ' num2str(conf(2)) ' sec' ],...
%         ['c = ' num2str(res1D.c) ' +/- ' num2str(conf(3)) ]},...
%         'Units', 'Normalized');
%     
%     save([readDir '_aspectRatio.mat'], 'res1D', 'conf')
% end
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