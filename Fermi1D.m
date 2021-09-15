function ret = Fermi1D( n0, q, x, R )
% calcualte the 1D Fermi profile using a smoothing spline interpolation of
% the polylog_5/2 function, with a precision of 5e-5

%% approximation
persistent polylog5halfSpline
%spline of the polylog function to speedup computation, in the range -19.63:1
if isempty(polylog5halfSpline)
    polylog5halfSpline = load('polylog5halfSpline.mat', 'polylog5halfSpline');
end
arg1 = -exp(q - x.^2./R^2.*((1+exp(q))./exp(q).*log(1+exp(q))) );
arg2 = -exp(q);
if any([arg1<-2017 ; arg1>1 ; arg2<-2017 ; arg2>1])
    disp('Warning: PolyLog argument out of spline approximation range (-2017 < x < 1)')
    return;
%     ret = n0*polylog(5/2, -exp(q - x.^2./R^2.*((1+exp(q))/exp(q)*log(1+exp(q))) ))./polylog(5/2, -exp(q));
else
    ret = n0*polylog5halfSpline.polylog5halfSpline( arg1 )./polylog5halfSpline.polylog5halfSpline( arg2 );
end
% assignin('base','arg1Temp', arg1);
% assignin('base','arg2Temp', arg2);
% evalin('base', 'arg1 = [arg1; arg1Temp]; arg2 = [arg2; arg2Temp]')
%% exact matlab function

% ret = n0*polylog(5/2, -exp(q - x.^2./R^2.*((1+exp(q))/exp(q)*log(1+exp(q))) ))./polylog(5/2, -exp(q));

%% generate polylog spline data

% clear;
% interpolationRange = 2300;
% x = 1:2000; %number of interpolation points
% dx = 1./x; % shift in interpolation vector
% dx = dx/sum(dx)*interpolationRange; % normalize to range
% x = cumsum(dx);
% x = x-x(end)+1;
% 
% % y = polylog(5/2,x);
% % cftool(x, y) % we generate interpolation using smoothing spline interpolation, saved into file 'polylog5halfSpline'
% 
% %% test interpolation quality
% % the values interpolationRange = 2300; and x = 1:2000; gave a maximum error
% % of 2.5e-7 between interpolation and exact function
% % load(polylog5halfSpline)
% xd = downsample(x,20)-0.003;
% dif = (polylog(5/2,xd)') - polylog5halfSpline(xd);
% figure; plot(xd, dif,'o')

end