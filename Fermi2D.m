function ret = Fermi2D( n0, q, x, Rx, y, Ry )
% calcualte the 2D Fermi profile using a smoothing spline interpolation of the polylog_2 function
% does not work well for very large q values (positive or negative), then
% arg1 or arg2 goes out of approximation range.
persistent polylog2Spline
%spline of the polylog function to speedup computation, in the range -2017:1
if isempty(polylog2Spline)
    polylog2Spline = load('polylog2Spline.mat', 'polylog2Spline');
end
arg1 = -exp(q - ( x.^2./Rx^2 + y.^2./Ry^2).*((1+exp(q))./exp(q).*log(1+exp(q))) );
arg2 = -exp(q); %only makes sure that the peak amplitude is correct, does not affect the value of the fugacity
if any([arg1(:)<-2017 ; arg1(:)>1 ; arg2(:)<-2017 ; arg2(:)>1])
    disp('Warning: PolyLog argument out of spline approximation range (-2017 < x < 1), no output returned')
    return
end

ret = n0.*polylog2Spline.polylog2Spline( arg1 )./polylog2Spline.polylog2Spline( arg2 ); %returns a vector for a matrix input, but is seems that it dosen't matter as we just need the sum

% validate results compared to excat function - seems to be a maximum difference around 5e-8 
% ret2 = n0*(abs(polylog(2, arg1 ))./abs(polylog(2, arg2 )))'; %returns a vector for a matrix input, but is seems that it dosen't matter as we just need the sum
% figure;
% plot(ret-ret2,'o')

% element return order for cfit structure:
% A = ([1:3;4:6;7:9;10:12]') - this array has the correct order of element (1:12)
% then using y=x fittedmodel function, we get back A by using reshape(fittedmodel(A),size(A))
% test:
% reshape(fittedmodel(A),size(A))-A
% ans =
%      0     0     0     0
%      0     0     0     0
%      0     0     0     0

% %% generate polylog 2 interpolation 
% clear;
% interpolationRange = 2300;
% x = 1:2000; %number of interpolation points
% dx = 1./x; % shift in interpolation vector
% dx = dx/sum(dx)*interpolationRange; % normalize to range
% x = cumsum(dx);
% x = x-x(end)+1;
% 
% y = polylog(2,x);
% cftool(x, y) % we generate interpolation using smoothing spline interpolation, saved into file 'polylog5halfSpline'
% 
% %% test polylog 2 interpolation quality
% % the values interpolationRange = 2300; and x = 1:2000; gave a maximum error
% % of 2e-7 between interpolation and exact function (at x=-2017). towards x=0 the error is around 1e-9
% if ~exist(polylog2Spline)
%     load(polylog2Spline)
% end
% 
% xd = downsample(x,20)+0.003;
% dif = abs(((polylog(2,xd)') - polylog2Spline.polylog2Spline(xd))./polylog(2,xd)');
% figure; plot(xd, dif,'o')
% set(gcf,'yscale', 'log')

end