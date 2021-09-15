function conf = calcBounds( res )
% calculates confint (fit errors), and checks cfit object if any parameters are bounded
conf = confint(res);
conf = (conf(2,:)-conf(1,:))/2;
if any(isnan(conf))
    names = coeffnames(res)';
    disp('Bounded fit parameters: ')
    display( names(isnan(conf)) )
end
end

