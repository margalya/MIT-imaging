function [out] = BoseFunc(n, z)
% Calculate the Bose function g_n(z)
% This is the same as the polylog(n,z) function!!
fun = @(x) x.^(n-1) ./ ( z^(-1) .* exp(x) -1 );
out = 1./gamma(n) .* integral( fun ,0, Inf, 'AbsTol', 1e-10, 'RelTol', 1e-10);

end

