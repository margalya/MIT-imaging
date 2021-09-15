function N = TwoThreeBodyLossTimeEvolution(t, N0, K1, K3T, a, b)
% Solve the time-depenednt euqation for N(t), given K1, and K3 loss parameters. We set K2=0 in this function
% Equation for time evoluation for the number of atoms: N(T)' = dN(t)/dt = -K1*N(t) - K2*<n>*N(t) - K3*<n^2>*N(t)
% Equation can also be written as: N(T)' = dN(t)/dt = -K1*N(t) - K2*N(t)^2/V - K3*N(t)^3/V^2
% Equation taken from "Inelastic Collisions of a Fermi Gas in the BEC-BCS Crossover" - Phys. Rev. Lett. 102, 250402 (2009) - https://www.physics.ncsu.edu/jet/publications/pdf/inelastic.pdf
% v2 - uses n2T2meanPolyFit
% Jacobian = K1 -K3T.*n2meanT - seems unecessary if the interpolation is ok
% a,b are the interpolation parameters for n2meanT2

% JacobiFun = @(t, N, K1, K3T, n2meanT2InversePolyFit) -K1 -K3T.*1./polyval(n2meanT2InversePolyFit, t);
options = odeset('RelTol',1e-9); %, 'Jacobian', @(t,N)JacobiFun(t, N, K1, K3T, n2meanT2InversePolyFit)); % 'Stats','on','OutputFcn',@odeplot);
[~, N] = ode45(@(t,N) ODEFunction(t, N, K1, K3T, a, b), t, N0, options);

end

function dNdt = ODEFunction(t, N, K1, K3T, a, b)
% if isstruct(n2T2meanPolyFit) %avoid problems with matlab's fittype function checker
%     n2meanT2 =  1./polyval(n2meanT2InversePolyFit, t); % <n^2>*T^2 at time t
% else
%     n2mean = n2T2meanPolyFit;
% end
n2meanT2 = a/(2*t+b)*1e34;
dNdt = (-K1*N - K3T*n2meanT2*N)';
% dNdt = (-K1*N - K2*N - K3T*n2meanT2*N)';

% dNdt = (-K1*N - K2*nmean*N - K3*n2mean*N)';
% dNdt = (-K1*N - K2*N^2./K3Volume - K3*N^3./K3Volume)';
% dNdt = (-K1*N - K2*N^2./V - K3*N^3./V.^(2/3))';
end