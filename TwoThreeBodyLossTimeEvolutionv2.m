function N = TwoThreeBodyLossTimeEvolution(t, N0, K1, K2, K3, n2meanInverseSpline) %#ok<DEFNU>
% Solve the time-depenednt euqation for N(t), given K1, K2, and K3 loss parameters.
% Equation for time evoluation for the number of atoms: N(T)' = dN(t)/dt = -K1*N(t) - K2*<n>*N(t) - K3*<n^2>*N(t)
% Equation can also be written as: N(T)' = dN(t)/dt = -K1*N(t) - K2*N(t)^2/V - K3*N(t)^3/V^2
% Equation taken from "Inelastic Collisions of a Fermi Gas in the BEC-BCS Crossover" - Phys. Rev. Lett. 102, 250402 (2009) - https://www.physics.ncsu.edu/jet/publications/pdf/inelastic.pdf
% v2
options = odeset('RelTol',1e-9); %,'Stats','on','OutputFcn',@odeplot);
[~, N] = ode45(@(t,N) ODEFunction(t, N, K1, K2, K3, n2meanInverseSpline), t, N0, options);

end

function dNdt = ODEFunction(t, N, K1, K2, K3, n2meanInverseSpline)
if isstruct(n2meanInverseSpline) %avoid problems with matlab's fittype function checker
    n2mean =  1./fnval(n2meanInverseSpline, t); % <n^2> at time t
else
    n2mean = n2meanInverseSpline;
end
dNdt = (-K1*N - K2*N - K3*n2mean*N)';
% dNdt = (-K1*N - K2*nmean*N - K3*n2mean*N)';
% dNdt = (-K1*N - K2*N^2./K3Volume - K3*N^3./K3Volume)';
% dNdt = (-K1*N - K2*N^2./V - K3*N^3./V.^(2/3))';
end