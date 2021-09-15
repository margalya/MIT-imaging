function [output] = calcXODTDensity(N, Vcontrol, TOF, SigmaX, SigmaY, varargin)
% calculate K3 volume for X ODT.
% N is the number of atoms
% Vcontrol is the X ODT control voltage
% TOF is time-of-flight is msec
% SigmaX and SigmaY are the TOF cloud sizes, in mm
% output is a structure, with fields containing relevant parameters
% K3 coefficient volume = 2*24*sqrt(3)*pi^3*sigma^6 [m^6]
% if varargin = 'nuemric', then the function computes the numerical density

kB = 1.3806504e-23; %J/K
hbar = 1.054571628e-34;
a0 = 0.52917720859e-10; %Bohr radius [m]
m = 6/6.022e26;
g = 2*pi*5.8724e6;
lambdaAtom = 670.977338e-9;
c = 2.99792458e8; % speed of light [m/s]
mu0 = 1.25663706212e-6; %Vacuum permeability [H/m]
muB = 9.27400915e-24; %Bohr magneton, J/T

Trans = 0.9936 * 0.999* 0.9926 * 0.995 * 0.9975; %Transmission of two lenses (AC508-500-C-ML and AC508-250-C-ML), one laser line mirror, and two surfaces of the vacuum window (from Aviv's thesis, Appendix C, p. 125).
P = Trans.*(1.7893.*Vcontrol + 0.0014); % taking into account optical transmissions [W]
w0 = 7.9786e-6;

lambdaIR = 1064e-9;
omega0 = 2*pi*c/lambdaAtom;
omega = 2*pi*c/lambdaIR;

zR = pi*w0*w0/lambdaIR * 2684/2019; %the extra factor is needed in order to explain the low measured axial frequency (2020-02-27)
U0 = 3*pi*c^2/(2*omega0^3)*( g/(omega0-omega) + g/(omega0+omega))*2.*P./(pi*w0*w0);

omega_r = sqrt(4.*U0./(m*w0*w0));
omega_z = sqrt(2.*U0./(m*zR^2));

Tr = (SigmaX*1e-3).^2./(1 + omega_r.^2.*(TOF.*1e-3).^2)* m ./ kB .* omega_r.^2; % radial temperature
Tz = (SigmaY*1e-3).^2./(1 + omega_z.^2.*(TOF.*1e-3).^2)* m ./ kB .* omega_z.^2; % axial temperature
T = mean([Tr; Tr; Tz], 1); %temperature as a function of hold time

% far-field temperaute, without the correction of the finite TOF
TFarField = mean([(SigmaX*1e-3./(TOF.*1e-3)).^2; (SigmaX*1e-3./(TOF.*1e-3)).^2; (SigmaY*1e-3./(TOF.*1e-3)).^2], 1) .* m ./ kB; %temperature as a function of hold time

% disp(['T = ' num2str(T*1e6) ' uk'])
% Thermal cloud numbers
sigma_r = sqrt(kB*Tr./(m*omega_r.^2));
sigma_z = sqrt(kB*Tz./(m*omega_z.^2));
nmeanTh = N./(8*pi^(3/2) * sigma_r .* sigma_r .* sigma_z); % density-weighted density, thermal atoms limit (Maxwell-Boltzmann) [1/m^3].
npeakTh = N./((2*pi)^(3/2) * sigma_r .* sigma_r .* sigma_z); % peak harmonic density, thermal atoms limit (Maxwell-Boltzmann) [1/m^3].
n2meanTh = N.^2./(24*sqrt(3)*pi^3 .* sigma_r.^2 .* sigma_r.^2 .* sigma_z.^2); % <n^2> [1/m^6]
% output.K3Volume = ((sigma_r).^2.*sigma_z).^2*24*sqrt(3)*pi^3; % K3 coefficient volume = 2*24*sqrt(3)*pi^3*sigma^6 [m^6]
K2Volume = 8 * pi^(3/2) * (sigma_r).^2.*sigma_z; % K3 coefficient volume = 2*24*sqrt(3)*pi^3*sigma^6 [m^6]
K3Volume = 24 * sqrt(3) * pi^3 * ((sigma_r).^2.*sigma_z).^2; % K3 coefficient volume = 2*24*sqrt(3)*pi^3*sigma^6 [m^6]
OD_r = 3*lambdaAtom.^2*N./( (2*pi)^2*sigma_r.*sigma_z ); %Resonant sigma+ light optical density, along the radial direction
OD_z = 3*lambdaAtom.^2*N./( (2*pi)^2*sigma_r.*sigma_r ); %Resonant sigma+ light optical density, along the axial direction


% calculate dipolar loss rate for lithium-6, state 6, using chromium paper equations
% We use sigma = 8*pi*a^2, since the 6-6 collision consits of indistinguishable atoms; 1+0.5 for fermions; and h(x=1)=1/2 (k_i=k_f) (i.e. almost not Zeeman energy)
F = 3/2;
J = 1/2; %ground state
I = 1;
gJ = 2.0023010;
gF =  gJ*( F*(F + 1) - I*(I + 1) + J*(J + 1))./(2*F*(F+1));
sigmaD_1 = 8*pi./15*F^3*(mu0 *(gF *muB)^2 * m / (4*pi*hbar^2)).^2*(1 + 0.5);
sigmaD_2 = 8*pi./15*F^2*(mu0 *(gF *muB)^2 * m / (4*pi*hbar^2)).^2*(1 + 0.5);
vrelTh = sqrt(16*kB*T./(pi*m)); % Mean relative velocity, thermal atoms
L2dipolar = (sigmaD_1 +2*sigmaD_2)* vrelTh; % Dipolar loss rate [m^3/sec]. One state 6 atom is lost in process sigmaD_1, and two state 6 atoms are lost in process sigmaD_2

% Degenerate FD numbers at T=0
Ef = hbar .* (omega_z .* omega_r .* omega_r .* 6 .* N).^(1/3); %Fermi energy [J]
Tf = Ef./kB; % Fermi temeprature [K]
% Rfr = sqrt(Ef/(m*omega_r^2)); %Fermi radius along r [m]
% Rfz = sqrt(Ef/(m*omega_z^2)); %Fermi radius along z [m]
% n0FD = 8*N/(pi^2 * Rfr * Rfr * Rfz); %peak density at T=0, FD statistics [1/m^3]
% nmeanFD = 4096 * N/(315 * pi^2 * Rfr * Rfr * Rfz); % density-weighted density, T=0 limit [1/m^3].

% calcualte p-wave collision rate
a_p = -35*a0; %p-wave scattering length for Li 6 (Mukaiyama PRA 2013)
V_p = a_p^3; %p-wave scattering volume
% gamma_pTh = (288*kB*m^3*N*T*V_p^2*omega_r*omega_r*omega_z)/(pi*4.1*hbar^4);


tau_pTh = (pi*4.1*hbar^4)./(288*kB*m^3*N.*T*V_p^2 .* omega_r .* omega_r .* omega_z); %p-wave thermalization time [sec]
% gamma_pFD = (169869312*sqrt(2)*kB*m^3*N*(Ef/kB)*V_p^2*omega_r*omega_r*omega_z)/(588245*pi^2*4.1*hbar^4);

% p-wave rate calcualtion using the simple sigma*n*v*2/alpha equation
vrel = sqrt(2*8*kB*T/(pi*m)); %mean relative velocity between colliding atoms
k = m*vrel/hbar;
sigma_p = 12*pi*k.^4*V_p^2; %(Ouerdane, Jamieson 2009)
snv = 2/4.1*sigma_p .* nmeanTh .* vrel;

% anharmonic corrections using integral results from Mathematica - omitted, it's done later inside matlab
% load anharmonicDensityCorrection, the input variable is kB*T/U0
% this was calculated using integration limits on r and z up to 95% of potential depths
% load('E:\Dropbox (MIT)\Pauli Blocking\p-wave data analysis\anaharmonicDensityCorrection.mat', 'anharmonicDensityCorrection');
% nmeanTh  = nmeanTh * anharmonicDensityCorrection(kB*T/U0); %anharmonic correction to the density
% disp( ['Anharmonic correction to the density = ' num2str(anharmonicDensityCorrection(kB*T/U0))])
% introduce the anharmonic corretion for the density:
% tau_pTh = tau_pTh./anharmonicDensityCorrection(kB.*T./U0); %p-wave thermalization time [sec]

eta = U0./(kB.*T);
beta = 1./(kB.*T);
mu = m/2; % Reduced mass in a collision

if ~isempty(varargin) && strcmp(varargin{1}, 'numeric')
    %% numerical intgration of the density, not using the harmonic approximation
    
    % UODT = - 3*pi*c^2/(2*omega0^3)*( g/(omega0-omega) + g/(omega0+omega))*2*P/(pi*w[z]^2) * exp(-2*r^2/w^2[z]) = - U0/( 1 + (z/zR)^2 ) * exp(-2*r^2/w^2[z])
    % w^2[z] = w0^2*( 1 + (z/zR)^2 )
    
    %integration limits on r and z, up to 98% of the potential - to avoid integrating outside of the trap:
    % IntSizeR = 1.39857;
    % IntSizeZ = 7;
    %integration limits on r and z, up to 95%
    IntSizeR = 1.22387;
    IntSizeZ = 4.3589;
    % The results at 1/eta = 3 (around minimum and beyond) depend on the specific integration limit used (i.e. 95% or 98% changes the value).
    % cylindrical coordiantes:
    fun1 = @(r,z,T,U0) r.*exp(1./(kB*T) * U0./( 1 + (z./zR).^2 ) .* exp(-2.*r.^2./ (w0^2*( 1 + (z./zR).^2 )) ) );
    fun2 = @(r,z,T,U0) r.*exp(2./(kB*T) * U0./( 1 + (z./zR).^2 ) .* exp(-2.*r.^2./ (w0^2*( 1 + (z./zR).^2 )) ) );
    fun3 = @(r,z,T,U0) r.*exp(3./(kB*T) * U0./( 1 + (z./zR).^2 ) .* exp(-2.*r.^2./ (w0^2*( 1 + (z./zR).^2 )) ) );
    
    % q = integral3(fun,xmin,xmax,ymin,ymax,zmin,zmax)
    q1 = zeros(size(T));
    q2 = q1;
    q3 = q1;
    if numel(U0)==1
        U0 = U0.*ones(size(T)); %extend U0 to the size of T
    end
    for i = 1 : length(T) %calcualte the numerical density integrals for each temperature and trap depth
        q1(i) = integral2(@(r,z)fun1(r,z,T(i),U0(i)), 0,IntSizeR*w0,-IntSizeZ*zR,IntSizeZ*zR,'AbsTol', 1e-20,'RelTol',1e-20);
        q2(i) = integral2(@(r,z)fun2(r,z,T(i),U0(i)), 0,IntSizeR*w0,-IntSizeZ*zR,IntSizeZ*zR,'AbsTol', 1e-20,'RelTol',1e-20);
        q3(i) = integral2(@(r,z)fun3(r,z,T(i),U0(i)), 0,IntSizeR*w0,-IntSizeZ*zR,IntSizeZ*zR,'AbsTol', 1e-20,'RelTol',1e-20);
    end
    nmeanThNumeric = N./(2*pi) .* q2./q1.^2; %<n>
    n2meanThNumeric = N.^2./((2*pi)^2) .* q3./q1.^3; %<n^2>
    npeakThNumeric = N.*exp(beta.*U0)./(2*pi*q1); %Peak density
    % clear fun1 fun2 q1 q2
    % % cartesian coordiantes:
    % fun1 = @(x,y,z) exp(2./(kB*T) * U0./( 1 + (z./zR).^2 ) .* exp(-2*(x.^2+y.^2)./ (w0^2*( 1 + (z./zR).^2 )) ) );
    % fun2 = @(x,y,z) exp(1./(kB*T) * U0./( 1 + (z./zR).^2 ) .* exp(-2*(x.^2+y.^2)./ (w0^2*( 1 + (z./zR).^2 )) ) );
    % q1 = integral3(fun1,-IntSizeR*w0,IntSizeR*w0,-IntSizeR*w0,IntSizeR*w0,-IntSizeZ*zR,IntSizeZ*zR,'AbsTol', 0,'RelTol',1e-9); % q = integral3(fun,xmin,xmax,ymin,ymax,zmin,zmax)
    % q2 = integral3(fun2,-IntSizeR*w0,IntSizeR*w0,-IntSizeR*w0,IntSizeR*w0,-IntSizeZ*zR,IntSizeZ*zR,'AbsTol', 0,'RelTol',1e-9);
    % nmeanThNumericCart = N * q1./q2.^2;
    
    % Optical density calculation
%     for i = 1 : length(T) %calcualte the numerical density integrals for each temperature and trap depth
%         q11D(i) = integral1(@(r,z)fun1(r,z,T(i),U0(i)), 0,IntSizeR*w0,-IntSizeZ*zR,IntSizeZ*zR,'AbsTol', 1e-20,'RelTol',1e-20);
%     end
%     OD_r = 3*lambdaAtom.^2./(2*pi) * ; %Resonant sigma+ light optical density, along the radial direction
%     OD_z = 3*lambdaAtom.^2*N./( (2*pi)^2*sigma_r*sigma_r ); %Resonant sigma+ light optical density, along the axial direction
    
    % calculating the p-wave thermalization rate, using numerical density integral, and lower integration limit for the scattering amplitude
    scatInt = (96*V_p^2*mu^2)./(beta.^5*hbar^4) - (4.*exp(-eta)*V_p^2*mu^2.*(24 + 24.*eta + 12.*eta.^2 + 4.*eta.^3 + eta.^4))./(beta.^5*hbar^4);
    gamma_pThNumeric = 2./(4.1*kB.*T).*nmeanThNumeric.*sqrt(8/pi/mu).*(kB.*T).^(-3/2)*12*pi.*scatInt;
    tau_pThNumeric = 1./gamma_pThNumeric;
    snvNumeric = 2/4.1 * sigma_p .* nmeanThNumeric .* vrel;
    
    output.nmeanThNumeric = nmeanThNumeric;
    output.n2meanThNumeric = n2meanThNumeric ;
    output.tau_pThNumeric = tau_pThNumeric;
    output.snvNumeric = snvNumeric;
    output.npeakThNumeric = npeakThNumeric;
end
% disp(['eta = beta*U0 = ' num2str(U0/(kB*T)) '; 1/eta = ' num2str(kB*T./U0) ';'])
% disp(['Harmonic approximation: ' num2str(nmeanTh, '%10.4e\n')])
% disp(['Numeric density: ' num2str(nmeanThNumeric, '%10.4e\n')])
% disp(['Cartesian integration: ' num2str(nmeanThNumericCart, '%10.2e\n')])

%% assign results to a structure
output.P = P;
output.U0 = U0;
output.omega_r = omega_r;
output.omega_z = omega_z;
output.Tr = Tr;
output.Tz = Tz;
output.T = T;
output.TFarField = TFarField;
output.Tf = Tf;
output.sigma_r = sigma_r;
output.sigma_z = sigma_z;
output.nmeanTh = nmeanTh;
output.n2meanTh = n2meanTh;
output.K2Volume = K2Volume;
output.K3Volume = K3Volume;
output.tau_pTh = tau_pTh;
output.eta = eta;
output.snv = snv;
output.L2dipolar = L2dipolar;
output.npeakTh = npeakTh;
output.OD_r = OD_r;
output.OD_z = OD_z;

end

% plot mean density vs 1/eta:
% sizeX = 0.08:0.02:1.2;output = calcXODTDensity(1e6, 3.4, 0.5, sizeX, sizeX, 'numeric');
% figure; plot(1./output.eta, output.npeakThNumeric./output.npeakTh,'o')

