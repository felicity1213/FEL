%%%%% the user's input parameters, similar to that in GENESIS 1.3
%%%%% these parameters are roughly divided into two parts as two structs.
%%%%% parameter.xxx & flag.xxx

%% Electron Beam Parameters
parameter.Ee = 13640.634e6;     % Electron beam energy [eV]
parameter.gamma0 = parameter.Ee/me/c^2*e0;     % Relativistic gamma factor
parameter.energyspread = 1.5e6;  % Absolute energy spread [eV]
parameter.deltagamma = parameter.energyspread/me/c^2*e0;
parameter.deltagammarel = parameter.deltagamma / parameter.gamma0;  % Relative energy spread dgamma/gamma
parameter.Np = 1024*1;            % # of macro electron in each slice
parameter.I = 4000;    % Peak current
parameter.emitx = 0.4e-6;  % Transverse emittance in x direction
parameter.emity = 0.4e-6;  % Transverse emittance in x direction
parameter.x0 = 0;  % initial transverse offset in x direction
parameter.y0 = 0;  % initial transverse offset in y direction
parameter.sigma_x = 1.2241e-5;  % beam radius in x direction (beam is assumed to be a circle)
parameter.sigma_y = 1.2241e-5;  % beam radius in y direction
parameter.betax = parameter.sigma_x^2*parameter.gamma0/parameter.emitx;  % beta fuction
parameter.betay = parameter.sigma_y^2*parameter.gamma0/parameter.emity;
parameter.A_e = 2*pi*parameter.sigma_x^2;     % beam cross section 

%% Radiation Parameters
parameter.lambda0 = 1.5*1e-10;           % radiation wavelength
parameter.k = 2*pi/parameter.lambda0;    % wavenumber in free space
P0 = 100; 
parameter.P0 = P0;                         % Seed power (W) 
parameter.zr = 5;                        % Rayleigh length of seed
parameter.waist = sqrt(parameter.zr*parameter.lambda0/pi);  % waist ofseed
parameter.A_mode = pi*parameter.waist^2/2;
parameter.E0 = sqrt(parameter.P0/c/eps0/parameter.A_mode);  % Assume circular polarization  

%% Undulator Parameters
parameter.lambdau = 3.0e-2;     % Undulator period
parameter.ku = 2.*pi./parameter.lambdau; % Undulator wavenumber
parameter.aw = sqrt(2*parameter.ku*parameter.gamma0^2/parameter.k-1);    % RMS undulator parameter
parameter.lwig = parameter.lambdau*2*0.5e3;  % Undulator length
% undulator type
flag.und = 0;  % 0 for helical, 1 for planar
if flag.und
    b = parameter.aw^2/(2+2*parameter.aw^2);
    parameter.JJ = besselj(0,b) - besselj(1,b);
else
    parameter.JJ = 1;
end
% taper option
flag.tapering = 0;  % tapering(0 for no taper, 1 for decelation rate , 2 for resonant, 3 for customization)
parameter.z0 = 10*parameter.lambdau;  % the place where the tapering starts
parameter.rate = 0.00018;
parameter.psir = pi/4; % resonant phase

%% Simulation Control
flag.phasespacemovie = 1; 
flag.itdp = 1;  % 0 for time-independent, 1 for time-dependent
flag.saveoutput = 1;
%
parameter.delz = 10; % Integration step size in measure of the undulator period length
parameter.zsep = parameter.delz; % Separation of beam slices in measures of radiation wavelength (n*delz)
parameter.stepsize = parameter.delz*parameter.lambdau; % Integration step size
parameter.Nsnap = round(parameter.lwig/parameter.stepsize);  % # of snapshots to take over the length of the undulator
flag.shotnoise = 1; % 1 -- with shotnoise, 0 -- without shotnoise
if(~flag.itdp)
    parameter.nslices = 1;
    flag.shotnoise =0;   % Note if you want to model time independent start-up from noise set P0 = pnoise
else
    parameter.nslices = round(8*parameter.Nsnap);  % Note you want more than 1 slippage length (Nsnap)
end

parameter.bunchlength = parameter.nslices*parameter.zsep*parameter.lambda0/c;

%% Simplifying constants
parameter.chi1 = mu0*c/2*parameter.I/parameter.A_e;
parameter.chi2 = e0/me/c^2;
% Constant for the resonant phase calculation   
parameter.const_resp = 1/parameter.chi2*(parameter.lambdau/2/parameter.lambda0);
