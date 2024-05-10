% calculate FEL parameters
rho1D = 1/parameter.gamma0*(1/8*parameter.I/IA*parameter.aw.^2/parameter.sigma_x^2/parameter.ku^2)^(1/3);
Lgain = parameter.lambdau/(4*sqrt(3)*pi*rho1D);
pnoise=parameter.gamma0*0.511e6*2*pi*3e8/parameter.lambda0*rho1D^2/2*1.6e-19;
Lsat = parameter.lambdau/rho1D;
Psat = 1.6*rho1D*parameter.Ee*parameter.I;