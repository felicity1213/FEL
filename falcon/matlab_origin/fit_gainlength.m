function Lgfit = fit_gainlength(Lgaintheory, z1, pow)
%%% this function is based on the L_gain defined by power
%%% P(z) = P(0)*exp((z-0)/L_gain) 
zmin = 5*Lgaintheory;
zmax = zmin + 10*Lgaintheory;
pgain = pow(z1>zmin & z1<zmax);
zgain = z1(z1>zmin & z1<zmax);

Lgfit = (zgain(end)-zgain(1))/(log(pgain(end)/pgain(1)));
fprintf(['Lg_theory = ', num2str(Lgaintheory), '\n']);
fprintf(['Lg_sim = ', num2str(Lgfit),'\n']);
% figure(2)
% semilogy(linspace(zmin, zmax, length(pgain)), pgain);
% hold on
% semilogy(linspace(zmin, zmax), pgain(1)*exp((linspace(zmin,zmax)-zmin)./Lgfit),'r--')
