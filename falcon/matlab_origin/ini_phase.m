%%%%% initialize phase space, radfield
%%%%% hammersley.m
tic
disp('Loading Particles ...');
nbins = 64; % # of bins for each slices
mpart = parameter.Np/nbins; % # of electrons for each bins
n_electron = parameter.I*parameter.lambda0*parameter.zsep/e0/c; % 
p1 = zeros(parameter.Np,1);

phasespace.gammap = zeros(parameter.Nsnap,parameter.nslices,parameter.Np);
phasespace.thetap = zeros(parameter.Nsnap,parameter.nslices,parameter.Np);
% phasespace.x = zeros(parameter.Nsnap,parameter.nslices,parameter.Np);
% phasespace.y = zeros(parameter.Nsnap,parameter.nslices,parameter.Np);
% phasespace.px = zeros(parameter.Nsnap,parameter.nslices,parameter.Np);
% phasespace.py = zeros(parameter.Nsnap,parameter.nslices,parameter.Np);
radfield = parameter.E0*ones(parameter.Nsnap,parameter.nslices);
bunching = zeros(1,parameter.nslices);

for islice = 1:1:parameter.nslices
    Hammersley = hammersley(2, parameter.Np);
    % gammap
    % uniform distribution
    phasespace.gammap(1,islice,:) = parameter.gamma0+parameter.deltagamma*Hammersley(1,:);
    % gauss distribution
%     gam1 = sqrt(-2*log(Hammersley(1,:))).*cos(2*pi*Hammersley(2,:));
%     phasespace.gammap(1,islice,:) = parameter.gamma0+parameter.deltagamma*gam1;
    % thetap
    auxtheta1 = hammersley(1,mpart)'*2*pi/nbins-pi; % the phase (-pi,pi) divided by nbins, each bin has a width of 2pi/nbins
    for jbin = 1:nbins
        for ipart = 1:mpart
            phasespace.thetap(1,islice,ipart+(jbin-1)*mpart)=auxtheta1(ipart)+2*(jbin-1)*pi/nbins;
        end
    end
    if(flag.shotnoise) % a random offset to each macro particle phase to generate the correct statistic for the bunching factor.
        an = 2*sqrt(-log(rand(1))/n_electron);
        phin = rand(1)*2*pi;
        for ipart = 1:parameter.Np;
            phasespace.thetap(1,islice,ipart) = phasespace.thetap(1,islice,ipart)-an*sin(phasespace.thetap(1,islice,ipart)+phin);
        end
    end
    % x,y
    % gauss distribution
%     xx = sqrt(-2*log(Hammersley(3,:))).*cos(2*pi*Hammersley(4,:));
%     yy = sqrt(-2*log(Hammersley(3,:))).*sin(2*pi*Hammersley(4,:));
%     phasespace.x(1,islice,:) = parameter.x0+parameter.sigma_x*xx;
%     phasespace.y(1,islice,:) = parameter.y0+parameter.sigma_y*yy;
    % px,py
    % remained to be solved
    % prebunching
%     if (param.prebunching)
%         % Double buncher scheme a la N. Sudar http://www.sciencedirect.com/science/article/pii/S0168900217301924
%         [thetap(1,islice,:),gammap(1,islice,:)]=prebunch_particles(squeeze(thetap(1,islice,:)),squeeze(gammap(1,islice,:)),param);
%         % Single buncher scheme
%         %[thetap(1,islice,:),gammap(1,islice,:)]=single_prebuncher_particles(squeeze(thetap(1,islice,:)),squeeze(gammap(1,islice,:)),param);
%     end
    % bunching factor
    bunching(islice) = (sum(exp(1i.*phasespace.thetap(1,islice,:))/parameter.Np));
end
disp('Particle Loaded Successfully!')
toc


