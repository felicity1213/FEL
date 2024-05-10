%%%%%%%%% This is edited for postprocess after calculation
close all
kw = 2*pi/parameter.lambdau;
%% Spectrum as a function of z (time-dependent)
zpos = [1:parameter.Nsnap]*parameter.stepsize;
zlocations = linspace(parameter.stepsize, parameter.lwig, 30);
fundpower = [];
sidebandpower = [];
zindices = round(zlocations/parameter.stepsize);
if flag.itdp
    omegamin = -5e-3;
    omegamax = 5e-3;
    h = figure(1);
    filename = 'radiation_movie.gif';
    for n = 1:1:length(zindices)
        [powerspec, omega] = spectrum_cal(radfield(zindices(n),:), parameter.lambda0, parameter.zsep);
        sidebandindex = omega>omegamin & omega< omegamax;
        fundspectrum = powerspec(sidebandindex);
        fundpower(n) = trapz(fundspectrum)/trapz(powerspec);
        figure(1)
        set(gcf, 'Color', 'w');
        subplot(1,2,1)
        semilogy(omega, abs(powerspec));
        xlabel('\delta\omega/\omega_0','FontSize',16);
        ylabel('P(\omega) [arb. units]','FontSize',16);
        xlim([-100,100].*rho1D);
        set(gca, 'FontSize',16);
        legend(sprintf(['z/L_u = ', num2str(zlocations(n)/parameter.lwig,'%.1f')]));
        legend boxoff
        subplot(1,2,2)
        plot([1:1:size(power,2)]*parameter.zsep, power(zindices(n),:));
        xlim([1,size(power,2)]*parameter.zsep);
        xlabel('ct/\lambda_0', 'FontSize',16);
        ylabel('Output Radiation Power [W]', 'FontSize',16);
        set(gca, 'FontSize',16)
        drawnow
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if n==1
            imwrite(imind,cm,filename,'gif','Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
    end
end
% sideband power fraction
if flag.itdp
    figure(2)
    set(gcf, 'Color', 'w');
    semilogy(zlocations, fundpower,'r')
    hold on
    semilogy(zlocations, 1-fundpower, 'b')
    legend('Fundamental Power(Fractional)','Sideband Power(Fractional)', 'location','best');
    legend boxoff
    ylim([0,1]);
    xlabel('z[m]')
    ylabel('P/P_{tolal}')
    enhance_plot
end



%% Radiation Power and spectrum at exit 
if ~flag.itdp    % (time-independent)
    Lgfit = fit_gainlength(Lgain, zpos, mean(power,2));
end
figure(3)
set(gcf, 'Units','Normalized', 'OuterPosition',[0,0,1,1], 'Color', 'w');
title('Simulation Output')
subplot(2,3,1)
semilogy(zpos./Lgain, mean(power,2)/parameter.Ee/parameter.I)
hold on
xlim([0,zpos(end)]/Lgain)
xlabel('z/L_g')
ylabel('P/P_{beam}')
title('Power Evolution')
legend(['P_{max} = ', num2str(max(mean(power,2)*1e-12), '%.2f'), 'TW'],'location','SouthEast');
legend boxoff
enhance_plot

if flag.itdp   % (time-depedent)
    subplot(2,3,2)
    plot([1:1:size(power,2)]*parameter.zsep, power(end,:))
    xlim([1,size(power,2)]*parameter.zsep)
    xlabel('ct/\lambda_0')
    ylabel('Power/[W]')
    title('Power along s')
    enhance_plot
    [powerspec, omega] = spectrum_cal(radfield(end,:),parameter.lambda0, parameter.zsep);
    subplot(2,3,3)
    semilogy((omega+1)*hbar*2*pi*c/parameter.lambda0, powerspec, 'b');
    xlim([omega(1)+1, omega(end)+1]*hbar*2*pi*c/parameter.lambda0);
    xlabel('Photon Energy [eV]');
    ylabel('P(\omega) [arb. units]');
    title('Outplot Spectrum');
    enhance_plot;
end

%% bunching factor and energy loss
figure(3)
subplot(2,3,4)
plot([1:1:parameter.Nsnap-1]*parameter.stepsize/Lgain, bunch);
xlim([0,parameter.Nsnap*parameter.stepsize]/Lgain);
xlabel('z/L_g');
ylabel('Bunching Factor')
title('Bunching Factor')
enhance_plot
for ij = 1:1:parameter.Nsnap
    meanenergy(ij) = mean(abs(mean(phasespace.gammap(ij,:,:),3)));
end
subplot(2,3,5)
plot([1:1:parameter.Nsnap]*parameter.stepsize/Lgain, (meanenergy/meanenergy(1)-1));
xlim([0, parameter.Nsnap]*parameter.stepsize/Lgain);
xlabel('z/L_g')
ylabel('\Delta\gamma/\gamma_0')
title('Energy Loss')
enhance_plot

%% tapering
if flag.tapering == 2
    const_resp=1/parameter.chi2*(parameter.lambdau/2/parameter.lambda0);
    meanfield = mean(radfield(:,:),2);
    dawdz = zeros(parameter.Nsnap, 1);
    dawdz(2:end) = abs(diff(parameter.awz)./parameter.delz./parameter.lambdau);
    psir = asin(const_resp.*(dawdz./abs(meanfield)));
    psirend = asin(const_resp.*(dawdz(end)./abs(radfield(end,:))));
    gammarofz = sqrt((parameter.lambdau/2/parameter.lambda0).*(1+parameter.awz.^2));
    dgammardz = -(parameter.lambdau/2/parameter.lambda0).*parameter.awz./gammarofz/2.*dawdz';
    
    options = optimset('Display','off');
    for i=1:length(psir)
        if psir(i)>0
            psr = psir(i);
            psi2(i) = pi-psr;
            psi1(i) = fsolve(@(x) cos(x)+x*sin(psr)-cos(psi2(i))-psi2(i)*sin(psr), -pi, options);
            dX = [psi1(i):pi/50:psi2(i)];
            Y = cos(psr)+cos(dX)-(pi-psr-dX)*sin(psr);
            % alpha(i)=sqrt(2)/8*trapz(dX,Y);
            alpha(i) = (1-sin(psr))/(1+sin(psr)); % an approx. by S.Y.Lee
        else
            psi(i)=pi;
            psi2(i)=-pi;
        end
    end
    figure(3)
    subplot(2,3,6)
    ksynch = parameter.ku*sqrt(cos(psir).*mean(parameter.chi2/parameter.k*abs(radfield),2)*parameter.awz/(1+parameter.awz.^2));
    semilogy(zpos./Lgain, 2*pi./ksynch)
    xlim([1,parameter.Nsnap]*parameter.stepsize./Lgain);
    xlabel('z/L_g');
    ylabel('\lambda_{synch}[m]');
    title('Synchronous Oscillation Wavelength')
    enhance_plot
    
    zoverlg = zpos./Lgain;
    figure(4)
    annotation('textbox',[0 0.9 1 0.1],...
        'String','Tapering Plots',...
        'EdgeColor','none',...
        'HorizontalAlignment','center',...
        'FontSize',30);
    set(gcf, 'Units','Normalized', 'OuterPosition',[0,0,1,1], 'Color','w');
    subplot(2,3,1)
    plot(zoverlg, psir*180/pi);
    xlabel('z/L_g');
    ylabel('\Psi_R [degree]');
    xlim([0,zoverlg(end)]);
    title('Resonant Phase')
    enhance_plot;
    subplot(2,3,2)
    plot(zoverlg,abs(meanfield)*1e-12);
    xlabel('z/L_g');
    ylabel('<E_{rad}> [TV/m]');
    xlim([0,zoverlg(end)]);
    title('Radfield')
    enhance_plot;
    subplot(2,3,3)
    plot(zoverlg,dgammardz*0.511);
    xlabel('z/L_g');
    ylabel('d\gamma/dz [MeV/m]');
    xlim([0,zoverlg(end)]);
    title('Energy descent rate')
    enhance_plot;
    if flag.itdp % time dependent resonant phase
        subplot(2,3,4)
        plot([1:1:size(radfield,2)]*parameter.zsep*parameter.lambda0*1e15/c, psirend*180/pi);
        xlim([1,size(radfield,2)]*parameter.zsep*parameter.lambda0*1e15/c);
        xlabel('t [fs]');
        ylabel('\Psi_R [degree]');
        title('Resonant Phase');
        enhance_plot;
    end
    subplot(2,3,6)
    plot(zoverlg,parameter.awz);
    xlabel('z/L_g');
    ylabel('RMS Undulator Parameter aw');
    title('RMS Undulator Parameter(aw)')
    xlim([0, zoverlg(end)]);
    enhance_plot;
end

%% Phasespace movie
zlocations = linspace(parameter.stepsize, parameter.lwig, 50);
zindices = round(zlocations/parameter.stepsize);
if flag.tapering == 2
    bucketheight = sqrt(parameter.chi2/parameter.ku.*parameter.awz./gammarofz.^2.*abs(mean(radfield(:,:),2))');
else
    bucketheight = sqrt(parameter.chi2/parameter.ku.*parameter.awz'./parameter.gamma0.^2.*abs(mean(radfield(:,:),2)));
end

if flag.phasespacemovie
    filename='particle_movie.gif';
    figure(5)
    set(gcf, 'Color','w')
    for i=1:1:length(zindices)
        x = linspace(-pi,pi,1e4);
        if flag.tapering==2
            sepa = bucketheight(zindices(i))'.*separatrix(x,res_phase(zindices(i)))+(gammarofz(zindices(i))-parameter.gamma0)/parameter.gamma0;
            sepamin = (gammarofz(zindices(i))-parameter.gamma0)/parameter.gamma0-bucketheight(zindices(i))'.*separatrix(x, res_phase(zindices(i)));
        else
            sepa = bucketheight(zindices(i))'.*separatrix(x, res_phase(zindices(i)));
            sepamin = -bucketheight(zindices(i))'.*separatrix(x,res_phase(zindices(i)));
            gammarofz = parameter.gamma0.*ones(length(zindices));
        end
        fieldphase(i)=mean(angle(radfield(zindices(i),:)));
        tp = squeeze(phasespace.thetap(zindices(i),:,:))+fieldphase(i)+pi/2;
        gp = squeeze(phasespace.gammap(zindices(i),:,:));
        tresh = reshape(tp, [1, size(tp,1)*size(tp,2)]);
        gresh = reshape(gp, [1, size(gp,1)*size(gp,2)]);
        % calculate trapping fraction
        [psi1, psi2] = bucket_parameter(parameter.psir);
        indi = (mod(tp,2*pi)-pi)>psi2 & (mod(tp, 2*pi)-pi)<psi1;
        g2 = (gp(indi)./(meanenergy(1))-1);
%         g2 = (gp(indi)./(parameter.gamma0)-1); % ~ the same as the above one
        angoli = (mod(tp(indi),2*pi)-pi);
        indi2 = g2<(bucketheight(zindices(i))'.*separatrix(angoli,parameter.psir)+(gammarofz(zindices(i))-parameter.gamma0)/parameter.gamma0);
        ftrap(i) = numel(g2(indi2))/numel(gp); % trapping fraction
        
        ind = x<psi1 & x>psi2;
        
        if flag.itdp
            tresh_mod=reshape(tp(1:1,:),[1,1*size(tp,2)]);
            gresh_mod=reshape(gp(1:1,:),[1,1*size(gp,2)]); % which slice do you want to draw
            plot((mod(tresh_mod,2*pi)-pi)./pi, (gresh_mod./(meanenergy(1))-1)*100,'.k','MarkerSize',1);
            drawnow
        else
            subplot(1,2,1);
            semilogy([1:zindices(i)]*parameter.stepsize/Lgain, power(1:zindices(i))/parameter.Ee/parameter.I,'k');
            xlabel('z/L_g');
            ylabel('P/P_{beam}');
            set(gca,'YTick',logspace(-4,-1,4));
            enhance_plot('Times',20);
            hold on;
            legend(['z/L_u = ',num2str(zlocations(i)/parameter.lwig)], 'location','SouthEast');
            legend boxoff
            subplot(1,2,2);
%             without the separatrix
%             plot((mod(tp,2*pi)-pi)./pi, (gp./(meanenergy(1))-1), '.k','MarkerSize',1)
            plot((mod(tp,2*pi)-pi)./pi, (gp./(meanenergy(1))-1)*100,'.k', x(ind)/pi,sepa(ind)*100,'r',x(ind)/pi,sepamin(ind)*100,'r','LineWidth',2);
            drawnow;
        end
        xlim([-1,1]);
        set(gca,'FontSize',20);
        xlabel('\Psi/pi');
        ylabel('\Delta \gamma/\gamma_0');
%         enhance_plot('FontSize',16);
%         if you want to save to GIF, uncomment the following line
        drawnow
        frame = getframe(5);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i==1
            imwrite(imind,cm,filename,'gif','Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
    end
    figure(4)
    subplot(2,3,5)
    plot(zlocations./parameter.lambdau, ftrap);
    xlabel('z/\lambda_u');
    ylabel('f_t[calculated]');
    title('Trapping Fraction')
    legend(sprintf(['f_t[Theory] = ', num2str((psi1-psi2)/2/pi)]))
    enhance_plot
end
            
            
        
        
        

























