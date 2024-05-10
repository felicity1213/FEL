%%%%% This code is used to calculate the integration step based on
%%%%% fundamental FEL differential equation. IN order to be in accordance
%%%%% with KMR style 1-D taper where you compute the undulator aw
%%%%% step-by-step by specifying a constant resonant phase according to the
%%%%% equation of electron energy evolution and resonance condition.
%%%%% The phasespace evolution is calculated by RK4 algorithm and radfield is
%%%%% calculated by 'predict-correct' method.
total_simtime = 0;
%% Time dependent case
if flag.itdp
    for ij = 1:1:(parameter.Nsnap-1) % takes Nsnap snapshots along length of undulator
        tstart = tic;
        firstslice = ij;
        for islice = firstslice:1:parameter.nslices
            [newphasespace, radfield(ij+1,islice)] = FEL_RK4([squeeze(phasespace.thetap(ij,islice,:)),squeeze(phasespace.gammap(ij,islice,:))],radfield(ij,islice),parameter,parameter.awz(ij));
            phasespace.thetap(ij+1,islice,:) = newphasespace(:,1);
            phasespace.gammap(ij+1,islice,:) = newphasespace(:,2);
        end
        % Slippage of radiation field
        radfield(ij+1,:) = circshift(radfield(ij+1,:).',1);
        radfield(ij+1,1:(ij-1)) = 0; % Set field which slips into the calculation window to zeros
        % Compute bunching
        bunch(ij) = mean(abs(mean(exp(1j.*phasespace.thetap(ij,:,:)),3)));
        % Compute undulator field at next step (constant resonance phase)
        parameter.awz(ij+1) = parameter.awz(ij)-parameter.stepsize/parameter.const_resp*abs(mean(radfield(ij,parameter.Nsnap:parameter.nslices),2)).*sin(res_phase(ij));
        
        % disp calculation information
        dispformat = '%.3f sec from z= %.3f to z= %.3f, total length: %.3f \n';
        fprintf(dispformat, toc(tstart), parameter.stepsize*(ij-1), parameter.stepsize*ij, parameter.Nsnap*parameter.stepsize);
    end
end
%% Time-indepent case
if ~flag.itdp
    for ij = 1:1:(parameter.Nsnap-1) % takes Nsnap snapshots along length of undulator
        [newphasespace, radfield(ij+1,1)] = FEL_RK4([squeeze(phasespace.thetap(ij,1,:)),squeeze(phasespace.gammap(ij,1,:))],radfield(ij,1),parameter,parameter.awz(ij));
        phasespace.thetap(ij+1,1,:) = newphasespace(:,1);
        phasespace.gammap(ij+1,1,:) = newphasespace(:,2);
        % Compute bunching
        bunch(ij) = mean(abs(mean(exp(1j.*phasespace.thetap(ij,:,:)),3)));
        % Compute undulator field at next step (constant resonance phase)
        parameter.awz(ij+1) = parameter.awz(ij)-parameter.stepsize/parameter.const_resp*mean(abs(radfield(ij,:)),2).*sin(res_phase(ij));
    end
end
%% Remove slices within one total slippage length for time-dependent case
if flag.itdp
    radfield(:,1:parameter.Nsnap) = [];
    phasespace.gammap(:,1:parameter.Nsnap,:) = [];
    phasespace.thetap(:,1:parameter.Nsnap,:) = [];
end
%% Calculation radiation power
power(:,:) = abs(radfield(:,:)).^2/377*parameter.A_e;
        
        
% power(:, :) = abs(radfield_save(:, :)).^2/377/2*A_e;




        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
