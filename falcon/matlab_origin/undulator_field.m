%%%%% this code is used to generate the undulator field -- aw depending on
%%%%% the taper.
%%%%% flag.tapering = 0 for no taper, aw(i) = parameter.aw
%%%%% 1 for decelation rate , 
%%%%% 2 for resonant, 
%%%%% 3 for customization

startindex=floor(parameter.z0/parameter.stepsize); % the index of the taper begin
if flag.tapering == 0
    parameter.awz = parameter.aw*ones(1,parameter.Nsnap);
    res_phase=zeros(1,parameter.Nsnap);
end
if flag.tapering == 1
    parameter.awz(1:startindex)=parameter.aw;
    ind = 0:1:(parameter.Nsnap-startindex);
    parameter.awz(startindex:parameter.Nsnap)=parameter.aw*(1-parameter.rate*(ind.*parameter.stepsize).^2);
end
if flag.tapering == 2
    res_phase(1:startindex)=0;
    res_phase(startindex:parameter.Nsnap)=parameter.psir;
    parameter.awz(1)=parameter.aw;
end
if flag.tapering ==3;
    parameter.awz=[0];
end

    





% if param.tapering==0
%     res_phase=zeros(1,param.Nsnap);
% else
%     res_phase(1:startindex)=0;
%     res_phase(startindex:param.Nsnap)=param.psir;
% end
% Kz(1)=param.K;
% plot([1:1:parameter.Nsnap]*parameter.stepsize,parameter.awz)
% xlabel('z [m]');ylabel('aw');enhance_plot
% figure(1)
% plot([1:1:parameter.Nsnap]*parameter.stepsize,res_phase)
% xlabel('z [m]');ylabel('\psi_r');enhance_plot
