%%%%% this code is edited to save the figure and .mat for FEL simulation
%%%%% .../Simulation_Results/timestr/mat  
%%%%% .../Simulation_Results/timestr/figure
%%%%% .../Simulation_Results/timestr/animation

%% create folder
currentdir = pwd;
time = datestr(now,31);
matpath = fullfile(currentdir,'Simulation_Results',time,'mat');mkdir(matpath);
figpath = fullfile(currentdir,'Simulation_Results',time,'figure');mkdir(figpath);
anipath = fullfile(currentdir,'Simulation_Results',time,'animation');mkdir(anipath);
%% .mat
% save(fullfile(matpath,'Workspace'))
averagepower = mean(power,2);
awz = parameter.awz;
save(fullfile(matpath,'average_power'),'averagepower')
save(fullfile(matpath,'average_bunching'),'bunch')
save(fullfile(matpath,'undulator_parameter(aw)'),'awz')
save(fullfile(matpath,'parameter'),'parameter')
if flag.itdp
    save(fullfile(matpath,'output_sqectrum'),'powerspec')
    save(fullfile(matpath,'output_frequency'),'omega')
    output_field = radfield(end,:);
    output_power = power(end,:);
    output_time = [1:1:size(power,2)]*parameter.zsep*parameter.lambda0;
    save(fullfile(matpath,'output_field'),'output_field')
    save(fullfile(matpath,'output_power'),'output_power')
    save(fullfile(matpath,'output_time'),'output_time')
end
if ~flag.itdp
    thetap=phasespace.thetap;
    gammap=phasespace.gammap;
    save(fullfile(matpath,'thetap'),'thetap')
    save(fullfile(matpath,'gammap'),'gammap')
    save(fullfile(matpath,'radfield'),'radfield')
end

%% figure
figure(2)
export_fig(fullfile(figpath,'Spectrum_fraction'),'-pdf','-r600')
figure(3)
export_fig(fullfile(figpath,'Simulation_result'),'-pdf','-r600')
figure(4)
export_fig(fullfile(figpath,'Tapering_plot'),'-pdf','-r600')

%% animation
copyfile('radiation_movie.gif',anipath);
copyfile('particle_movie.gif',anipath);
delete('radiation_movie.gif')
delete('particle_movie.gif')










