%%%%% Simple FEL simulation code
%%%%% omit the transverse effect and use RK4 algorithm
clear
clc
%% Physical Constants
physical_constants
%% Input setting
FEL_sim_input
%% Calculate 1-D FEL parameters
FEL_parameter_cal
%% Calculate undulator field depending on taper
undulator_field

% figure(1)
% plot([1:1:parameter.Nsnap]*parameter.stepsize,res_phase)
% xlabel('z [m]');ylabel('\psi_r');enhance_plot

% figure(2)
% plot([1:1:parameter.Nsnap]*parameter.stepsize,parameter.awz)
% xlabel('z [m]');ylabel('aw');enhance_plot
%% Generate Initial Phasespace & Radfield
ini_phase

% figure(2)
% plot(squeeze(phasespace.thetap(1,1,:)),squeeze(phasespace.gammap(1,1,:)),'.','MarkerSize',10)
% xlabel('\psi');ylabel('\gamma');enhance_plot
% figure(3)
% plot(squeeze(phasespace.x(1,1,:)),squeeze(phasespace.y(1,1,:)),'.','MarkerSize',10)
% xlabel('x [m]');ylabel('y [m]');enhance_plot
%% Main Integration Routine
t0 = tic;
FEL_run_integration
disp(['Simulation time = ', num2str(toc(t0)./60), 'min'])
%% Postprocess
FEL_postprocess_output
if flag.saveoutput
    FEL_save_output
end



















