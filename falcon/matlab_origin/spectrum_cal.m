function [power_spectrum, omega] = spectrum_cal(field, xlamds, zsep)
%% FFT
ft = fftshift(fft(field));
power_spectrum = abs(ft).^2;
%% Calculate the frequency range
nslice = size(field,2);
omegas = 2*pi/(xlamds/2.99792458e8);
dt = nslice * zsep * xlamds / 2.99792458e8;
df = 2*pi/dt;
omega = df*(1:length(ft))';
omega = omega - median(omega) +omegas;
omega = omega/omegas;
%% Center the spetrum so you have delta omega/omega on x-axis
omega = omega - 1;
end