clear
close all
clc

% Create signal with three different frequency components
x = linspace(-5,5,501);
freq1 = 1; freq2 = 3; freq3 = 6;
f1 = sin(2*pi*freq1*x) + sin(2*pi*freq2*x) + sin(2*pi*freq3*x);
f2 = exp(-1i*2*pi*freq1*x) + exp(-1i*2*pi*freq2*x) + exp(-1i*2*pi*freq3*x);
f3 = x;

% Plot signal
figure
hold on
% plotFFT(x,f1)
% plotFFT(x,f2)
plotFFT(x,f3)
hold off