%Pattern plotter
clc;
clear all;

%Parameters
lambda = 2;
d = lambda/2;
beta_constant = 2*pi/lambda;
resolution = 1000;
N = 16;
Phases = zeros(1,N); %The phases are absolute phases
s = -beta_constant*d*cos(pi/4);
Phases = linspace(0,(N-1)*s,N);
N = length(Phases);
C_k = exp(1j*Phases);
PhaseMatrix = [628.4442  628.4147  628.3892  628.3676  628.3499  628.3362  628.3264  628.3205  628.3185  628.3205  628.3264  628.3362  628.3499  628.3676  628.3892  628.4147];
EffectMatrix = exp(1j*PhaseMatrix);

%Data structures
Phi = linspace(0,pi,resolution);
Theta = beta_constant*d .* cos(Phi);
Basis = exp(linspace(0,N-1,N).' * Theta .* 1j);
AF = (C_k .* EffectMatrix) * Basis ./N;

% Mask specification
MASK_L = zeros(1, resolution);
MASK_H = zeros(1, resolution);
MASK_H(1:end) = 0.25;
MASK_H(resolution/4 - floor(0.17/pi*resolution):resolution/4 + floor(0.17/pi*resolution)) = 1;
MASK_L(resolution/4 - floor(0.0298/pi*resolution):resolution/4 + floor(0.0298/pi*resolution)) = 0.9;

%Plotting
figure;
plot(Phi,abs(AF));
hold on;
plot(Phi,MASK_H,'k');
hold on;
plot(Phi,MASK_L,'k');
grid on;
axis([0 pi 0 1]);
set(gca,'xtick',0:pi/8:pi);
set(gca,'xticklabel',{'0','pi/8','pi/4','3 pi/8','pi/2','5 pi/6','3 pi/4', '7pi/8', 'pi'});
xlabel('Far-Field Angle (Radians)');
ylabel('Pattern Magnitude');
title('Simple Pattern Plotter');


