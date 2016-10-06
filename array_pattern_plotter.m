% Jianwei Sun
% 1000009821
function array_pattern_plotter(phase_delays, lambda, d, number_points, color)
% This function generates the magnitude plot of the array factor determined
% from an array of N antenna sources, each with a different phase
% difference from the first element
beta = 2*pi*d / lambda;

% phi is the angle between the observation point and the axis containing
% the point sources. We will sweep from 0 to 2*pi to get the full circle
phi = transpose(linspace(0, 2*pi, number_points));
N = length(phase_delays) + 1;
% The first phase term is 0
phase = zeros(number_points, 1);
array_factor = exp(1i*phase);
% Sum up the remaining point sources
for index = 2:N
    % phase is accumulated since it is measured relative to the first point
    % source
    phase = phase + beta*d*cos(phi) + phase_delays(index - 1);
    array_factor = array_factor + exp(1i*phase);
end
% The array_factor is normalized
array_factor = array_factor / N;
clearvars phase index;


% subplot(2,1,1);
% plot(phi, abs(array_factor));
% title('Normalized Array Factor');
% xlabel('Angle Phi (rads)');
% ylabel('|A(phi)|');
% grid on;
% xlim([0 2*pi]);
% set(gca,'Xtick',[0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4,2*pi]);
% 
% subplot(2,1,2);

plot(phi, abs(array_factor), color);

% polar(phi, abs(array_factor), color);
% title('Normalized Array Factor (Polar)');
% xlabel('Angle Phi (degrees)');
% ylabel('|A(phi)|');
% grid on;


end

