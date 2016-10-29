clc;
clear all;

%% Generating the phase matrix and effect matrix
lambda = 2;
d = lambda / 2;
beta = 2*pi/lambda;
% All distances are in terms of lambda
normal_feed_distance = 100 * lambda;

Sources = [0];
Reflect_elements = linspace(-8/2,7/2,16);
%Reflect_elements = linspace(-100/2, 100/2, 201);

N = length(Reflect_elements);
M = length(Sources);

fprintf('Generating phase matrix\n');
% Each row corresponds to a source
PHASE_MATRIX = zeros(length(Sources),length(Reflect_elements));

for i = 1:length(Sources)
    PHASE_MATRIX(i,:) = sqrt((Reflect_elements - Sources(i)).^2 + ones(1,length(Reflect_elements)).*normal_feed_distance^2).*(2*pi)/lambda;
end

EFFECT_MATRIX = exp(PHASE_MATRIX .* 1j);
%% Optimizing
resolution = 1000;
Phi = linspace(0,pi,resolution);
Theta = 2*pi*d/lambda .* cos(Phi);
Basis = exp(linspace(0,length(Reflect_elements)-1,length(Reflect_elements)).' * Theta .* 1j);
% Mask specification
MASK_L = zeros(1, resolution);
MASK_H = zeros(1, resolution);
MASK_H(1:end) = 0.3;
MASK_H(resolution/2+resolution/20:resolution/2+resolution/5+resolution/8-resolution/20) = 1;
MASK_L(resolution/2 + resolution/8:resolution/2 + resolution/10 + resolution/10) = 0.9;
MASK_H = fliplr(MASK_H);
MASK_L = fliplr(MASK_L);

% C_k has to be the same length as Reflect_elements
C_k = exp(1j*zeros(length(Reflect_elements),1));
%C_k = exp(-1j*PHASE_MATRIX(1,:).');
multibeam_error_sumsqr_points_outside_mask(MASK_L, MASK_H, C_k.', ones(1,N), Basis);


%% Plot results
fprintf('Plotting results at resolution of %d points\n',resolution);
PLOT_ARRAY_FACTOR = zeros(length(Sources),resolution);
for i = 1:length(Sources)
    PLOT_ARRAY_FACTOR(i,:) = (C_k.' .* EFFECT_MATRIX(i,:)) * Basis;
    %PLOT_ARRAY_FACTOR(i,:) = (C_k.' ) * Basis;
end
PLOT_ARRAY_FACTOR = PLOT_ARRAY_FACTOR .* (1/length(Reflect_elements));

plot_colors = distinguishable_colors(length(Sources));
figure;
for i = 1:length(Sources)
    plot(Phi, abs(PLOT_ARRAY_FACTOR(i,:)), 'Color',plot_colors(i,:));
    hold on;
end

plot(fliplr(acos(Theta/(beta*d))),fliplr(MASK_H),'k');
hold on;
plot(fliplr(acos(Theta/(beta*d))),fliplr(MASK_L),'k');
hold on;
axis([0 pi 0 1]);
set(gca,'xtick',0:pi/8:pi);
set(gca,'xticklabel',{'0','pi/8','pi/4','3 pi/8','pi/2','5 pi/6','3 pi/4', '7pi/8', 'pi'});
xlabel('Far-Field Angle (Radians)');
ylabel('Pattern Magnitude');
grid on;
title('Multibeam Pattern Synthesis');

% figure;
% for i = 1:length(Sources)
%     plot(Phi, 0.5*mag2db(abs(PLOT_ARRAY_FACTOR(i,:))), 'Color',plot_colors(i,:));
%     hold on;
% end
% axis([0 pi -inf inf]);
% set(gca,'xtick',0:pi/8:pi);
% set(gca,'xticklabel',{'0','pi/8','pi/4','3 pi/8','pi/2','5 pi/6','3 pi/4', '7pi/8', 'pi'});
% xlabel('Far-Field Angle (Radians)');
% ylabel('Pattern Power Magnitude (dB)');
% grid on;
% title('Multibeam Pattern Synthesis');


