clc;
clear all;

%% Generating the phase matrix and effect matrix
lambda = 2;
d = lambda / 2;
% All distances are in terms of lambda
normal_feed_distance = 10 * lambda;

Sources = [0];
Reflect_elements = [-5/2 -4/2 -3/2 -2/2 -1/2 0 1/2 2/2 3/2 4/2 5/2];

fprintf('Generating phase matrix\n');
% Each row corresponds to a source
PHASE_MATRIX = zeros(length(Sources),length(Reflect_elements));

for i = 1:length(Sources)
    PHASE_MATRIX(i,:) = sqrt((Reflect_elements - Sources(i)).^2 + ones(1,length(Reflect_elements)).*normal_feed_distance^2).*(2*pi);
end

EFFECT_MATRIX = exp(PHASE_MATRIX .* 1j);
%% Optimizing
% C_k has to be the same length as Reflect_elements
C_k = exp(1j*zeros(length(Reflect_elements),1));

%% Plot results
resolution = 1000;
fprintf('Plotting results at resolution of %d points\n',resolution);
Phi = linspace(0,pi,resolution);
Theta = 2*pi*d^2/lambda .* cos(Phi);
Basis = exp(linspace(0,length(Reflect_elements)-1,length(Reflect_elements)).' * Theta .* 1j);
PLOT_ARRAY_FACTOR = zeros(length(Sources),resolution);
for i = 1:length(Sources)
    for j = 1:resolution
        PLOT_ARRAY_FACTOR(i,j) = (C_k.' .* EFFECT_MATRIX(i,:)) * Basis(:,j);
    end
end
PLOT_ARRAY_FACTOR = PLOT_ARRAY_FACTOR .* (1/length(Reflect_elements));

plot_colors = distinguishable_colors(length(Sources));
figure;
for i = 1:length(Sources)
    plot(Phi, abs(PLOT_ARRAY_FACTOR(i,:)), 'Color',plot_colors(i,:));
    hold on;
end
%axis([0 pi 0 1]);
set(gca,'xtick',0:pi/8:pi);
set(gca,'xticklabel',{'0','pi/8','pi/4','3 pi/8','pi/2','5 pi/6','3 pi/4', '7pi/8', 'pi'});
xlabel('Far-Field Angle (Radians)');
ylabel('Pattern Magnitude');
grid on;

% figure
% 
% %Calculates a set of colours for the variable number of motors to plot
% plot_colors = distinguishable_colors(quantity);
% 
% %Generate the legend labels
% plot_handle = zeros(quantity, 1);
% legend_entries = cell(quantity,1);
% for i = 1:quantity
%     kv_string = num2str(KV(i));
%     padding = repmat(' ', 1, 4 - length(kv_string));
%     legend_entries(i) = cellstr(['Kv = ', padding, kv_string]);
% end
% 
% %Plot the quantities
% subplot(3,1,1);
% for i = 1:quantity
%     plot_handle(i) = plot(F(i,1:end-1), P_IN_F(i,1:end-1), 'Color',plot_colors(i,:));
%     hold on;
%     plot(F(i,end), P_IN_F(i,end), 'Marker','o', 'MarkerSize',7, 'Color',plot_colors(i,:));
%     hold on;
% end
% xlabel('Force(N)');
% ylabel('Input Power(W)');
% grid on;
% title('Power - Force Curves for Varying Motors for a given Propeller');
% legend(plot_handle, legend_entries);
