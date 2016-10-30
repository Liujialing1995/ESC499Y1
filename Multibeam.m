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
rng(1);
resolution = 1000;
Phi = linspace(0,pi,resolution);
Theta = 2*pi*d/lambda .* cos(Phi);
Basis = exp(linspace(0,length(Reflect_elements)-1,length(Reflect_elements)).' * Theta .* 1j);
% Mask specification
MASK_L = zeros(1, resolution);
MASK_H = zeros(1, resolution);
MASK_H(1:end) = 0.25;
MASK_H(resolution/4 - floor(0.17/pi*resolution):resolution/4 + floor(0.17/pi*resolution)) = 1;
MASK_L(resolution/4 - floor(0.0298/pi*resolution):resolution/4 + floor(0.0298/pi*resolution)) = 0.9;

% C_k has to be the same length as Reflect_elements
C_k = exp(1j*zeros(length(Reflect_elements),1));
%C_k = exp(-1j*PHASE_MATRIX(1,:).');
multibeam_error_sumsqr_points_outside_mask(MASK_L, MASK_H, C_k.', ones(1,N), Basis);

% Set up parameters and configurations
generation = 0;
generation_limit = 250;
population_size = 100;
number_elites = 1;
sub_population_size = population_size - number_elites; % Elitism preserves
p_cross = 0.8;
p_mutation = 0.5;
number_of_mutation_stages = 1;

% Initializing population
P = uint16(zeros(N,population_size));  %Each column is a chromosome
for i = 1:population_size
    for j = 1:N
        P(j,i) = encode_gene(2*pi*rand);
    end
end
% Evaluate initial fitness
Fitness_history = evaluatePopulationFitness(P, MASK_L, MASK_H, EFFECT_MATRIX, Basis);
Best_fitness = max(Fitness_history);
% Evolution cycle
while 1
    if((Best_fitness(end) == 1) || (generation > generation_limit))
        break;
    end
    generation = generation + 1; 
    fprintf('Generation count: %d. ',generation);
    % Perform selection
    Sub_population = uint16(zeros(N, sub_population_size));
    for i = 1:floor(sub_population_size/2)*2;  % Ensures an even number for parents
       % Select N = 2 participants, and include the fitter one in the selection 
       ind_1 = ceil(rand*population_size);
       ind_2 = ind_1;
       while(ind_2 == ind_1)
           ind_2 = ceil(rand*population_size);
       end
       if(evaluateFitness(P(:,ind_1),MASK_L,MASK_H,EFFECT_MATRIX,Basis) >= evaluateFitness(P(:,ind_2),MASK_L,MASK_H,EFFECT_MATRIX,Basis))
           Sub_population(:,i) = P(:,ind_1);
       else
           Sub_population(:,i) = P(:,ind_2);
       end
    end
    
    % Perform crossover
    P_prime = uint16(zeros(N,sub_population_size));
    for i = 1:2:(sub_population_size-1)
        if(rand < p_cross)
           cross_index = ceil(N*rand);
           parent_1 = Sub_population(:,i);
           parent_2 = Sub_population(:,i+1);
           child_1 = parent_1(1:cross_index);
           child_1 = [child_1; parent_2(cross_index+1 : end)];
           child_2 = parent_2(1:cross_index);
           child_2 = [child_2; parent_1(cross_index+1 : end)];
           P_prime(:,i) = child_1;
           P_prime(:,i+1) = child_2;
        else
           % Children do not have cross_over properties
           P_prime(:,i) = Sub_population(:,i);
           P_prime(:,i+1) = Sub_population(:,i+1);
        end
    end
    P_prime(:,end) = Sub_population(:,end);
    
    %Perform mutation
    for i = 1:sub_population_size
        if(rand < p_mutation)
            for j = 1:number_of_mutation_stages
                mutation_index = ceil(N*rand);
                % Bitwise-xor with a 1 in the position of the mutation flips
                % the bit
                mask = bitshift(uint16(1),mutation_index-1);
                P_prime(:,i) = bitxor(P_prime(:,i),mask);
            end
        end
    end
    
    %Preserve elitism
    [Fittest, fittestIndices] = sort(Fitness_history(:,generation),'descend');
    P_prime = [P_prime P(:,fittestIndices(1:number_elites))];
    
    
    % Update variables
    Fitness_history = [Fitness_history evaluatePopulationFitness(P, MASK_L, MASK_H, EFFECT_MATRIX, Basis)];
    Best_fitness = [Best_fitness max(Fitness_history(:,generation))];
    fprintf('Best fitness: %d\n',Best_fitness(end));
    P = P_prime;
end
Latest_fitness = Fitness_history(:,generation);
max_index = find(Latest_fitness == max(Latest_fitness),1);
C_k_bar = exp(1j*arrayfun(@decode_gene,P(:,max_index)));
if(Best_fitness(end) == 1)
    fprintf('\nGlobal optimum found\n');
end
if(generation > generation_limit)
    fprintf('\nGeneration limit of %d generations reached\n',generation_limit);
    fprintf('Latest fitness score: %d\n\n',max(Latest_fitness));
end


%% Plot results

fprintf('Plotting results at resolution of %d points\n',resolution);
PLOT_ARRAY_FACTOR = zeros(length(Sources),resolution);
for i = 1:length(Sources)
    PLOT_ARRAY_FACTOR(i,:) = (C_k_bar.' .* EFFECT_MATRIX(i,:)) * Basis;
    %PLOT_ARRAY_FACTOR(i,:) = (C_k_bar.' ) * Basis;
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

% Plot errors
figure;
plot(Best_fitness);
grid on;

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


