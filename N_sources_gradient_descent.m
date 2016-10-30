% Jianwei Sun 
% 1000009821
% This script illustrates the gradient descent method used to minimize a
% cost function associated with determining a target phase. For the case of
% only two radiating sources, the problem is easily solvable with gradient
% descent.

clc
clear all

lambda = 1;
% d is the distance between adjacent sources
d = lambda / 2;
beta = 2*pi*d / lambda;
% array of phase delays relative to the previous adjacent element. There
% are N-1 elements because the first delay is 0 since it is relative to
% itself
X = [-beta*d*cos(pi/3) -beta*d*cos(pi/3) -beta*d*cos(pi/3) -beta*d*cos(pi/3) -beta*d*cos(pi/3)];
X = zeros(1,15);
number_phases = length(X);
N = length(X) + 1;

%% Gradient descent method
iterations_limit = 10000;
cost_threshold = 10^-20;
% Initial phase condition = zeros
X_bar = zeros(number_phases, 1);
grad_X_bar = zeros(number_phases, 1);
learning_rate = 0.01;
rate_up = 1.01;
rate_down = 1.01;
% Number of data points is represented by variable D
D = 1000;
phi = linspace(0, 2*pi, D);

fprintf('Beginning pre-calculation for unchanging cos and sin term\n');
thetas = zeros(number_phases, D);
for i = 1:D
    thetas(:,i) = beta*d*cos(phi(i)) + X;
end

sum_thetas = zeros(number_phases, D);
sum_thetas(1,:) = thetas(1,:);
for i = 2:number_phases
    sum_thetas(i,:) = sum_thetas(i-1,:) + thetas(i,:);
end

exp_j_sum_thetas = exp(1j*sum_thetas);
exp_mj_sum_thetas = exp(-1j*sum_thetas);

fprintf('Pre-calculation complete. Beginning gradient descent\n');
for iteration = 1:iterations_limit
    thetas_bar = zeros(number_phases, D);
    for i = 1:D
        thetas_bar(:,i) = beta*d*cos(phi(i)) + X_bar;
    end  
    
    sum_thetas_bar = zeros(number_phases, D);
    sum_thetas_bar(1,:) = thetas_bar(1,:);
    for i = 2:number_phases
        sum_thetas_bar(i,:) = sum_thetas_bar(i-1,:) + thetas_bar(i,:);
    end

    exp_j_sum_thetas_bar = exp(1j*sum_thetas_bar);
    exp_mj_sum_thetas_bar = exp(-1j*sum_thetas_bar);
    
    B = zeros(1,D);
    B_bar = zeros(1,D);
    for i = 1:D
        for k = 1:number_phases
            B(i) = B(i) + exp_j_sum_thetas(k,i) - exp_j_sum_thetas_bar(k,i);
            B_bar(i) = B_bar(i) + exp_mj_sum_thetas(k,i) - exp_mj_sum_thetas_bar(k,i);
        end
    end
    B = B / N;
    B_bar = B_bar / N;
    
    d_B = zeros(number_phases,D);
    d_B_bar = zeros(number_phases,D);
    d_B(end,:) = exp_j_sum_thetas_bar(end,:);
    d_B_bar(end,:) = exp_mj_sum_thetas_bar(end,:);
    for i = (number_phases-1):-1:1
        d_B(i,:) = d_B(i + 1,:) + exp_j_sum_thetas_bar(i,:);
        d_B_bar(i,:) = d_B_bar(i + 1,:) + exp_mj_sum_thetas_bar(i,:);
    end
    d_B = d_B * (-1j/N);
    d_B_bar = d_B_bar * (1j/N);
    
    grad_X_bar = zeros(number_phases, 1);
    for i = 1:D
        grad_X_bar = grad_X_bar + d_B(:,i) .* B_bar(i) + B(i) .* d_B_bar(:,i);
    end
    
    X_bar = X_bar - learning_rate * grad_X_bar;
        
    iteration_cost = sumsqr(B.*B_bar);
%     if iteration > 1
%        if COST(iteration) > COST(iteration - 1)  %Becomes more expensive
%            learning_rate = learning_rate / rate_down;
%        end
%        if COST(iteration) < COST(iteration - 1)  %Heading in the right direction
%            learning_rate = learning_rate * rate_up;
%        end
%     end
    
    fprintf('Calculated %d iterations\n', iteration);
    if(iteration_cost < cost_threshold)
        break;
    end
end

figure;
array_pattern_plotter(X, lambda, d, D, 'ro');
hold on;
array_pattern_plotter(X_bar, lambda, d, D, 'b');
legend('Target radiation pattern','Pattern generated from gradient descent');

% figure;
% subplot(2,1,1);
% scatter(linspace(1, iterations, iterations), calculated_phases, 'r');
% hold on;
% plot([1 iterations], [phase_delay phase_delay]);
% title('Determining target phase');
% xlabel('Iteration number');
% ylabel('Phase');
% grid on;
% legend('Phase determined from gradient descent','Target phase');
% 
% subplot(2,1,2);
% plot(COST);
% title('Cost function per iteration');
% xlabel('Iteration number');
% ylabel('Cost');
% grid on;

