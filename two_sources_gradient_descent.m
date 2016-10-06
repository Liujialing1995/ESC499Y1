% Jianwei Sun 
% 1000009821
% This script illustrates the gradient descent method used to minimize a
% cost function associated with determining a target phase. For the case of
% only two radiating sources, the problem is easily solvable with gradient
% descent.

lambda = 1;
% d is the distance between adjacent sources
d = lambda / 2;
beta = 2*pi*d / lambda;
% array of phase delays relative to the previous adjacent element. There
% are N-1 elements because the first delay is 0 since it is relative to
% itself
phase_delay = -beta*d*cos(pi/3);

%% Gradient descent method
iterations = 10;
gradients = 0;
% Initial phase = 0
phase = 0;
learning_rate = 0.1;
rate_up = 0;
rate_down = 0;

calculated_phases = zeros(iterations, 1);
COST = zeros(iterations, 1);
for iteration = 1:iterations
    gradients = 2*(cos(phase_delay)-cos(phase))*sin(phase) - 2*(sin(phase_delay)-sin(phase))*cos(phase);
    phase = phase - learning_rate * gradients;
    calculated_phases(iteration) = phase;
    
    COST(iteration) = sumsqr([(cos(phase_delay)-cos(phase)) (sin(phase_delay)-sin(phase))]);
    if iteration > 1
       if COST(iteration) > COST(iteration - 1)  %Becomes more expensive
           learning_rate = learning_rate - rate_down;
       end
       if COST(iteration) < COST(iteration - 1)  %Heading in the right direction
           learning_rate = learning_rate + rate_up;
       end
    end
end

figure;
subplot(2,1,1);
scatter(linspace(1, iterations, iterations), calculated_phases, 'r');
hold on;
plot([1 iterations], [phase_delay phase_delay]);
title('Determining target phase');
xlabel('Iteration number');
ylabel('Phase');
grid on;
legend('Phase determined from gradient descent','Target phase');

subplot(2,1,2);
plot(COST);
title('Cost function per iteration');
xlabel('Iteration number');
ylabel('Cost');
grid on;

