

clc
clear all

lambda = 2;
% d is the distance between adjacent sources
d = lambda / 2;
beta = 2*pi*d / lambda;

number_samples = 1000;

%phi = linspace(0, 2*pi, number_samples);
%theta = beta*d*cos(phi);
theta = linspace(-beta*d,beta*d, number_samples);

%series_coefficients = ones(N,1);
% series_coefficients = zeros(N,1);
series_coefficients = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];
direction = pi/2;
% series_coefficients = [1; ...
%                        exp(1j*-beta*d*cos(direction)); ...
%                        exp(1j*2*-beta*d*cos(direction)); ...
%                        exp(1j*3*-beta*d*cos(direction)); ...
%                        exp(1j*4*-beta*d*cos(direction)); ...
%                        exp(1j*5*-beta*d*cos(direction))];
                   
N = length(series_coefficients);

%Input requirements for the MASK
MASK_L = zeros(1, number_samples);
MASK_H = zeros(1, number_samples);


MASK_H(1:end) = 0.3;
MASK_H(number_samples/2+number_samples/20:number_samples/2+number_samples/5+number_samples/8-number_samples/20) = 1;
MASK_L(number_samples/2 + number_samples/8:number_samples/2 + number_samples/10 + number_samples/10) = 0.9;

MASK_H = fliplr(MASK_H);
MASK_L = fliplr(MASK_L);

A = zeros(1,number_samples);
for i = 0:(N-1)
    % Computes the array factor for each sample point
    %A(i) = (1/N)*transpose(series_coefficients)*exp(1j*theta(i)*transpose(linspace(0,N-1,N)));
    A = A + series_coefficients(i+1)*exp(1j*i*theta);
end
A = A/N;

figure;
plot(fliplr(acos(theta/(beta*d))),fliplr(abs(A)),'r');
hold on;
plot(fliplr(acos(theta/(beta*d))),fliplr(MASK_H),'k');
hold on;
plot(fliplr(acos(theta/(beta*d))),fliplr(MASK_L),'k');
hold on;
axis([0 pi 0 1]);

%% Method to iteratively generate the array factor

iterations_limit = 2000;
iterations_count = 1;
acceptance_threshold = 0.005;
valid_yet = false;

all_errors_1 = zeros(1,iterations_limit);
all_errors_2 = zeros(1,iterations_limit);

Derived_AF = zeros(1,number_samples);
while(~valid_yet&&(iterations_count < iterations_limit))
    Derived_AF = zeros(1,number_samples);
    for i = 0:(N-1)
        % Computes the array factor for each sample point
        Derived_AF = Derived_AF + series_coefficients(i+1)*exp(1j*i*theta);
    end
    Derived_AF = Derived_AF/N;
    
    %Check for satisfying mask requirements
    if((any(abs(Derived_AF) > MASK_H) == 0)&&(any(abs(Derived_AF) < MASK_L) == 0))
        valid_yet = true;
    end
    
    error_1 = sumsqr(abs(Derived_AF(abs(Derived_AF) < MASK_L) - MASK_L(abs(Derived_AF) < MASK_L))) + ...
            sumsqr(abs(Derived_AF(abs(Derived_AF) > MASK_H) - MASK_H(abs(Derived_AF) > MASK_H)));
      

    Basis = exp(linspace(0,N-1,N).' * theta .* 1j);
    error_1 = multibeam_error_sumsqr_points_outside_mask(MASK_L, MASK_H, series_coefficients.', ones(1,N), Basis);

    error_2 = (sumsqr(abs(Derived_AF(abs(Derived_AF) < MASK_L)).^2 - abs(MASK_L(abs(Derived_AF) < MASK_L)).^2) + ...
            sumsqr(abs(Derived_AF(abs(Derived_AF) > MASK_H)).^2 - abs(MASK_H(abs(Derived_AF) > MASK_H)).^2));
        
    all_errors_1(iterations_count) = error_1;
    all_errors_2(iterations_count) = error_2;
    if(error_2 < acceptance_threshold)
        break;
    end
%     if(iterations_count > 1 && (error > all_errors(iterations_count - 1)))
%         break;
%     end
    
    %Doesn't satisfy, then compute the series coefficients
    Derived_AF(abs(Derived_AF) > MASK_H) = MASK_H(abs(Derived_AF) > MASK_H);
    Derived_AF(abs(Derived_AF) < MASK_L) = MASK_L(abs(Derived_AF) < MASK_L);
    % Map to N series and plot
    for i = 0:(N-1)
        integrand = Derived_AF.*exp(-1j*i*theta);
        integral_sum = 0;
        for m = 1:number_samples-1
            integral_sum = integral_sum + beta*d/number_samples*(integrand(m)+integrand(m+1));
        end
        series_coefficients(i+1) = integral_sum*N/(2*beta*d);
    end
    
    
    iterations_count = iterations_count + 1;
    if(mod(iterations_count,100)==0)
        fprintf('Elapsed iterations: %d\n',iterations_count);
    end
   
    
end
%%

fprintf('Number of iterations taken: %d\n',iterations_count);
fprintf('Final error_1: %d\n',error_1);
fprintf('Final error_2: %d\n',error_2);

plot(fliplr(acos(theta/(beta*d))),fliplr(abs(Derived_AF)),'b');

figure;
subplot(2,1,1);
plot(all_errors_1(1:iterations_count));
title('Error 1');
subplot(2,1,2);
plot(all_errors_2(1:iterations_count));
title('Error 2');

for i = 1:iterations_limit - 1
    if(all_errors_1(i+1) > all_errors_1(i))
        fprintf('Index %d has increasing error_1\n', i);
        break;
    end
end
for i = 1:iterations_limit - 1
    if(all_errors_2(i+1) > all_errors_2(i))
        fprintf('Index %d has increasing error_2\n', i);
        break;
    end
end