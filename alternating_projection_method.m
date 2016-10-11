

clc
clear all

lambda = 2;
% d is the distance between adjacent sources
d = lambda / 2;
beta = 2*pi*d / lambda;

N = 6;

number_samples = 10000;

%phi = linspace(0, 2*pi, number_samples);
%theta = beta*d*cos(phi);
theta = linspace(-beta*d,beta*d, number_samples);

%series_coefficients = ones(N,1);
% series_coefficients = zeros(N,1);
series_coefficients = [1;1;1;1;1;1];
direction = pi/2;
series_coefficients = [1; ...
                       exp(1j*-beta*d*cos(direction)); ...
                       exp(1j*2*-beta*d*cos(direction)); ...
                       exp(1j*3*-beta*d*cos(direction)); ...
                       exp(1j*4*-beta*d*cos(direction)); ...
                       exp(1j*5*-beta*d*cos(direction))];

%Input requirements for the MASK
MASK_L = zeros(1, number_samples);
MASK_H = zeros(1, number_samples);


MASK_H(1:number_samples/2) = 0.3;
MASK_H(number_samples/2:end) = 1;
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

%all_errors = zeros(1,iterations_limit);

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
    
    error = sumsqr(abs(Derived_AF(abs(Derived_AF) < MASK_L) - MASK_L(abs(Derived_AF) < MASK_L))) + ...
            sumsqr(abs(Derived_AF(abs(Derived_AF) > MASK_H) - MASK_H(abs(Derived_AF) > MASK_H)));
    if(error < acceptance_threshold)
        break;
    end
    %all_errors(iterations_count) = error;
    
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
% figure;
% plot(theta,abs(A),'r');
% hold on;

% % Map to N series and plot
% for i = 0:(N-1)
%     integrand = A.*exp(-1j*i*theta);
%     integral_sum = 0;
%     for m = 1:number_samples-1
%         integral_sum = integral_sum + beta*d/number_samples*(integrand(m)+integrand(m+1));
%     end
%     series_coefficients(i+1) = integral_sum*N/(2*beta*d);
% end
% Derived_AF = zeros(1,number_samples);
% for i = 0:(N-1)
%     % Computes the array factor for each sample point
%     %A(i) = (1/N)*transpose(series_coefficients)*exp(1j*theta(i)*transpose(linspace(0,N-1,N)));
%     Derived_AF = Derived_AF + series_coefficients(i+1)*exp(1j*i*theta);
% end
% Derived_AF = Derived_AF/N;

% % Map to any number of series and plot
% N_augmented = 2*N;
% derived_coeffs = zeros(N_augmented,1);
% for i = 0:(N_augmented-1)
%     integrand = A.*exp(-1j*i*theta);
%     integral_sum = 0;
%     for m = 1:number_samples-1
%         integral_sum = integral_sum + beta*d/number_samples*(integrand(m)+integrand(m+1));
%     end
%     derived_coeffs(i+1) = integral_sum*N_augmented/(2*beta*d);
% end
% Derived_AF = zeros(1, number_samples);
% for i = 0:(N_augmented - 1)
%     Derived_AF = Derived_AF + derived_coeffs(i+1)*exp(1j*i*theta);
% end
% Derived_AF = Derived_AF/N_augmented;

% fprintf('Error is %d\n',sumsqr(abs(A-Derived_AF)));
fprintf('Number of iterations taken: %d\n',iterations_count);

plot(fliplr(acos(theta/(beta*d))),fliplr(abs(Derived_AF)),'b');