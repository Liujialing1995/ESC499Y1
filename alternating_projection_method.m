

clc
clear all

lambda = 2;
% d is the distance between adjacent sources
d = lambda / 2;
beta = 2*pi*d / lambda;

N = 5;

number_samples = 1000;

%phi = linspace(0, 2*pi, number_samples);
%theta = beta*d*cos(phi);
theta = linspace(-beta*d,beta*d, number_samples);

%series_coefficients = ones(N,1);
series_coefficients = zeros(N,1);
series_coefficients(1) = 1; 
series_coefficients(2) = exp(1j*pi/3); 
series_coefficients(3) = exp(1j*2*pi/3);
series_coefficients(4) = exp(1j*3*pi/3); 
series_coefficients(5) = exp(1j*4*pi/3);


A = zeros(1,number_samples);
for i = 0:(N-1)
    % Computes the array factor for each sample point
    %A(i) = (1/N)*transpose(series_coefficients)*exp(1j*theta(i)*transpose(linspace(0,N-1,N)));
    A = A + series_coefficients(i+1)*exp(1j*i*theta);
end
A = A/N;

figure;
plot(theta,abs(A),'r');
hold on;

% % Map to N series and plot
% for i = 0:(N-1)
%     integrand = A.*exp(-1j*i*theta);
%     integral_sum = 0;
%     for m = 1:number_samples-1
%         integral_sum = integral_sum + beta*d/number_samples*(integrand(m)+integrand(m+1));
%     end
%     series_coefficients(i+1) = integral_sum*N/(2*beta*d);
% end
% B = zeros(1,number_samples);
% for i = 0:(N-1)
%     % Computes the array factor for each sample point
%     %A(i) = (1/N)*transpose(series_coefficients)*exp(1j*theta(i)*transpose(linspace(0,N-1,N)));
%     B = B + series_coefficients(i+1)*exp(1j*i*theta);
% end
% B = B/N;

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

fprintf('Error is %d\n',sumsqr(abs(A-Derived_AF)));

plot(theta,abs(Derived_AF),'b');