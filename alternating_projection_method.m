

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

for i = 0:(N-1)
%     series_coefficients(i) = (1/number_samples)*A*exp(-1j*(i-1)*transpose(theta));
%     series_coefficients(i) = series_coefficients(i)/abs(series_coefficients(i));
    
    integrand = A.*exp(-1j*i*theta);
    integral_sum = 0;
    for m = 1:number_samples-1
        integral_sum = integral_sum + beta*d/number_samples*(integrand(m)+integrand(m+1));
    end
    series_coefficients(i+1) = integral_sum*N/(2*beta*d);
end

B = zeros(1,number_samples);
for i = 1:number_samples
    % Computes the array factor for each sample point
    B(i) = (1/N)*transpose(series_coefficients)*exp(1j*theta(i)*transpose(linspace(0,N-1,N)));
end

fprintf('Error is %d\n',sumsqr(abs(A-B)));

plot(theta,abs(B),'b');