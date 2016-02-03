function [ output_args ] = plot_fixation( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

S = 0:0.01:0.5
for i = 1:length(S)
s = S(i);
alpha = 5; %partial immune population
%alpha = 10; %neutral population
%alpha = 1000; %naive population
%mu = 2*(10^(-9));
mu = 0.2*(10^(-9));
N = 3.3*(10^(7));
lambda(i) = expected_b_mutations( s, alpha, mu, N );
Pr_fix(i) = pi_survive(s)*exp(-lambda(i));
end

end

