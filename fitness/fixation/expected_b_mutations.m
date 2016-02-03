function [ exp_b ] = expected_b_mutations( s, alpha, mu, N )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

  exp_b = mu/s*(N)*log(N)*exp(-alpha*s)* pi_survive(s+1/alpha);
    
end

