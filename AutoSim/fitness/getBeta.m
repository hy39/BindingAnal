function [ B_Trans ] = getBeta( sk,delta,V )
%Subplot1 f(k,V)
%Subplot2 g(V) = exp(-a*V.^b)
%Subplot3 Beta(k,V)
%Subplot4 dBeta(k,V)/dV
%Author: Hsiang-Yu Yuan
%1st version: Jul 16, 2012
%Rho = f x g
%2nd verison: Aug 01, 2013
%Rho = 1 - R0^-1
%R0 = f x g x n

%Binding avidity range
%V = 0.8;

%Transmission parameters
%p = 2;
%r = 1;
%r = 2;
%b = 3;
%a = 0.7;
p = 4;
r = 70;
b = 3;
a = 0.7;
c = 0.7; % contact rate


c = 0.5; % contact rate
nv = 4; % average copies number of each virion
   k = (1:length(sk)) - 1;
   P_Ab = exp(-p*(V+1));
   P_Trans = (1-P_Ab).^(r*(k-delta)); 
   
   P_Rep = exp(-a*V.^b);

   R0_Trans = P_Trans.*P_Rep.*nv;
   Rho_Trans = 1 - R0_Trans.^-1;
   Rho_Trans(find(Rho_Trans<0))=0;
   B_Trans = c.*Rho_Trans;
end

   

