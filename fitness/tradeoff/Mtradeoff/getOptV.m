function [ R ] = getOptV( v, delta, sk )
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

%Sk
X = normpdf(1:20,8);
%X = normpdf(1:20,4);
%X = normpdf(1:20,1);
X0 = 1-sum(X);
X(1) = X(1) + X0;
if exist('sk')
else
    sk = X;
end
%Xnorm = X/sum(X);
%sk = Xnorm;
%Sk = [0.0000    0.0000    0.0000    0.0001    0.0044     0.0540    0.2420    0.3989    0.2420    0.0540  0.0044    0.0001    0.0000    0.0000    0.0000  0.0000    0.0000    0.0000    0.0000    0.0000]

%Binding avidity range
if exist('v')
V = v
else
V = repmat(0.8,1,length(sk));
end

%Transmission parameters
p = 4;
%r = 1;
r = 70;
b = 3;
a = 0.7;
c = 0.5; % contact rate
nv = 4; % average copies number of each virion
gamma = 1/5;
sk = sk./sum(sk); %normalize
for i=1:length(sk)
   k = i-1;
   j = k - delta;
   j(find(j < 0)) = 0; 
   P_Ab = exp(-p*(V+1));
   P_Trans = (1-P_Ab).^(r*(j)); 
   P_Rep = exp(-a*V.^b);

   R0_Trans = P_Trans.*P_Rep.*nv;
   Rho_Trans = 1 - R0_Trans.^-1;
   Rho_Trans(find(Rho_Trans<0))=0
   B_Trans(i,:) = c.*Rho_Trans; 
end
   rsk = repmat(sk',1,length(delta));
   B_Trans = sum(B_Trans.*rsk)
   R = B_Trans*c/gamma;
end
   

