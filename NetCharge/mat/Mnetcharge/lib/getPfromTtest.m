function [slope P se T] = getPfromTtest(X,y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
beta = ((X'*X)^-1)*(X'*y);
slope = beta(length(beta(:,1)));
x = X(:,2);
n = length(x);
r = y-(beta(1)+beta(2).*x);
ssr = sum(r.^2);
xi = x-mean(x);
se = sqrt(ssr/((n-2).*(xi'*xi)));
T = beta(2)/se;
%if T>= 0 
    P = 1-tcdf(T,n-2);
%else
%    P = tcdf(T,n-2);
%end

end

