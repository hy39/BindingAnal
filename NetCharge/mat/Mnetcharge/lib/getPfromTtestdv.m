function [slope P se T r2] = getPfromTtestdv(X,y,dv)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
beta = ((X'*X)^-1)*(X'*y);
slope = beta(length(beta(:,1)));
x = X(:,end);
n = length(x);
df = n - length(X(1,:));
r = y-(dv*beta(1:end-1)+beta(end).*x);
ssr = sum(r.^2); %the residual sum of squares 
xi = x-mean(x);
se = sqrt(ssr/((df).*(xi'*xi)));
T = beta(end)/se;
%P = 1-tcdf(T,n-1);
P = 1-tcdf(T,df);
sst = sum((y-mean(y)).^2); %the total sum of squares 
r2 = 1-ssr./sst;
end

