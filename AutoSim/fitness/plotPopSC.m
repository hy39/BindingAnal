% plot the selection coefficient with mutational effects under different Sk
% scenario
function [  ] = plotPopSC( )


X1 = normpdf(1:12,1);
X2 = normpdf(1:12,4);
X3 = normpdf(1:12,8);
delta = 1:15;
X = X1;
X0 = 1-sum(X);
X(1) = X(1) + X0;
sk = X;
s = getPopSC( delta, sk)
figure;
plot(s);
hold on;

X = X2;
X0 = 1-sum(X);
X(1) = X(1) + X0;
sk = X;
s = getPopSC( delta, sk)
plot(s);

X = X3;
X0 = 1-sum(X);
X(1) = X(1) + X0;
sk = X;
s = getPopSC( delta, sk)
plot(s);
end